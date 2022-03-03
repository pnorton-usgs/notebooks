# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python [conda env:bandit]
#     language: python
#     name: conda-env-bandit-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import numpy as np
from collections import OrderedDict
from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile

import pandas as pd
import pydot
import networkx as nx

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_CONUS'
workdir = '/Users/pnorton/tmp/check_paramdb_v11'
pdb = ParamDb(paramdb_dir=workdir, verbose=True, verify=True)
params = pdb.parameters

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/pipestem'
# filename = '{}/myparam.param'.format(workdir)

# # Load parameter file
# pfile = ParameterFile(filename, verbose=True, verify=True)
# params = pfile.parameters

# %%
hru_segment_nhm = params['hru_segment_nhm'].tolist()
hru_segment = params['hru_segment'].tolist()
hru_area = params['hru_area'].tolist()

nhm_id = params['nhm_id'].tolist()
nhm_seg = params['nhm_seg'].tolist()
tosegment_nhm = params['tosegment_nhm'].tolist()

# %% [markdown]
# ## Compute the HRU area for each segment

# %%
# Sum up the HRU area for each segment
area_by_seg = OrderedDict()

for seg, harea in zip(hru_segment_nhm, hru_area):
    if seg in area_by_seg:
        area_by_seg[seg] += harea
    else:
        area_by_seg[seg] = harea

# %% [markdown]
# ## Build the stream network

# %%
num_outlets = 0
dag_streamnet = nx.DiGraph()

for ii, vv in enumerate(tosegment_nhm):
    if vv == 0 or vv not in nhm_seg:
#     if vv == 0 or tosegment[ii] == 0:
        dag_streamnet.add_node(nhm_seg[ii], style='filled', fontcolor='white', fillcolor='blue')
        dag_streamnet.add_node(vv, style='filled', fontcolor='white', fillcolor='grey')
        dag_streamnet.add_edge(nhm_seg[ii], vv)

        num_outlets += 1
    else:
        dag_streamnet.add_edge(nhm_seg[ii], tosegment_nhm[ii])
        if nhm_seg[ii] > 56460:
            # NOTE: This will only work correctly prior to renumbering the NHM segments
            dag_streamnet.nodes[nhm_seg[ii]]['fillcolor'] = 'deeppink'
            dag_streamnet.nodes[nhm_seg[ii]]['style'] = 'filled'

# Output any cycles/loops
# Also update attributes for nodes which are part of a cycle
for xx in nx.simple_cycles(dag_streamnet):
    for yy in xx:
        dag_streamnet.nodes[yy]['style'] = 'filled'
        dag_streamnet.nodes[yy]['fillcolor'] = 'orange'
    print(xx)
    
print("Number of subset nodes: {}".format(dag_streamnet.number_of_nodes()))
print("Number of subset edges: {}".format(dag_streamnet.number_of_edges()))
print("Number of outlet nodes: {}".format(num_outlets))
print("Is a DAG: {}".format(nx.is_directed_acyclic_graph(dag_streamnet)))    

print("Isolate nodes:")
print(list(nx.isolates(dag_streamnet)))

# %% [markdown]
# ## Compute cumulative drainage area by stream segment

# %%
fhdl = open(f'{workdir}/seg_cum_area.csv', 'w')
fhdl.write('$id,seg_cum_area\n')

# Create upstream graph
dag_us = dag_streamnet.reverse()

# dsmost_segs = list(dag_streamnet.nodes)
dsmost_segs = nhm_seg
print(f'Number of nodes: {len(dsmost_segs)}')

for dsmost_seg in dsmost_segs:
    uniq_seg_us = set()
    
    try:
        pred = nx.dfs_predecessors(dag_us, dsmost_seg)
        uniq_seg_us = uniq_seg_us.union(set(pred.keys()).union(set(pred.values())))
    except KeyError:
        print('\nKeyError: Segment {} does not exist in stream network'.format(dsmost_seg))

    dag_ds_subset = dag_streamnet.subgraph(uniq_seg_us).copy()

    node_outlets = [ee[0] for ee in dag_ds_subset.edges()]
#     true_outlets = []
    true_outlets = set([dsmost_seg]).difference(set(node_outlets))

    # Add the downstream segments that exit the subgraph
    for xx in true_outlets:
        nhm_outlet = list(dag_streamnet.neighbors(xx))[0]
        dag_ds_subset.add_node(nhm_outlet, style='filled', fontcolor='white', fillcolor='grey')
        dag_ds_subset.add_edge(xx, nhm_outlet)

    toseg_idx = list(set(xx[0] for xx in dag_ds_subset.edges))
    
    csum = 0.0
    for xx in toseg_idx:
        csum += area_by_seg.get(xx, 0.0)
    fhdl.write(f'{dsmost_seg},{csum}\n')

fhdl.close()

# %%

# %%

# %%

# %%
