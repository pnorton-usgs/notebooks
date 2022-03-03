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
#     display_name: Python [conda env:bandit_py3]
#     language: python
#     name: conda-env-bandit_py3-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import pandas as pd
import numpy as np
import pydot
import networkx as nx

from pyPRMS.ParameterFile import ParameterFile

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190916_Pipestem_ST'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20180809_pipestem'

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190208_DelawareRiver'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190916_Delaware_CONUS'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190913_Delware_streamtemp'
filename = '{}/myparam.param'.format(workdir)
plot_filename = '{}/Delaware_base.png'.format(workdir)

# Load parameter file
pfile = ParameterFile(filename)

# %%
# cparam = 'hru_up_id'
# var_df = pfile.parameters.get_dataframe(cparam)

# %%
# var_df.head()

# %%
# Get hru_segment parameter
tosegment = pfile.parameters['tosegment'].data
hru_segment = pfile.parameters['hru_segment'].data
hru_segment_nhm = pfile.parameters['hru_segment_nhm'].data
nhm_id = pfile.parameters['nhm_id'].data
nhm_seg = pfile.parameters['nhm_seg'].data
tosegment_nhm = pfile.parameters['tosegment_nhm'].data

# Build the stream network
dag_streamnet = nx.DiGraph()
for ii, vv in enumerate(tosegment):
    if vv == 0:
        dag_streamnet.add_edge(nhm_seg[ii], 'Out_{}'.format(ii+1))
#         dag_streamnet.add_edge(ii + 1, 'Out_{}'.format(ii + 1))
    else:
        dag_streamnet.add_edge(nhm_seg[ii], tosegment_nhm[ii])
#         dag_streamnet.add_edge(ii + 1, vv)

# Now add the HRUs
nr_cnt = 0

for ii, vv in enumerate(hru_segment_nhm):
    hru_node = 'H_{}'.format(nhm_id[ii])
#     if hru_segment_nhm[ii] in [4184, 2686, 2696, 2691]:
#         print(hru_node)
#         dag_streamnet.add_node(hru_node, style='filled', fillcolor='purple')
#     else:
    dag_streamnet.add_node(hru_node, style='filled', fillcolor='yellow')
    
    if vv == 0:
        nr_cnt += 1
        strm_node = 'NR_{}'.format(nr_cnt)
        
        dag_streamnet.add_node(strm_node, fontcolor='white', style='filled', fillcolor='red')
        dag_streamnet.add_edge(hru_node, strm_node)
    else:
        dag_streamnet.add_edge(hru_node, vv)
        

# %%
    
print("Number of subset nodes: {}".format(dag_streamnet.number_of_nodes()))
print("Number of subset edges: {}".format(dag_streamnet.number_of_edges()))
print("Is a DAG: {}".format(nx.is_directed_acyclic_graph(dag_streamnet)))


# %%

# %%
# Create plot of the graph using pydot
pd = nx.nx_pydot.to_pydot(dag_streamnet)
pd.write_png(plot_filename)

# %%

# %% [markdown]
# ### Test subsetting the cascade network

# %%
# Create upstream graph
dag_rev = dag_uphru.reverse()

dsmost_nodes = [3547]

# Get all unique segments u/s of the starting segment
uniq_seg = set()
    
for xx in dsmost_nodes:
    try:
        pred = nx.dfs_predecessors(dag_rev, xx)
        uniq_seg = uniq_seg.union(set(pred.keys()).union(set(pred.values())))
    except KeyError:
        print('\nKeyError: Segment {} does not exist in stream network'.format(xx))

# Get a subgraph in the dag_ds graph and return the edges
dag_uphru_subset = dag_uphru.subgraph(uniq_seg).copy()

# 2018-02-13 PAN: It is possible to have outlets specified which are not truly
#                 outlets in the most conservative sense (e.g. a point where
#                 the stream network exits the study area). This occurs when
#                 doing headwater extractions where all segments for a headwater
#                 are specified in the configuration file. Instead of creating
#                 output edges for all specified 'outlets' the set difference
#                 between the specified outlets and nodes in the graph subset
#                 which have no edges is performed first to reduce the number of
#                 outlets to the 'true' outlets of the system.
# node_outlets = [ee[0] for ee in dag_ds_subset.edges()]
# true_outlets = set(dsmost_seg).difference(set(node_outlets))
# print('node_outlets: {}'.format(','.join(map(str, node_outlets))))
# print('true_outlets: {}'.format(','.join(map(str, true_outlets))))

# Add the downstream segments that exit the subgraph
# for xx in true_outlets:
#     dag_ds_subset.add_edge(xx, 'Out_{}'.format(xx))
    
print("Number of subset nodes: {}".format(dag_uphru_subset.number_of_nodes()))
print("Number of subset edges: {}".format(dag_uphru_subset.number_of_edges()))

# %%
# Create plot of the graph using pydot
pd = nx.nx_pydot.to_pydot(dag_uphru_subset)
pd.write_png('{}/sagehen_subset.png'.format(workdir))

# %%
