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
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
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
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190507_leavesley/13139510'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190913_Delware_streamtemp'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190916_Delaware_CONUS'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20200121_pipestem_v11'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20200122_696_v11'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20200326_cape_fear'

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190916_Pipestem_ST'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20180809_pipestem'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20191007_pipestem'

filename = '{}/myparam.param'.format(workdir)
plot_filename = '{}/model_cape_fear.pdf'.format(workdir)

# Load parameter file
pfile = ParameterFile(filename)

# Get parameters
# tosegment = pfile.parameters['tosegment'].tolist()
# hru_segment = pfile.parameters['hru_segment'].tolist()
hru_segment_nhm = pfile.parameters['hru_segment_nhm'].tolist()
nhm_id = pfile.parameters['nhm_id'].tolist()
nhm_seg = pfile.parameters['nhm_seg'].tolist()
tosegment_nhm = pfile.parameters['tosegment_nhm'].tolist()

poi = {}
if pfile.parameters.exists('poi_gage_segment'):
    poi_gage_id = pfile.parameters['poi_gage_id'].tolist()
    poi_gage_segment = pfile.parameters['poi_gage_segment'].tolist()
    
    for gg, ss in zip(poi_gage_id, poi_gage_segment):
        poi[ss] = gg


# %% [markdown]
# ## Generate stream network

# %%
# Build the stream network
num_outlets = 0
include_hrus = False

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
            
                                
    if ii + 1 in poi:
        dag_streamnet.nodes[nhm_seg[ii]]['shape'] = 'box'
        dag_streamnet.nodes[nhm_seg[ii]]['label'] = '{}\n POI: {}'.format(nhm_seg[ii], poi[ii + 1])
        
if include_hrus:
    # Add the HRUs
    nr_cnt = 0

    for ii, vv in enumerate(hru_segment_nhm):
        hru_node = 'H_{}'.format(nhm_id[ii])

        dag_streamnet.add_node(hru_node, style='filled', fillcolor='yellow')

        if vv == 0:
            nr_cnt += 1
            strm_node = 'NR_{}'.format(nr_cnt)

            dag_streamnet.add_node(strm_node, fontcolor='white', style='filled', fillcolor='red')
            dag_streamnet.add_edge(hru_node, strm_node)
        else:
            dag_streamnet.add_edge(hru_node, vv)

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



# node_outlets = [ee[0] for ee in dag_ds_subset.edges()]
# true_outlets = set(dsmost_seg).difference(set(node_outlets))
# print('node_outlets: {}'.format(','.join(map(str, node_outlets))))
# print('true_outlets: {}'.format(','.join(map(str, true_outlets))))

# # Add the downstream segments that exit the subgraph
# for xx in true_outlets:
#     nhm_outlet = list(dag_streamnet.neighbors(xx))[0]
#     dag_ds_subset.nodes[xx]['style'] = 'filled'
#     dag_ds_subset.nodes[xx]['fontcolor'] = 'white'
#     dag_ds_subset.nodes[xx]['fillcolor'] = 'blue'
#     dag_ds_subset.add_node(nhm_outlet, style='filled', fontcolor='white', fillcolor='grey')
#     dag_ds_subset.add_edge(xx, nhm_outlet)

# %% [markdown]
# ### Output stream network with HRUs using pydot

# %%
# Create plot of the graph using pydot
pd = nx.nx_pydot.to_pydot(dag_streamnet)
pd.write_pdf(plot_filename)

# %%

# %%
type(dag_streamnet)

# %%
crap = dag_streamnet.reverse()

# %%
type(crap)

# %%
