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
#     display_name: Python [conda env:bandit_38]
#     language: python
#     name: conda-env-bandit_38-py
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
# workdir = '/Users/pnorton/tmp'
# # filename = '{}/sagehen_gsflow.params'.format(workdir)
# # plot_filename = '{}/sagehen_gsflow.png'.format(workdir)
# filename = '{}/sagehen_grid.param'.format(workdir)
# plot_filename = '{}/sagehen_grid.png'.format(workdir)

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/tmp/nvwsc'
# filename = '{}/Desch.params6'.format(workdir)
# plot_filename = '{}/Desch.png'.format(workdir)

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/gsflow_examples/sagehen_2lay_nolakes/input/prms'
# filename = f'{workdir}/sagehen_grid_test.params'
# plot_filename = '/Users/pnorton/tmp/sagehen_grid_test.png'

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/sagehen/input/prms'
filename = f'{workdir}/prms_orig_v2.params'
plot_filename = f'{workdir}/sagehen_prms_orig_grid.png'

# Load parameter file
pfile = ParameterFile(filename)

# %%
# cparam = 'hru_up_id'
# var_df = pfile.parameters.get_dataframe(cparam)

# %%
# var_df.head()

# %%
# Get cascade-related parameters
hru_up_id = pfile.parameters['hru_up_id'].data
hru_down_id = pfile.parameters['hru_down_id'].data
hru_pct_up = pfile.parameters['hru_pct_up'].data
hru_strmseg_down_id = pfile.parameters['hru_strmseg_down_id'].data

# Build the network
dag_uphru = nx.DiGraph()
for ii, vv in enumerate(hru_up_id):
    # ii => cascade number
    # vv => hru id
    dag_uphru.add_node(vv)
    
    if hru_strmseg_down_id[ii] == 0 and hru_pct_up[ii] > 0.0:
        dag_uphru.add_node(hru_down_id[ii])
        dag_uphru.add_edge(vv, hru_down_id[ii])
    else:
        # HRU drains to a stream segment
        strm_node = 'S_{}'.format(hru_strmseg_down_id[ii])
        dag_uphru.add_node(strm_node, fontcolor='white', style='filled', fillcolor='blue')
        dag_uphru.add_edge(vv, strm_node)

# %%
    
print("Number of subset nodes: {}".format(dag_uphru.number_of_nodes()))
print("Number of subset edges: {}".format(dag_uphru.number_of_edges()))
print("Is a DAG: {}".format(nx.is_directed_acyclic_graph(dag_uphru)))


# %%
# Output any cycles/loops
# Also update attributes for nodes which are part of a cycle
for xx in nx.simple_cycles(dag_uphru):
    for yy in xx:
        dag_uphru.nodes[yy]['style'] = 'filled'
        dag_uphru.nodes[yy]['fillcolor'] = 'orange'
    print(xx)

# %%
# Create plot of the graph using pydot
pd = nx.nx_pydot.to_pydot(dag_uphru)
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
