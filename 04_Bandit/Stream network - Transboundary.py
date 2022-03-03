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

# %%
import pandas as pd
import numpy as np
import networkx as nx
import pydot

import matplotlib.pyplot as plt


# %%
def read_csv(filename):
    fhdl = open(filename)
    rawdata = fhdl.read().splitlines()
    fhdl.close()
    it = iter(rawdata)

    data = []

    # Skip header
    next(it)

    # Read the parameter values
    for rec in it:
        data.append(int(rec.split(',')[1]))
    return data


# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/20190927_transboundary'

tosegment_nhm = read_csv('{}/tosegment_nhm.csv'.format(workdir))
hru_segment_nhm = read_csv('{}/hru_segment_nhm.csv'.format(workdir))
seg_id_nhm = read_csv('{}/seg_id_nhm.csv'.format(workdir))
hru_id_nhm = read_csv('{}/hru_id_nhm.csv'.format(workdir))

# %%
# # Output a new nhm_seg file
# fdl = open('{}/nhm_seg.csv'.format(workdir), 'w')

# fdl.write('$id,nhm_seg\n')

# for xx in nhm_seg_out:
#     fdl.write('{0},{0}\n'.format(xx))
# fdl.close()

# %%
# Build the stream network
include_hrus = False

dag_streamnet = nx.DiGraph()
for ii, vv in enumerate(tosegment_nhm):
    if vv == 0:
        dag_streamnet.add_edge(seg_id_nhm[ii], 'Out_{}'.format(ii+1))
    else:
        dag_streamnet.add_edge(seg_id_nhm[ii], tosegment_nhm[ii])

    if include_hrus:
        # Now add the HRUs
        nr_cnt = 0

        for ii, vv in enumerate(hru_segment_nhm):
            hru_node = 'H_{}'.format(hru_id_nhm[ii])
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

# Output any cycles/loops
# Also update attributes for nodes which are part of a cycle
for xx in nx.simple_cycles(dag_streamnet):
    for yy in xx:
        dag_streamnet.nodes[yy]['style'] = 'filled'
        dag_streamnet.nodes[yy]['fillcolor'] = 'orange'
    print(xx)        

# %%
print("Number of subset nodes: {}".format(dag_streamnet.number_of_nodes()))
print("Number of subset edges: {}".format(dag_streamnet.number_of_edges()))
print("Is a DAG: {}".format(nx.is_directed_acyclic_graph(dag_streamnet)))

# %%
# Create upstream graph
dag_us = dag_streamnet.reverse()

# %%
dsmost_seg = [58242]

# Get all unique segments u/s of the starting segment
uniq_seg_us = set()
    
for xx in dsmost_seg:
    try:
        pred = nx.dfs_predecessors(dag_us, xx)
        uniq_seg_us = uniq_seg_us.union(set(pred.keys()).union(set(pred.values())))
    except KeyError:
        print('\nKeyError: Segment {} does not exist in stream network'.format(xx))

# Get a subgraph in the dag_ds graph and return the edges
dag_ds_subset = dag_streamnet.subgraph(uniq_seg_us).copy()

node_outlets = [ee[0] for ee in dag_ds_subset.edges()]
true_outlets = set(dsmost_seg).difference(set(node_outlets))
# print('node_outlets: {}'.format(','.join(map(str, node_outlets))))
# print('true_outlets: {}'.format(','.join(map(str, true_outlets))))

# Add the downstream segments that exit the subgraph
for xx in true_outlets:
    dag_ds_subset.add_edge(xx, 'Out_{}'.format(xx))
    
print("Number of subset nodes: {}".format(dag_ds_subset.number_of_nodes()))
print("Number of subset edges: {}".format(dag_ds_subset.number_of_edges()))

# %%
# Create plot of the subset
plot = nx.nx_pydot.to_pydot(dag_ds_subset)
plot.write_pdf('{}/tb_seg_{}.pdf'.format(workdir, dsmost_seg[0]))

# %%
# Plot entire stream network
plot = nx.nx_pydot.to_pydot(dag_streamnet)
plot.write_pdf('{}/tb_streamnet.pdf'.format(workdir))

# %% [markdown]
# ### Subset Cycle/Loop problems

# %%
# Output any cycles/loops
for xx in nx.simple_cycles(dag_streamnet):
    print('='*40)
    print(xx)
        
    dsmost_seg = xx

    # Get all unique segments u/s of the starting segment
    uniq_seg_us = set()

    for xx in dsmost_seg:
        try:
            pred = nx.dfs_predecessors(dag_us, xx)
            uniq_seg_us = uniq_seg_us.union(set(pred.keys()).union(set(pred.values())))
        except KeyError:
            print('\nKeyError: Segment {} does not exist in stream network'.format(xx))

    # Get a subgraph in the dag_ds graph and return the edges
    dag_ds_subset = dag_streamnet.subgraph(uniq_seg_us).copy()

    node_outlets = [ee[0] for ee in dag_ds_subset.edges()]
    true_outlets = set(dsmost_seg).difference(set(node_outlets))
    # print('node_outlets: {}'.format(','.join(map(str, node_outlets))))
    # print('true_outlets: {}'.format(','.join(map(str, true_outlets))))

    # Add the downstream segments that exit the subgraph
    for xx in true_outlets:
        dag_ds_subset.add_edge(xx, 'Out_{}'.format(xx))

    print("Number of subset nodes: {}".format(dag_ds_subset.number_of_nodes()))
    print("Number of subset edges: {}".format(dag_ds_subset.number_of_edges()))
    
    # Create plot of the subset
    plot = nx.nx_pydot.to_pydot(dag_ds_subset)
    plot.write_pdf('{}/tb_seg_CYCLE_{}.pdf'.format(workdir, dsmost_seg[0]))

# %%
