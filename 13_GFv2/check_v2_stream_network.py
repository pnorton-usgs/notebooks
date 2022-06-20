# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python [conda env:bandit_38]
#     language: python
#     name: conda-env-bandit_38-py
# ---

# %%

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

    seg_data = []
    toseg_data = []
    
    # Skip header
    next(it)

    # Read the parameter values
    for rec in it:
        crec = rec.split(',')
        
        seg_data.append(int(crec[0]))
        try:
            toseg_data.append(int(crec[1]))
        except ValueError:
            toseg_data.append(0)

    return seg_data, toseg_data


# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/FY2021_projects/modeling_fabric'

segid, tosegment = read_csv(f'{workdir}/GFv2_stream_network.csv')

# tosegment_nhm = read_csv('{}/tosegment_nhm.csv'.format(workdir))
# hru_segment_nhm = read_csv('{}/hru_segment_nhm.csv'.format(workdir))
# seg_id_nhm = read_csv('{}/seg_id_nhm.csv'.format(workdir))
# hru_id_nhm = read_csv('{}/hru_id_nhm.csv'.format(workdir))

# %%

# %%
# Build the stream network
include_hrus = False

dag_streamnet = nx.DiGraph()
for ii, vv in enumerate(tosegment):
    if vv == 0:
        dag_streamnet.add_edge(segid[ii], 'Out_{}'.format(ii+1))
    else:
        dag_streamnet.add_edge(segid[ii], tosegment[ii])

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
    print(f'{yy},1')     

# %%
print("Number of subset nodes: {}".format(dag_streamnet.number_of_nodes()))
print("Number of subset edges: {}".format(dag_streamnet.number_of_edges()))
print("Is a DAG: {}".format(nx.is_directed_acyclic_graph(dag_streamnet)))

# %%
# Create upstream graph
dag_us = dag_streamnet.reverse()

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
    plot.write_pdf(f'{workdir}/seg_CYCLE_{dsmost_seg[0]}.pdf')

# %%

# %%
