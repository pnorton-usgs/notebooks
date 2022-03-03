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
import networkx as nx
import pydot
from collections import Counter, OrderedDict

import matplotlib.pyplot as plt

from pyPRMS.ParamDb import ParamDb


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
def check_for_disconnected_graphs(g):
    # Check for disconnected graphs within the graph object
    tot_weak_comp = 0
    node_size_cnt = Counter()

    for xx in nx.weakly_connected_components(g):
        out_node = [ii for ii in list(xx) if 'Out_' in str(ii)]
        tot_weak_comp += 1
        node_size_cnt[len(xx)] += 1

        try:
            print('{}: nodes {}'.format(out_node[0], len(xx)))
        except IndexError:
            # If this occurs it most likely indicates a loop or cycle
            print('WEIRD: ', xx)
    print('Total weakly connected components: {}'.format(tot_weak_comp))

    print('Top 10 node sizes (size, count)')
    print(node_size_cnt.most_common(10))


# %%
def output_cycles(g, out_dir):
    # Output any cycles/loops as a subset using pydot
    
    # Create upstream graph
    dag_us = g.reverse()
    
    for xx in nx.simple_cycles(g):
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
        plot.write_pdf('{}/seg_CYCLE_{}.pdf'.format(out_dir, dsmost_seg[0]))


# %%
workdir_merged = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/20190927_transboundary/paramdb_merged'
workdir_tb = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/20190927_transboundary'

mrg_hru_segment_nhm = read_csv('{}/hru_segment_nhm.csv'.format(workdir_merged))
mrg_tosegment_nhm = read_csv('{}/tosegment_nhm.csv'.format(workdir_merged))
mrg_nhm_id = read_csv('{}/nhm_id.csv'.format(workdir_merged))
mrg_nhm_seg = read_csv('{}/nhm_seg.csv'.format(workdir_merged))


tb_cutoff_segs = read_csv('{}/seg_id_nat_v1.csv'.format(workdir_tb))
# tb_tosegment_nhm = read_csv('{}/tosegment_nhm.csv'.format(workdir_tb))
# tb_hru_segment_nhm = read_csv('{}/hru_segment_nhm.csv'.format(workdir_tb))
# tb_seg_id_nhm = read_csv('{}/seg_id_nhm.csv'.format(workdir_tb))
# tb_hru_id_nhm = read_csv('{}/hru_id_nhm.csv'.format(workdir_tb))
tb_v1repl_nhm_id = read_csv('{}/hru_v1repl_nhm_id.csv'.format(workdir_tb))
tb_type_nca = read_csv('{}/type_nca.csv'.format(workdir_tb))

# new_nhm_seg_map = OrderedDict((val, idx) for idx, val in enumerate(tb_seg_id_nhm))

# %% [markdown]
# ### Create NHM stream network and remove replaced segments

# %%
# %time pdb = ParamDb(paramdb_dir=workdir_merged, verbose=True, verify=True)

# Must use tolist() to have fastest access to these arrays
tosegment_nhm = pdb.parameters['tosegment_nhm'].data.tolist()
nhm_seg = pdb.parameters['nhm_seg'].data.tolist()

orig_nhm_seg_map = pdb.parameters['nhm_seg'].index_map

# %%
# # %%time
# Build the stream network
num_outlets = 0
include_hrus = False

dag_streamnet = nx.DiGraph()

for ii, vv in enumerate(tosegment_nhm):
    if vv == 0:
        dag_streamnet.add_edge(nhm_seg[ii], 'Out_{}'.format(nhm_seg[ii]))
        num_outlets += 1
    else:
        dag_streamnet.add_edge(nhm_seg[ii], tosegment_nhm[ii])

    if nhm_seg[ii] in tb_cutoff_segs:
        dag_streamnet.nodes[nhm_seg[ii]]['style'] = 'filled'
        dag_streamnet.nodes[nhm_seg[ii]]['fillcolor'] = 'red'
        
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
    
print("Number of subset nodes: {}".format(dag_streamnet.number_of_nodes()))
print("Number of subset edges: {}".format(dag_streamnet.number_of_edges()))
print("Number of outlet nodes: {}".format(num_outlets))
print("Is a DAG: {}".format(nx.is_directed_acyclic_graph(dag_streamnet)))    

print("Isolate nodes:")
print(list(nx.isolates(dag_streamnet)))

# %%
check_for_disconnected_graphs(dag_streamnet)

# %% [markdown]
# ### Connect the TB stuff to the old NHM

# %%
tb_in = {30479: 58873, 
         48360: 61096,
         440: 56613,
         899: 56613,
         48417: 61096}

# Connect the new segments flowing from NHM to TB
for nhm_node, tb_node in tb_in.items():
    dag_streamnet.nodes[nhm_node]['style'] = 'filled'
    dag_streamnet.nodes[nhm_node]['fillcolor'] = 'green'
    dag_streamnet.add_edge(nhm_node, tb_node)

# %%
check_for_disconnected_graphs(dag_streamnet)

# %%
# Create upstream graph
dag_us = dag_streamnet.reverse()

# %%
# Remove the segments that will be replaced by the transboundary segments
for xx in tb_cutoff_segs:
    try:
        dag_us.remove_node(xx)
    except KeyError:
        print('WARNING: nhm_segment {} does not exist in stream network'.format(xx))
        
# Update the downstream graph
dag_streamnet = dag_us.reverse()

print("Number of subset nodes: {}".format(dag_us.number_of_nodes()))
print("Number of subset edges: {}".format(dag_us.number_of_edges()))
print("Is a DAG: {}".format(nx.is_directed_acyclic_graph(dag_us)))    

# %%
dag_ds_subset = dag_us.reverse()

# Create list of toseg ids for the model subset
# NOTE: These are really nhm_seg ids whereas xx[1] would be tosegment_nhm
toseg_idx = list(set(xx[0] for xx in dag_ds_subset.edges))
toseg_idx0 = [xx-1 for xx in toseg_idx]  # 0-based version of toseg_idx

# %% [markdown]
# ## Example subset

# %%
dsmost_seg = [33298]

uniq_seg_us = set()
if dsmost_seg:
    for xx in dsmost_seg:
        try:
            pred = nx.dfs_predecessors(dag_us, xx)
            uniq_seg_us = uniq_seg_us.union(set(pred.keys()).union(set(pred.values())))
        except KeyError:
            print('\nKeyError: Segment {} does not exist in stream network'.format(xx))

dag_ds_subset = dag_streamnet.subgraph(uniq_seg_us).copy()

node_outlets = [ee[0] for ee in dag_ds_subset.edges()]
true_outlets = set(dsmost_seg).difference(set(node_outlets))
# print('node_outlets: {}'.format(','.join(map(str, node_outlets))))
# print('true_outlets: {}'.format(','.join(map(str, true_outlets))))

# Add the downstream segments that exit the subgraph
for xx in true_outlets:
    nhm_outlet = list(dag_streamnet.neighbors(xx))[0]
    dag_ds_subset.nodes[xx]['style'] = 'filled'
    dag_ds_subset.nodes[xx]['fontcolor'] = 'white'
    dag_ds_subset.nodes[xx]['fillcolor'] = 'blue'
    dag_ds_subset.add_node(nhm_outlet, style='filled', fontcolor='white', fillcolor='grey')
    dag_ds_subset.add_edge(xx, nhm_outlet)

toseg_idx = list(set(xx[0] for xx in dag_ds_subset.edges))
toseg_idx0 = [xx-1 for xx in toseg_idx]  # 0-based version of toseg_idx
# print(toseg_idx)
# print(toseg_idx0)

plot = nx.nx_pydot.to_pydot(dag_ds_subset)
plot.write_pdf('{}/mrg_seg_{}.pdf'.format(workdir_merged, dsmost_seg[0]))

# %%
# list(dag_streamnet.neighbors(30119))[0]

# %%

# %%
dag_ds_subset = dag_us.reverse()

# Create list of toseg ids for the model subset
# NOTE: These are really nhm_seg ids whereas xx[1] would be tosegment_nhm
toseg_idx = list(set(xx[0] for xx in dag_ds_subset.edges))
toseg_idx0 = [xx-1 for xx in toseg_idx]  # 0-based version of toseg_idx

# %%
len(toseg_idx)

# %%
# NOTE: what happens if a true_outlet is actually zero (e.g. coastal segment)?
# for xx in true_outlets:
#     dag_ds_subset.add_edge(xx, list(dag_streamnet.neighbors(30119))[0])

# %%
list(dag_streamnet.neighbors(33298))

# %%
dag_ds_subset.edges

# %%
# Create a mapping of nhm_seg (ee[0]) to tosegment_nhm (ee[1])
aa = OrderedDict((ee[0], ee[1]) for ee in dag_ds_subset.edges)
# print(aa)

# Use the mapping to create subsets of nhm_seg, tosegment_nhm, and tosegment
new_nhm_seg = [kk for kk in aa.keys()]
# print('nhm_seg: {}'.format(','.join(map(str, new_nhm_seg))))

new_tosegment_nhm = [kk for kk in aa.values()]
# print('tosegment_nhm: {}'.format(','.join(map(str, new_tosegment_nhm))))

# Using a dictionary for speed
new_nhm_seg_dict = {}
for ii, ss in enumerate(new_nhm_seg):
    new_nhm_seg_dict[ss] = ii+1

# Generate the renumbered NHM tosegments
new_tosegment = [new_nhm_seg_dict[kk] if kk in new_nhm_seg_dict else 0 for kk in aa.values()]
# new_tosegment = [new_nhm_seg.index(kk)+1 if kk in new_nhm_seg else 0 for kk in aa.values()]
# print('tosegment: {}'.format(','.join(map(str, new_tosegment))))


# %%
# Generate the renumbered NHM segment ids
nhm_seg_to_renum = [new_nhm_seg_dict[kk] for kk in new_nhm_seg]
print(nhm_seg_to_renum)


# %%
# Basic idea for using _map from original and new datasets to build new arrays for the merged dataset.

crap = []

for nseg in aa.keys():
    if nseg > 56460:
        pass
        crap.append(new_nhm_seg_map[nseg])
    else:
        # Lookup index from old
        crap.append(orig_nhm_seg_map[nseg])
print(crap)

# %% [markdown]
# ### Check for Cycle/Loop problems

# %%
output_cycles(dag_streamnet, workdir_nhm)

# %% [markdown]
# # Merge NHM and TB HRUs

# %%
# paramdb_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_tb_work'

# %%
# print('Loading NHM ParamDb')
# pdb = ParamDb(paramdb_dir)

# %%
nhm_params = pdb.parameters
nhm_global_dimensions = pdb.dimensions

hru_segment_nhm = nhm_params.get('hru_segment_nhm').data.tolist()
nhm_id = nhm_params.get('nhm_id').data.tolist()

print('Number of nhm_id: {}'.format(len(nhm_id)))
print('Number of hru_segment_nhm: {}'.format(len(hru_segment_nhm)))

# %%
# Remove the HRUs being replaced by the TB HRUs
# NOTE: this is a very slow way to do this -- but it works.
addl_removals = [46760, 65421, 99973, 100528, 102300]
nhm_id_new = []
hru_segment_nhm_new = []

removed_hru_idx = []

for idx, xx in enumerate(nhm_id):
    if xx not in tb_v1repl_nhm_id:
        if xx not in addl_removals:
            nhm_id_new.append(xx)
            hru_segment_nhm_new.append(hru_segment_nhm[idx])
        else:
            removed_hru_idx.append(idx)
    else:
        removed_hru_idx.append(idx)
        
# for xx in tb_v1repl_nhm_id:
#     idx = nhm_id.index(xx)
#     nhm_id.remove(xx)
#     hru_segment_nhm.remove(hru_segment_nhm[idx])

print('Number of nhm_id: {}'.format(len(nhm_id_new)))
print('Number of hru_segment_nhm: {}'.format(len(hru_segment_nhm_new))) 

# %%
# Use removed_hru_idx to remove entries from other hru-based parameters
removed_hru_idx

# %%
# Create dictionary mapping of nhm_id to index position and hru_segment_nhm
nhm_id_to_hru_segment = {}

for idx, xx in enumerate(nhm_id_new):
    nhm_id_to_hru_segment[xx] = [idx, hru_segment_nhm_new[idx]]
    
print(nhm_id_to_hru_segment)

# %%
len(nhm_id_to_hru_segment)

# %%
nhm_id_to_idx = {}
for ii, vv in enumerate(nhm_id_new):
    # keys are 1-based, values are 0-based
    nhm_id_to_idx[vv] = ii

# print('Number of NHM hru_segment entries: {}'.format(len(hru_segment)))

# Create a dictionary mapping segments to HRUs
seg_to_hru = {}
for ii, vv in enumerate(hru_segment_nhm_new):
    # keys are 1-based, values in arrays are 1-based
    seg_to_hru.setdefault(vv, []).append(ii + 1)

# Get HRU ids ordered by the segments in the model subset - entries are 1-based
hru_order_subset = []
for xx in toseg_idx:
    if xx in seg_to_hru:
        for yy in seg_to_hru[xx]:
            hru_order_subset.append(yy)
    else:
        print('Stream segment {} has no HRUs connected to it.'.format(xx))
        # raise ValueError('Stream segment has no HRUs connected to it.')


# %%
len(hru_order_subset)

# %%
cnt = 0

for idx, xx in enumerate(hru_segment_nhm_new):
    if xx == 0:
        cnt += 1
    elif xx not in toseg_idx:
        print('{}: {}'.format(idx, xx))
print('Count of missing stream segments: {}'.format(cnt))

# %%
print('Number of hru_segment_nhm: {}'.format(len(hru_segment_nhm_new)))  

# %%
nhm_id[1236]

# %%
hru_segment_nhm[1]

# %%
nhm_id[1]

# %%
nhm_id_to_hru_segment[1236]

# %%
nhm_id.index(1236)

# %%
tb_v1repl_nhm_id[7144]

# %%
nhm_id_to_idx[109950]

# %%
toseg_idx

# %%
