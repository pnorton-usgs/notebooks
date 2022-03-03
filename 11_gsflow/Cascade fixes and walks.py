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
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/sagehen/input/prms'
filename = f'{workdir}/prms_orig_v2.params'
plot_filename = f'{workdir}/sagehen_prms_orig_grid.png'

# Load parameter file
pfile = ParameterFile(filename)

# %%
# Get cascade-related parameters
hru_up_id = pfile.parameters['hru_up_id'].data
hru_down_id = pfile.parameters['hru_down_id'].data
hru_pct_up = pfile.parameters['hru_pct_up'].data
hru_strmseg_down_id = pfile.parameters['hru_strmseg_down_id'].data
tosegment = pfile.parameters['tosegment'].data

# Build the network
dag_uphru = nx.DiGraph()

# Generate the stream network
for ii, vv in enumerate(tosegment):
    if vv == 0:
        dag_uphru.add_node(f'S_{ii+1}', style='filled', fontcolor='white', fillcolor='blue')
        dag_uphru.add_node(f'S_{vv}', style='filled', fontcolor='white', fillcolor='grey')
        dag_uphru.add_edge(f'S_{ii+1}', f'S_{vv}')
    else:
        dag_uphru.add_edge(f'S_{ii+1}', f'S_{vv}')

# Add the hru cascades and stream network connections
for ii, vv in enumerate(hru_up_id):
    # ii => cascade number
    # vv => hru id
    dag_uphru.add_node(vv)
    
#     if hru_pct_up[ii] > 0.0:
    if hru_strmseg_down_id[ii] == 0:
#     if hru_strmseg_down_id[ii] == 0 and hru_pct_up[ii] > 0.0:
        dag_uphru.add_node(hru_down_id[ii])
        dag_uphru.add_edge(vv, hru_down_id[ii])
#     if hru_strmseg_down_id[ii] > 0:
    else:
        # HRU drains to a stream segment
        strm_node = 'S_{}'.format(hru_strmseg_down_id[ii])
        dag_uphru.add_node(strm_node, fontcolor='white', style='filled', fillcolor='blue')
        dag_uphru.add_edge(vv, strm_node)
        
        
# Output any cycles/loops
# Also update attributes for nodes which are part of a cycle
print('Cycles and loops')
for xx in nx.simple_cycles(dag_uphru):
    for yy in xx:
        dag_uphru.nodes[yy]['style'] = 'filled'
        dag_uphru.nodes[yy]['fillcolor'] = 'orange'
    print(xx)

# %%
print("Number of subset nodes: {}".format(dag_uphru.number_of_nodes()))
print("Number of subset edges: {}".format(dag_uphru.number_of_edges()))
print("Is a DAG: {}".format(nx.is_directed_acyclic_graph(dag_uphru)))

# %%
# Create plot of the graph using pydot
out_plot = nx.nx_pydot.to_pydot(dag_uphru)
out_plot.write_png(plot_filename)

# %%
aa = {'hru_up_id': hru_up_id, 
      'hru_down_id': hru_down_id,
      'hru_pct_up': hru_pct_up,
      'hru_strmseg_down_id': hru_strmseg_down_id}

# %%
new_df = pd.DataFrame.from_dict(aa)

# %%
new_df

# %%
# Where ever hru_strmseg_down_id is connected to a stream segment (e.g. > 0) 
# make sure the hru_down_id is set to zero. This is just to avoid confusion.
new_df['hru_down_id'] = np.where(new_df['hru_strmseg_down_id'].gt(0), 0, new_df['hru_down_id'])

# %%
new_df[new_df['hru_up_id'] == 2]

# %%
# For each hru_up_id find the row with the maximum hru_pct_up value.
try2 = new_df.loc[new_df.reset_index().groupby(['hru_up_id'])['hru_pct_up'].idxmax()]

# %%

# %%
# Get cascade-related parameters
hru_up_id = try2['hru_up_id'].tolist()
hru_down_id = try2['hru_down_id'].tolist()
hru_pct_up = try2['hru_pct_up'].tolist()
hru_strmseg_down_id = try2['hru_strmseg_down_id'].tolist()
tosegment = pfile.parameters['tosegment'].data

# Build the network
dag_uphru = nx.DiGraph()

# Generate the stream network
# for ii, vv in enumerate(tosegment):
#     if vv == 0:
#         dag_uphru.add_node(f'S_{ii+1}', style='filled', fontcolor='white', fillcolor='blue')
#         dag_uphru.add_node(f'S_{vv}', style='filled', fontcolor='white', fillcolor='grey')
#         dag_uphru.add_edge(f'S_{ii+1}', f'S_{vv}')
#     else:
#         dag_uphru.add_edge(f'S_{ii+1}', f'S_{vv}')
        
# Add the hru and stream network connections      
for ii, vv in enumerate(hru_up_id):
    # ii => cascade number
    # vv => hru id
    dag_uphru.add_node(vv)
    
#     if hru_pct_up[ii] > 0.0:
    if hru_strmseg_down_id[ii] == 0:
#     if hru_strmseg_down_id[ii] == 0 and hru_pct_up[ii] > 0.0:
        dag_uphru.add_node(hru_down_id[ii])
        dag_uphru.add_edge(vv, hru_down_id[ii])
#     if hru_strmseg_down_id[ii] > 0:
    else:
        # HRU drains to a stream segment
        strm_node = f'S_{hru_strmseg_down_id[ii]}'
        dag_uphru.add_node(strm_node, fontcolor='white', style='filled', fillcolor='blue')
        dag_uphru.add_edge(vv, strm_node)
        
        
# Output any cycles/loops
# Also update attributes for nodes which are part of a cycle
print('Cycles and loops')
for xx in nx.simple_cycles(dag_uphru):
    for yy in xx:
        dag_uphru.nodes[yy]['style'] = 'filled'
        dag_uphru.nodes[yy]['fillcolor'] = 'orange'
    print(xx)

# %%
print("Number of subset nodes: {}".format(dag_uphru.number_of_nodes()))
print("Number of subset edges: {}".format(dag_uphru.number_of_edges()))
print("Is a DAG: {}".format(nx.is_directed_acyclic_graph(dag_uphru)))

# %%
out_plot = nx.nx_pydot.to_pydot(dag_uphru)
out_plot.write_png(f'{workdir}/sagehen_prms_orig_grid_FIX.png')

# %%
ohdl = open('/Users/pnorton/tmp/hru_up_id.txt', 'w')

for xx in hru_up_id:
    ohdl.write(f'{xx}\n')
    
ohdl.close()

# %%
# Create array of hru_segment
dag_us = dag_uphru.reverse()
arr = np.zeros((128), dtype=np.int)

for xx in range(15):
    try:
        uniq_seg_us = set()
        pred = nx.dfs_predecessors(dag_us, f'S_{xx+1}')
        uniq_seg_us = uniq_seg_us.union(set(pred.keys()).union(set(pred.values())))

        dd = list(uniq_seg_us)
        dd.remove(f'S_{xx+1}')

        for yy in dd:
            arr[yy-1] = xx+1
    except KeyError:
        print(f'Segment {xx+1} has no HRUs connected to')

# %%
ohdl = open('/Users/pnorton/tmp/hru_segment.param', 'w')
# ohdl.write('nhm_id,hru_segment\n')

for nn,xx in enumerate(arr):
#     ohdl.write(f'{nn+1},{xx}\n')
    ohdl.write(f'{xx}\n')
    
ohdl.close()

# %%

# %%

# %%
arr

# %%

# %%
