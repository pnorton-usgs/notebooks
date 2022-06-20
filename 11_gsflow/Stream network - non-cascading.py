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

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import pandas as pd
import numpy as np
import pydot
import networkx as nx

from pyPRMS.ParameterFile import ParameterFile

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/sagehen/prms6'
filename = f'{workdir}/prms_grid_v2.param'
plot_filename = f'{workdir}/sagehen_prms_grid_v2_nocascades.png'

# Load parameter file
pfile = ParameterFile(filename)

# %%
# Get cascade-related parameters
tosegment = pfile.parameters['tosegment'].data
hru_segment = pfile.parameters['hru_segment'].data

# Build the network
dag_uphru = nx.DiGraph()

# Generate the stream network
for ii, vv in enumerate(tosegment):
    dag_uphru.add_node(f'S_{ii+1}', style='filled', fontcolor='white', fillcolor='blue')
    
    if vv == 0:
        dag_uphru.add_node(f'S_{vv}', style='filled', fontcolor='white', fillcolor='grey')

    dag_uphru.add_edge(f'S_{ii+1}', f'S_{vv}')


# Add the hru connections and stream network connections
for ii, vv in enumerate(hru_segment):
    # ii => hru
    # vv => hru_segment
    
    dag_uphru.add_node(ii+1)
    if vv > 0:
        dag_uphru.add_edge(ii+1, f'S_{vv}')
    else:
        dag_uphru.nodes[ii+1]['style'] = 'filled'
        dag_uphru.nodes[ii+1]['fillcolor'] = 'gray'
        
        
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

# %%
