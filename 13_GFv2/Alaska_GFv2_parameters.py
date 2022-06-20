# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %%
import glob
import os
import pandas as pd
import pydot
import networkx as nx
import numpy as np
import warnings

from collections import Counter, OrderedDict
from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterSet import ParameterSet
from pyPRMS.constants import DATATYPE_TO_DTYPE

warnings.filterwarnings('ignore')

# %%
base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/FY2022_projects/modeling_fabric'
params_dir = f'{base_dir}/gfv2.0_ak'
pdb_dir = f'{base_dir}/paramdb_v20_daymet_alaska'
plot_dir = f'{base_dir}/20220616_bad_parameters_plots'

# HRU polygons
hru_geodatabase = f'{base_dir}/NHM_19.gdb'
hru_layer_name = 'nhru'
hru_shape_key = 'hru_id'

# Segment lines
seg_geodatabase = f'{base_dir}/NHM_19.gdb'
seg_layer_name = 'nsegment'
seg_shape_key = 'segment_id'


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
    
    
def _data_it(filename):
    """Get iterator to a parameter db file.

    :returns: iterator
    """

    # Read the data
    fhdl = open(filename)
    rawdata = fhdl.read().splitlines()
    fhdl.close()
    return iter(rawdata)


def read_data(filename):
    it = _data_it(filename)
    next(it)

    data = []

    # Read the parameter values
    for rec in it:
        idx, val = rec.split(',')
        data.append(val)

    return data


# %%
# Get the list of parameter files to read
gf_files = []

file_it = glob.glob(f'{params_dir}/*.csv')
for kk in file_it:
    gf_files.append(kk)
    
gf_files.sort()

# %% [markdown]
# ## Read the geospatial fabric parameter files and create a new parameterSet

# %% tags=[] jupyter={"outputs_hidden": true}
new_ps = ParameterSet(verbose=True)

# Some dimensions are derived from particular parameters
derived_dimensions = {'nhru': 'nhm_id',
                      'nsegment': 'nhm_seg',
                      'npoigages': 'poi_gage_id',
                      'ndeplval': 'snarea_curve'}

for dd, pp in derived_dimensions.items():
    tmp_data = read_data(f'{params_dir}/{pp}.csv')
    new_ps.dimensions.add(dd, size=len(tmp_data))

# Add the constant-size dimensions
new_ps.dimensions.add('ndays', size=366)
new_ps.dimensions.add('nmonths', size=12)
new_ps.dimensions.add('one', size=1)

global_dims = new_ps.dimensions
print('='*10, 'Dimensions', '='*10)
print(global_dims)

# =======================================================
# Load the parameters
print('='*10, 'Parameters', '='*10)
for ff in gf_files:
    cname = os.path.basename(os.path.splitext(ff)[0])
    print(cname)

    try:
        if new_ps.master_parameters is not None:
            new_ps.parameters.add(cname, info=new_ps.master_parameters[cname])
        else:
            new_ps.parameters.add(cname)

        # Add the dimension sizes
        master_info = new_ps.master_parameters[cname]
        for dd in master_info.dimensions.values():
            new_ps.parameters[cname].dimensions.add(dd.name, size=global_dims[dd.name].size)

        tmp_data = pd.read_csv(ff, skiprows=0, usecols=[1],
                               dtype={1: DATATYPE_TO_DTYPE[new_ps.parameters.get(cname).datatype]}).squeeze('columns')                  

        new_ps.parameters[cname].data = tmp_data
    except ValueError:
        print(f'{cname} is not a valid PRMS parameter; skipping.')

# %% [markdown]
# ## Run basic consistency checks for all parameters

# %% tags=[]
new_ps.parameters.check()

# %%

# %%
# %%time
# Load the shapefiles for the model HRUs and segments
new_ps.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)
new_ps.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %%
pdict = dict(carea_max={'fillValue': -9999},
             covden_sum={'fillValue': -999.99},
             covden_win={'fillValue': -999.99},
             dprst_frac={'fillValue': -9999.0, 'limits': (0.0, 0.2)},
             hru_area={'fillValue': -999.99, 'limits': (0.0, 340000.0)},
             hru_deplcrv={'fillValue': 0, 'cmap': 'Pastel1'},
             hru_elev={'fillValue': -99999.0, 'limits': (0.0, 4000.0)},
             hru_percent_imperv={'fillValue': -9999.0, 'limits': (0.0, 0.3)},
             hru_slope={'fillValue': -999.99, 'limits': (0.0, 0.005)},
             rad_trncf={'fillValue': -99999.0},
             seg_slope={'fillValue': -999.99, 'limits': (0.0, 0.37), 'cmap': 'YlGnBu'},
             smidx_coef={'fillValue': -9999.0},
             soil_moist_max={'fillValue': -99999.0, 'limits': (0.0, 910.0)},
             soil_type={'fillValue': -99999, 'cmap': 'Pastel1'},
             sro_to_dprst_imperv={'fillValue': -9999.0, 'limits': (0.0, 0.4)},
             sro_to_dprst_perv={'fillValue': -9999.0})

# Some parameters have values outside the valid range; convert those values to fillValue for plotting
pdict_nofill = {'covden_win', 'hru_area', 'hru_slope', 'seg_slope'}

for pp in pdict_nofill:
    aa = new_ps.parameters[pp]
    aa.data = np.where(aa.data < aa.minimum, -999.99, aa.data)

# hru_deplcrv has fillValue and zero-values; valid range is 1 to 9
# When hru_deplcrv values are bounded we use the default value as the minimum for comparison
aa = new_ps.parameters['hru_deplcrv']
aa.data = np.where(aa.data < aa.default, 0, aa.data)

# %%
# Plot multiple parameters to files
for cparam, args in pdict.items():
    print(new_ps.parameters[cparam].stats())

    # noinspection PyTypeChecker
    cmap = args.setdefault('cmap', 'terrain')
    # noinspection PyTypeChecker
    limits = args.setdefault('limits', 'valid')
    
    try:
        new_ps.parameters[cparam].default = args['fillValue']
        new_ps.parameters.plot(cparam, limits=limits, cmap=cmap, mask_defaults='deeppink', output_dir=plot_dir)
    except KeyError:
        new_ps.parameters.plot(cparam, limits=limits, cmap=cmap, output_dir=plot_dir)

# %%
# Single plot to the screen
cparam = 'seg_slope'
cmap = pdict[cparam].setdefault('cmap', 'terrain')

print(new_ps.parameters[cparam].stats())

try:
    new_ps.parameters[cparam].default = pdict[cparam]['fillValue']
    new_ps.parameters.plot(cparam, limits=pdict[cparam]['limits'], cmap=cmap, mask_defaults='deeppink')
except KeyError:
    new_ps.parameters.plot(cparam, limits=pdict[cparam]['limits'], cmap=cmap)

# %%
# TODO: Need to check for missing parameters needed by selected PRMS modules

# %%

# %%

# %% [markdown]
# ## Create a parameter database


# %% jupyter={"outputs_hidden": true} tags=[]
# %%time
new_ps.write_paramdb(f'{pdb_dir}')

# %% [markdown]
# ## Load the parameter database and check the stream network

# %%
pdb = ParamDb(paramdb_dir=pdb_dir, verbose=True, verify=True)

# %%
nhm_id = pdb.parameters['nhm_id'].tolist()
nhm_seg = pdb.parameters['nhm_seg'].tolist()

tosegment = pdb.parameters['tosegment'].tolist()
hru_segment = pdb.parameters['hru_segment'].tolist()
# hru_segment = pdb.parameters['hru_segment_nhm'].tolist()
# tosegment = pdb.parameters['tosegment_nhm'].tolist()

poi = {}
if pdb.parameters.exists('poi_gage_segment'):
    poi_gage_id = pdb.parameters['poi_gage_id'].tolist()
    poi_gage_segment = pdb.parameters['poi_gage_segment'].tolist()
    
    for gg, ss in zip(poi_gage_id, poi_gage_segment):
        poi[ss] = gg

# %%
# %%time
# Build the stream network
num_outlets = 0
include_hrus = False

dag_streamnet = nx.DiGraph()

for ii, vv in enumerate(tosegment):
    if vv == 0 or vv not in nhm_seg:
#     if vv == 0 or tosegment[ii] == 0:
        # Color the oulet segment nodes
        dag_streamnet.add_node(nhm_seg[ii], style='filled', fontcolor='white', fillcolor='blue')
    
        # The next two lines would add a node and edge for the '0' segment
        # dag_streamnet.add_node(vv, style='filled', fontcolor='white', fillcolor='grey')
        # dag_streamnet.add_edge(nhm_seg[ii], vv)

        num_outlets += 1
    else:
        dag_streamnet.add_edge(nhm_seg[ii], tosegment[ii])
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

    for ii, vv in enumerate(hru_segment):
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
    
print(f'Number of subset nodes: {dag_streamnet.number_of_nodes()}')
print(f'Number of subset edges: {dag_streamnet.number_of_edges()}')
print(f'Number of outlet nodes: {num_outlets}')
print(f'Number of isolate nodes: {len(list(nx.isolates(dag_streamnet)))}')
print(f'Is a DAG: {nx.is_directed_acyclic_graph(dag_streamnet)}')   

# print('Isolate nodes:')
# print(list(nx.isolates(dag_streamnet)))

# %% jupyter={"outputs_hidden": true} tags=[]
check_for_disconnected_graphs(dag_streamnet)

# %%
# %%time
# Create plot of the graph using pydot
pd = nx.nx_pydot.to_pydot(dag_streamnet)
pd.write_pdf(f'{base_dir}/stream_network_ak.pdf')

# %%
