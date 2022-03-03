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
from collections import OrderedDict
from pyPRMS.ParamDb import ParamDb
# from pyPRMS.ParameterFile import ParameterFile

import networkx as nx

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_CONUS'
workdir = '/Users/pnorton/tmp/check_paramdb_v111'
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

seg_to_hru = params.seg_to_hru
if 0 in seg_to_hru:
    # We don't want any non-routed HRUs
    del seg_to_hru[0]

# %%
len(seg_to_hru)

# %% [markdown]
# ## Compute the HRU area for each segment

# %%
# Sum up the HRU area for each segment
# This is used as a lookup table when computing the cumulative
# area by segment
area_by_seg = OrderedDict()

for seg, harea in zip(hru_segment_nhm, hru_area):
    if seg in area_by_seg:
        area_by_seg[seg] += harea
    else:
        area_by_seg[seg] = harea

# %% [markdown]
# ## Build the stream network

# %%
dag_streamnet = params.stream_network(tosegment='tosegment_nhm', seg_id='nhm_seg')

print(f'Number of NHM downstream nodes: {dag_streamnet.number_of_nodes()}')
print(f'Number of NHM downstream edges: {dag_streamnet.number_of_edges()}')
print(f'Is a DAG: {nx.is_directed_acyclic_graph(dag_streamnet)}')    

print("Isolate nodes:")
print(list(nx.isolates(dag_streamnet)))

# %% [markdown]
# ## Compute cumulative drainage area by stream segment

# %%
# Maximum area allowed for headwater in square kilometers
# Included headwaters will have an area less than this.
max_area = 3000

# Roland's original script used 0.00404685642 for conversion
acre_to_sq_km = 0.0040468564224

# Create upstream graph
dag_us = dag_streamnet.reverse()

dsmost_segs = nhm_seg
print(f'Number of segments: {len(dsmost_segs)}')

out_seg_cum_area = OrderedDict()
out_seg_hw = OrderedDict()

for dsmost_seg in dsmost_segs:
    uniq_seg_us = set()
    
    try:
        pred = nx.dfs_predecessors(dag_us, dsmost_seg)
        uniq_seg_us = uniq_seg_us.union(set(pred.keys()).union(set(pred.values())))
    except KeyError:
        print('\nKeyError: Segment {} does not exist in stream network'.format(dsmost_seg))

    dag_ds_subset = dag_streamnet.subgraph(uniq_seg_us).copy()

    node_outlets = [ee[0] for ee in dag_ds_subset.edges()]
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
        
    # Convert acres to square kilometers
    csum *= acre_to_sq_km
    
    if 0.0 < csum <= max_area:
        out_seg_cum_area[dsmost_seg] = csum
        out_seg_hw[dsmost_seg] = toseg_idx


# %%
# Build a list of HW segments that are ordered from largest headwater to smallest
largest_hw = 0

for kk, vv in out_seg_hw.items():
    if largest_hw < len(vv):
        largest_hw = len(vv)

print(f'Maximum number of segments in a headwater: {largest_hw}')

seg_by_size = []

while largest_hw != 0 :
    for kk, vv in out_seg_hw.items():
        if len(vv) == largest_hw:
            seg_by_size.append(kk)
    largest_hw -= 1


# %%
# Remove headwaters that are a subset of larger headwaters
# NOTE: This is an inefficient approach
# for xx in seg_by_size:
#     remove_list = []
#     for kk, vv in out_seg_hw.items():
#         if kk == xx:
#             continue
#         if xx in out_seg_hw:
#             if set(vv).issubset(set(out_seg_hw[xx])):
#                 remove_list.append(kk)

#     for rr in remove_list:
#         del out_seg_hw[rr]
        

# %%
# Remove headwaters that are a subset of larger headwaters
for xx in seg_by_size:
    remove_list = []
    
    if xx in out_seg_hw:
        for vv in out_seg_hw[xx]:
            if xx != vv and vv in out_seg_hw:
                remove_list.append(vv)

        for rr in remove_list:
            del out_seg_hw[rr]
            
print(f'Number of unique headwater areas less than or equal to {max_area} sq km: {len(out_seg_hw)}')

# %%
# Output two files:
# 1) The headwater segment file (used for headwater model extractions)
# 2) Headwater IDs by HRU (used for mapping purposes)

hru_by_hw = OrderedDict()

fhdl_hw_segs = open(f'{workdir}/hw_segs.csv', 'w')
fhdl_hw_segs.write('hw_id,start_seg,child_segs\n')

fhdl_hru_hw = open(f'{workdir}/hw_hrus.csv', 'w')
fhdl_hru_hw.write('nhm_id,hw_id\n')

cnt = 1
for kk, vv in out_seg_hw.items():
    out_seg_str = ','.join([str(xx) for xx in vv])
    fhdl_hw_segs.write(f'{cnt},{kk},{out_seg_str}\n')
    
    for ss in vv:
        if ss in seg_to_hru:
            # Some segments have no HRUs. If they have HRUs we'll process them.
            for hh in seg_to_hru[ss]:
                hru_by_hw[hh] = cnt
                fhdl_hru_hw.write(f'{hh},{cnt}\n')
    cnt += 1
    
fhdl_hw_segs.close()
fhdl_hru_hw.close()

# %%
len(out_seg_hw)

# %%
