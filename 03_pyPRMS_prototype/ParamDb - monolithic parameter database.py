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
import numpy as np
from collections import OrderedDict
from pyPRMS.ParamDb import ParamDb

# %%
# workdir = '/Users/pnorton/tmp/tmp_paramdb'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_CONUS'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/paramdb_v10_daymet_CONUS'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_lauren_daymet_byHRU_mizu'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/streamtemp_paramdb/paramdb'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/20190927_transboundary/paramdb_merged'

workdir = '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHWobs/20211129_gm_byHWobs/paramdb_headwaters_merge'

# %%
# %%time
pdb = ParamDb(paramdb_dir=workdir, verbose=True, verify=True)

# %%
# List the dimensions
for kk in pdb.dimensions.values():
    print(f'{kk.name} = {kk.size}')

# %%
pdb.parameters.check()

# %%
pdb.parameters['snarea_thresh'].stats()

# %%
pdb.parameters['adjmix_rain'].data.ndim

# %%
print(pdb.parameters['dday_slope'].data.max())
print(pdb.parameters['dday_slope'].data.min())

# %%
pdb.parameters['soil_type'].maximum

# %%
print(pdb.parameters['albset_rna'])

# %%
param1 = pdb.parameters['gwflow_coef']
param1.stats()

# %%
pdb.parameters.get_dataframe('K_coef').head()

# %%
pdb.parameters['cecn_coef'].dimensions.keys()

# %%
pdb.parameters['cecn_coef'].dimensions.ndims

# %%
remove_list = [0, 1, 2, 4]
pdb.parameters.remove_hrus_by_index(remove_list)

# Update global dimensions
pdb.dimensions['nhru'].size -= len(remove_list)
pdb.dimensions['nssr'].size -= len(remove_list)
pdb.dimensions['ngw'].size -= len(remove_list)

# %%
print(pdb.parameters['tmax_cbh_adj'].dimensions)

# %%
nhm_id = pdb.parameters['nhm_id'].tolist()
hru_segment_nhm = pdb.parameters['hru_segment_nhm'].tolist()

# %%
# Given a subset of nhm_seg ids
new_nhm_seg = [30113,30114,30115,30116,30117,30118,30119]

# Create an ordered dictionary to which contains the reduced mapping of nhm_id to hru_segment_nhm
cc = OrderedDict()
for kk, vv in zip(nhm_id, hru_segment_nhm):
    if vv in new_nhm_seg:
        cc[kk] = vv
print(cc)  

# %%
# Use the dictionary from above to create subsets of nhm_id, hru_segment_nhm, and hru_segment

new_nhm_id = [xx for xx in cc.keys()]
print('nhm_id: {}'.format(','.join(map(str, new_nhm_id))))

new_hru_segment_nhm = [kk if kk in new_nhm_seg else 0 if kk == 0 else -1 for kk in cc.values()]
print('hru_segment_nhm: {}'.format(','.join(map(str, new_hru_segment_nhm))))

new_hru_segment = [new_nhm_seg.index(kk)+1 if kk in new_nhm_seg else 0 if kk == 0 else -1 for kk in cc.values()]
print('hru_segment: {}'.format(','.join(map(str, new_hru_segment))))

# %%
  

# %%
print(pdb.parameters['nhm_id'].data)

# %%
aa = np.delete(pdb.parameters['tmax_cbh_adj'].data, 0, axis=0)

# %%
aa[0]

# %%
pdb.parameters['tmax_cbh_adj'].dimensions.get_position('nmonth')

# %%
aa = pdb.parameters['nhm_id'].data
bb = pdb.parameters['hru_segment_nhm'].data

cc = np.column_stack((aa, bb))

# %%
cc

# %%
nhm_id_to_hru_segment_nhm = OrderedDict((nhm, hseg) for nhm, hseg in cc)

# %%
nhm_id_to_hru_segment_nhm

# %%
type(aa)

# %%
aa.ndim

# %%
type(nhm_id_to_hru_segment_nhm)

# %%
remove_list = [1, 2, 4, 10, 109951]
# %time pdb.remove_by_global_id(hrus=remove_list)

# %%
pdb.parameters['hru_deplcrv'].data.shape

# %%
pdb.parameters['snarea_curve'].data.shape

# %%
pdb.parameters['snarea_curve'].data.reshape(-1, 11)

# %%
pdb.parameters['snarea_curve'].data.reshape(-1, 11)[4, :]

# %%
with np.nditer(data_copy, op_flags=['readwrite']) as it:
    for xx in it:
#         print(xx)
        xx[...] = uniq_dict[int(xx)]

# %%
data_copy

# %%
pdb.parameters['hru_deplcrv'].data

# %%
print(pdb.dimensions)

# %%
print(pdb.parameters.get('cecn_coef').dimensions)

# %%
pdb.dimensions['ndeplval'].size / 11

# %%
possible_idx = [xx for xx in range(1, int(pdb.dimensions['ndeplval'].size / 11)+1)]

# %%
len(possible_idx)

# %%
curr_idx = pdb.parameters['hru_deplcrv'].data.tolist()
len(curr_idx)

# %%
set(possible_idx).difference(set(curr_idx))

# %%
pdb.parameters['snarea_curve'].data.shape

# %%
pdb.parameters['nhm_seg'].index_map

# %%
nhm_id = [57874, 57875, 57878, 57881, 57868, 57873, 57879, 57880, 57882, 57883, 57869, 57870, 57864, 57865]

aa = pdb.parameters.get_subset('hru_segment_nhm', nhm_id)
print(aa)
dd = pdb.parameters.get_subset('hru_segment', nhm_id)
print(dd)

# %%
nhm_seg = [30113, 30114, 30115, 30116, 30117, 30118, 30119]

bb = pdb.parameters.get_subset('tosegment_nhm', nhm_seg)
print(bb)
cc = pdb.parameters.get_subset('tosegment', nhm_seg)
print(cc)

# %%
# An empty list of IDs should return an empty array
cc = pdb.parameters.get_subset('tosegment_nhm', [])
print(cc)

# %%
aa = iter(pdb.parameters['tmax_cbh_adj'].dimensions)

# %%
for xx in aa:
    print(xx)

# %%
for aa in pdb.parameters['tmax_cbh_adj'].dimensions.names:
    print(aa)

# %%
bb = iter(pdb.parameters['tmax_cbh_adj'].dimensions.keys())

# %%
for cc in bb:
    print(cc)

# %%
pdb.parameters['sat_threshold'].check_values()

# %%
pdb.parameters['soil_type'].data.min()

# %%
bb[0][1].size

# %%
pdb.master_parameters.keys()

# %%
aa = np.array(0.8, ndmin=1)

print(aa)
print(type(aa))
print(aa.size)
print(aa.shape)
print(aa.ndim)

# %%
sorted(list(pdb.parameters.keys()))

# %%
# I[1] a = np.ones(100000)
# I[2] timeit (a == a[0]).all()
# 1000 loops, best of 3: 203 us per loop
# I[3] timeit a.min() == a.max()
# 10000 loops, best of 3: 106 us per loop
# I[4] timeit np.ptp(a)
# 10000 loops, best of 3: 106 us per loop

td = pdb.parameters['width_alpha'].data

# %timeit (td == td[0]).all()



# %%
# %timeit td.min() == td.max()

# %%
# %timeit np.ptp(td)

# %%
hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
hru_layer_name = 'nhruv11_sim30'
hru_shape_key='nhru_v11'


# %%
# %%time
pdb.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)

# pdb.parameters.shapefile_hrus('/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp', 
#                               layer_name=None, shape_key='hru_id_nat')

# %%
# %%time
pdb.parameters.plot('carea_max', output_dir=None, cmap='tab20')

# %%
pdb.parameters['dday_intcp'].dimensions.get_dimsize_by_index(0)

# %%
import pandas as pd

param_stats = []

for pp in pdb.parameters.values():
    param_stats.append(pp.stats())
    
df = pd.DataFrame.from_records(param_stats, columns=['name', 'min', 'max', 'mean', 'median'])
df.to_csv(f'/Users/pnorton/tmp/pdb_v11_stats.csv', index=False)

# %%
pdb.parameters.get('tmax_cbh_adj').data[82940,:]

# %%
pdb.degenerate_parameters()

# %%
print(pdb.parameters.get('snarea_curve'))

# %%
pdb.expand_parameter('hru_deplcrv')

# %%
# Parameter, type, units, description, minimum, maximum, default, dimensions, source
# out_str = ''
out_list = []
for pp in pdb.parameters.values():
    
#     out_str += f'{pp.name},{pp.datatype},{pp.units},{pp.description},{pp.minimum},{pp.maximum},{pp.default}'
#     out_str += ',['
    
    dims = []
    for dd in pp.dimensions.values():
        dims.append(dd.name)
#         out_str += f'{dd.name},'
#         print(f'\t{dd.name}')
    out_list.append([pp.name, pp.datatype, pp.units, pp.description, pp.minimum, pp.maximum, pp.default, dims])
#     dim_str = ','.join(dims)
#     out_str += f',[{dim_str}]\n'
#     print(pp.name, pp.datatype, pp.units, pp.description, pp.minimum, pp.maximum, pp.default)
# print(out_str)

# %%
print(pdb.parameters['K_coef'].dimensions)

# %%
import pandas as pd
col_names = ['parameter', 'datatype', 'units', 'description', 'minimum', 'maximum',
             'default', 'dimensions']
df = pd.DataFrame.from_records(out_list, columns=col_names)
df.to_csv(f'/Users/pnorton/tmp/parameters_v11.csv', sep='\t', index=False)

# %%
df

# %%
out_list

# %%
aa = pdb.parameters['fastcoef_lin']

# %%
print(aa)

# %%
aa.data[99]

# %%
idx0 = pdb.parameters['nhm_id']._value_index(100)

# %%
idx0[0]

# %%
aa.update_element(idx0[0], 0.5)
aa.data[idx0[0]]

# %%
aa.data[98:102]

# %%

# %%
