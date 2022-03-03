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
#     display_name: Python [conda env:gis]
#     language: python
#     name: conda-env-gis-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
# # %load_ext cython

# import csv
import geopandas as gpd
# import json
import numpy as np
import pandas as pd
import xarray as xr

# import matplotlib as mpl
# import matplotlib.pyplot as plt

from numpy.ma import masked
# from shapely.geometry import Polygon

# import multiprocessing as mp

# %%
def read_paramdb_file(filename):
    fhdl = open(filename)
    rawdata = fhdl.read().splitlines()
    fhdl.close()
    it = iter(rawdata)

    data = []

    # Skip header
    next(it)

    # Read the parameter values
    for rec in it:
        try:
            data.append(int(rec.split(',')[1]))
        except ValueError:
            data.append(int(float(rec.split(',')[1])))
    return data


# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/ERA5'
geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_wgs84/nhruNationalIdentifier.shp'
paramdb_dir = '/Users/pnorton/tmp/tmp_paramdb'


src_file = f'{workdir}/tmax_all_ymon_2007-2011.nc'
weights_file = f'{workdir}/era5_weights.csv'

# %%
nhm_ids = read_paramdb_file(f'{paramdb_dir}/nhm_id.csv')
# print(nhm_id)

# %%
# Read the geodatabase
# gdf = gpd.read_file(geodatabase)
# gdf.head()

# %%
# prev = 0

# for index, row in gdf.iterrows():
#     curr_id = row[wgt_key]
    
#     if curr_id - prev != 1:
#         print(f'ID skip at: {curr_id}')
#     prev = curr_id

# %% [markdown]
# ### Create weighted averages by HRU

# %%
def np_get_wval(ndata, wghts, **kwargs):
    nan_default = kwargs.get('nan', np.nan)
    cgrid_ids = wghts['grid_ids'].tolist()
    cwgts = wghts['w'].tolist()
    mdata = np.ma.masked_array(ndata[cgrid_ids], np.isnan(ndata[cgrid_ids]))

    tmp = np.ma.average(mdata, weights=cwgts)
    if tmp is masked:
        return nan_default

    return tmp

# %%
# # %%cython -a
# import numpy as np
# from numpy.ma import masked


# def np_get_wval2(ndata, weights):
# #     cdef double _tmp = 0.0
# #     cdef long[:] cgrid_ids
# #     cdef double[:] cwgts
# #     cdef double[:] mdata
    
#     nan_default = 273.15
# #     nan_default = kwargs.get('nan', 0)
#     cgrid_ids = weights['grid_ids'].to_numpy()
#     cwgts = weights['w'].to_numpy()
#     mdata = np.ma.masked_array(ndata[cgrid_ids], np.isnan(ndata[cgrid_ids]))

#     _tmp = np.ma.average(mdata, weights=weights)
#     if _tmp is masked:
#         return nan_default

#     return _tmp

# %%
# Read the weights files
wgt_df = pd.read_csv(weights_file)
wgt_key = wgt_df.columns[1]
hru_ids = wgt_df.groupby(wgt_key)
nhru = len(hru_ids)

# Read the source data file
src_data = xr.open_dataset(src_file, mask_and_scale=True)

# wgt_df.head()
# print(wgt_key)
# print(hru_ids.head())
# src_data.head()

# %%
# # %%time
# allsnow = np.zeros((nhru, 12))
# temp_default = 273.15  # Default should match units of src_data

# for month in range(12):
#     print(f'Month: {month+1}')
#     aa = src_data['t2m_snow'][month].values.flatten(order="K")
    
#     for index, row in gdf.iterrows():
#         print(f'{index}, {row[wgt_key]}')
#         break
#         # curr_id is the current nhm_id in the weights file; wgt_key is name of the column
#         c_nhm_id = row[wgt_key]
        
#         # get the rows associated with the current nhm_id
#         curr_wgts = hru_ids.get_group(c_nhm_id)
        
# #         # get array of grid cells to use from src_data
# #         cgrid_ids = curr_wgts['grid_ids'].to_numpy()
        
# #         # get array of weights associated with the grid cells
# #         cwgts = curr_wgts['w'].to_numpy()

# #         allsnow[curr_id-1, month] = np_get_wval2(aa, cgrid_ids, cwgts)
#         allsnow[curr_id-1, month] = np_get_wval(aa, curr_wgts, nan=temp_default)

#         if index % 10000 == 0:
#             print(index, curr_id)
            
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Write the data to a file
# # Convert Kelvin to Celsius
# allsnow -= 273.15
# allsnow = allsnow * 1.8 + 32.0

# # Write the output file
# varname = 'tmax_allsnow'

# with open(f'{workdir}/{varname}.csv', 'w') as ff:
#     outstr = f'$id,{varname}\n'

#     for ii, dd in enumerate(allsnow.ravel(order='F').tolist()):
#         tmp = '{:<20f}'.format(dd).rstrip('0 ')
#         if tmp[-1] == '.':
#             tmp += '0'
#         outstr += '{},{}\n'.format(ii+1, tmp)

#     ff.write(outstr)            

# %%
# %%time
allsnow = np.zeros((nhru, 12))
temp_default = 273.15  # Default should match units of src_data
# nhm_ids = gdf[wgt_key].to_list()


for month in range(12):
    print(f'Month: {month+1}')
    aa = src_data['t2m_snow'][month].values.flatten(order="K")
    
    for c_nhm_id in nhm_ids:
        # get the rows associated with the current nhm_id
        curr_wgts = hru_ids.get_group(c_nhm_id)

        allsnow[c_nhm_id-1, month] = np_get_wval(aa, curr_wgts, nan=temp_default)

        if c_nhm_id % 10000 == 0:
            print(c_nhm_id)
            
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write the data to a file
# Convert Kelvin to Celsius
allsnow -= 273.15
allsnow = allsnow * 1.8 + 32.0

# Write the output file
varname = 'tmax_allsnow'

with open(f'{workdir}/{varname}.csv', 'w') as ff:
    outstr = f'$id,{varname}\n'

    for ii, dd in enumerate(allsnow.ravel(order='F').tolist()):
        tmp = '{:<20f}'.format(dd).rstrip('0 ')
        if tmp[-1] == '.':
            tmp += '0'
        outstr += '{},{}\n'.format(ii+1, tmp)

    ff.write(outstr) 

# %%
# %%time
allrain = np.zeros((nhru, 12))
temp_default = 274.65  # Default should match units of src_data

for month in range(12):
    print(f'Month: {month+1}')
    aa = src_data['t2m_rain'][month].values.flatten(order="K")
    
    for c_nhm_id in nhm_ids:
        # get the rows associated with the current nhm_id
        curr_wgts = hru_ids.get_group(c_nhm_id)

        allrain[c_nhm_id-1, month] = np_get_wval(aa, curr_wgts, nan=temp_default)

        if c_nhm_id % 10000 == 0:
            print(c_nhm_id)

# Replace any nans with the default minimum allrain temperature
# np.nan_to_num(allrain, nan=274.65, copy=False)

# Convert Kelvin to Fahrenheit 
allrain -= 273.15
allrain = allrain * 1.8 + 32.0

# %%
# Create the allrain_offset variable
allrain_offset = allrain - allsnow

# %%
# np.nan_to_num(allrain_offset, nan=1.0, copy=False)
# allrain_offset

# %%
# Write the output file
varname = 'tmax_allrain_offset'

with open(f'{workdir}/{varname}.csv', 'w') as ff:
    outstr = f'$id,{varname}\n'

    for ii, dd in enumerate(allrain_offset.ravel(order='F').tolist()):
        tmp = '{:<20f}'.format(dd).rstrip('0 ')
        if tmp[-1] == '.':
            tmp += '0'
        outstr += '{},{}\n'.format(ii+1, tmp)

    ff.write(outstr) 

# %% [markdown]
# ### Parallel work

# %%
# import SharedArray

# def multiprocess_data(fake_data):
#     data = SharedArray.create('data', (x_dim, y_dim))
    
#     def calc_row(i):
#         for j, data_point in enumerate(fake_data[i]):
#             if data_point == None:
#                 data[i, j] = 0
#             data[i, j] = math.exp(math.sqrt(data_point))
    
#     processes = []
#     for i in range(len(fake_data)):
#         process = multiprocessing.Process(target=calc_row, args=(i,))
#         processes.append(process)
#         process.start()
    
#     for process in processes:
#         process.join()
    
#     return data

# %%
# Parallel
# import SharedArray

# def multiprocess_data(nhru):
#     # Read the weights files
#     wgt_df = pd.read_csv(weights_file)
#     wgt_key = wgt_df.columns[1]
#     hru_ids = wgt_df.groupby(wgt_key)
#     nhru = len(hru_ids)
    
#     # Read the source data file
#     src_data = xr.open_dataset(src_file, mask_and_scale=True)

#     temp_default = 273.15  # Default should match units of src_data
    
#     data = SharedArray.create('data', (nhru, 12))
    
#     def np_get_wval(ndata, wghts, curr_id, month, **kwargs):
#         nan_default = kwargs.get('nan', 0)
#         cgrid_ids = wghts['grid_ids'].tolist()
#         cwgts = wghts['w'].tolist()
#         mdata = np.ma.masked_array(ndata[cgrid_ids], np.isnan(ndata[cgrid_ids]))

#         data[curr_id-1, month] = np.ma.average(mdata, weights=cwgts)
# #         if tmp is masked:
# #             return nan_default
# #         return tmp

#     cpus = 6
#     shp_chunk = np.array_split(gdf, cpus)
#     pool = mp.Pool(processes=cpus)
#     chunk_processes = [pool.apply_async(outer_loop, args=(chunk, ncfcells)) for chunk in shp_chunk]
    
#     for month in range(12):
#         print(f'Month: {month+1}')
#         aa = src_data['t2m_snow'][month].values.flatten(order="K")

#         for index, row in gdf.iterrows():
#             curr_id = row[wgt_key]
#             curr_wgts = hru_ids.get_group(curr_id)

#             data[curr_id-1, month] = np_get_wval(aa, curr_wgts, curr_id, month, nan=temp_default)

#             if index % 10000 == 0:
#                 print(index, curr_id)



# %%
# Convert Kelvin to Celsius
allsnow -= 273.15
allsnow = allsnow * 1.8 + 32.0

# Write the output file
varname = 'tmax_allsnow'

with open(f'{workdir}/{varname}.csv', 'w') as ff:
    outstr = f'$id,{varname}\n'

    for ii, dd in enumerate(allsnow.ravel(order='F').tolist()):
        tmp = '{:<20f}'.format(dd).rstrip('0 ')
        if tmp[-1] == '.':
            tmp += '0'
        outstr += '{},{}\n'.format(ii+1, tmp)

    ff.write(outstr)

# %% [markdown]
# ## Scratch stuff

# %%
# This is base on https://swanlund.space/parallelizing-python

# def outer_loop(shp_chunk, src_netcdf):
#     results = []
    
#     # Have to set spatial_index here instead of passing it in as an argument,
#     # otherwise the intersections don't work (all return empty sets)
#     # see: https://github.com/geopandas/geopandas/issues/1259
#     spatial_index = src_netcdf.sindex
    
#     for index, row in shp_chunk.iterrows():
#         count = 0

#         possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))
        
#         if len(possible_matches_index) != 0:
#             possible_matches = src_netcdf.iloc[possible_matches_index]
#             precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]

#             if len(precise_matches) != 0:
#                 res_intersection = gpd.overlay(shp_chunk.loc[[index]], precise_matches, how='intersection')

#                 for nindex, row in res_intersection.iterrows():
#                     tmpfloat = np.float(res_intersection.area.iloc[nindex] / shp_chunk.loc[[index], 'geometry'].area)
#                     results.append([np.int(precise_matches.index[count]), np.int(row['hru_id_nat']), tmpfloat])
#                     count += 1
#     return results

# def parallelize():
#     cpus = 6
#     # cpus = mp.cpu_count()

#     isect_dist = []
    
#     spatial_index = ncfcells.sindex
    
#     shp_chunk = np.array_split(gdf, cpus)
    
#     pool = mp.Pool(processes=cpus)
    
#     chunk_processes = [pool.apply_async(outer_loop, args=(chunk, ncfcells)) for chunk in shp_chunk]

#     intersect_results = [chunk.get() for chunk in chunk_processes]
#     isect_dist.append(intersect_results)
    
#     return isect_dist

# %%
from datetime import datetime
aa = datetime(1998, 1, 1)
bb = datetime(2016, 1, 31)
aa + (bb - aa) / 2.

# %%
