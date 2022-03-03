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
#     display_name: Python [conda env:bandit]
#     language: python
#     name: conda-env-bandit-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import csv
import geopandas as gpd
import json
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib as mpl
import matplotlib.pyplot as plt

from numpy.ma import masked
from shapely.geometry import Polygon

import multiprocessing as mp
# import dask
# from dask.distributed import Client

# %%
np.show_config()

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/ERA5'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/snodas/unmasked'

geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
# geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_wgs84/nhruNationalIdentifier.shp'
# geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp'

# %%
ds = xr.open_dataset(f'{workdir}/SWE_median_of_max_yearly_inch.nc', mask_and_scale=False)
# ds = xr.open_dataset(f'{workdir}/land_sea_mask_upk.nc', mask_and_scale=False)
print(ds)

# %%
# Read the geodatabase
gdf = gpd.read_file(geodatabase)
# gdf.reset_index(drop=True).head()
gdf.head()

# %%

# %%
# Create some variables from the netcdf file for use later
# print('The meta data is: \n', json.dumps(ds.attrs, indent=4))

# lathandle = ds['latitude']
# lonhandle = ds['longitude']
lathandle = ds['lat']
lonhandle = ds['lon']
datahandle = ds['SWE']
# timehandle = ds['time']
# datahandle = ds['lsm']

# %%
print('Data attributes, sizes, and coords \n')
print('Data sizes are: \n', datahandle.sizes)
print('\nData coords are: \n', datahandle.coords)

# %%
df = pd.DataFrame({'grid': datahandle.isel(time=0).values.flatten()})
df.info()

# %%
df.head()

# %%
datahandle.isel(time=0).plot()

# %%
lon, lat = np.meshgrid(lonhandle, lathandle)
df = pd.DataFrame({'grid': datahandle.isel(time=0).values.flatten()})

# res is the 'resolution' (e.g. gridspacing in degrees)
# gridspacing = 0.25
gridspacing = 0.04166666
res = gridspacing / 2.0
numcells = lon.size
# numcells = np.shape(lat)[0] * np.shape(lon)[1]
poly = []
# index = np.zeros(numcells)
index = []
count = 0

print(f'res = {res}')
print(f'numcells = {numcells}')
print(f'index = {index}')
print(f'count = {count}')

for i in range(np.shape(lat)[0]):
    for j in range(np.shape(lon)[1]):
        lat_point_list = [lat[i, j] - res, lat[i, j] + res, lat[i, j] + res, lat[i, j] - res]
        lon_point_list = [lon[i, j] + res, lon[i, j] + res, lon[i, j] - res, lon[i, j] - res]
        poly.append(Polygon(zip(lon_point_list, lat_point_list)))
        index.append(count)
        count += 1
        
# ncfcells = gpd.GeoDataFrame(df, index=index, crs=ds['crs'], geometry=poly)
ncfcells = gpd.GeoDataFrame(df, index=index, geometry=poly)
ncfcells.head()

# %%

# %%
print(lat_point_list)
print(lon_point_list)
print(lonhandle)


# %%
# spatial_index is built from the source netcdf file
spatial_index = ncfcells.sindex
tcount = 0

with open(f'{workdir}/tmp_weights2.csv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for index, row in gdf.iterrows():
        # Iterate through all rows in the NHM HRU-shapefile
        count = 0

        if tcount == 0:
            writer.writerow(['grid_ids', 'hru_id_nat', 'w'])
#             writer.writerow(['grid_ids', 'nhm_id', 'hru_id_nat', 'w'])
            
        possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))
#         print(possible_matches_index)
        
        if len(possible_matches_index) != 0:
            possible_matches = ncfcells.iloc[possible_matches_index]
            precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]

            if len(precise_matches) != 0:
                res_intersection = gpd.overlay(gdf.loc[[index]], precise_matches, how='intersection')
                
                for nindex, row in res_intersection.iterrows():
                    tmpfloat = np.float(res_intersection.area.iloc[nindex] / gdf.loc[[index], 'geometry'].area)
                    writer.writerow([np.int(precise_matches.index[count]), np.int(row['hru_id_nat']), tmpfloat])
#                     writer.writerow([np.int(precise_matches.index[count]), np.int(row['nhm_id']), np.int(row['hru_id_nat']), tmpfloat])
                    count += 1
                tcount += 1
                
                if tcount%1000 == 0:
                    print(tcount, index)
        else:
            print('no intersection: ', index, np.int(row['nhm_id']))
            
        print(f'index: {index}')
        print(row)
        break

# %% [markdown]
# ## Begin parallel experiment

# %%
# possible_matches_index (uses spatial_index and row)
# ncfcells


# def compute_weights(spatial_index, index, row, src_gdf, shp_gdf):
#     possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))
#     results = None
#     count = 0
    
#     if len(possible_matches_index) != 0:
#         possible_matches = ncfcells.iloc[possible_matches_index]
#         precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]

#         if len(precise_matches) != 0:
#             res_intersection = gpd.overlay(gdf.loc[[index]], precise_matches, how='intersection')
#             results = []
            
#             for nindex, row in res_intersection.iterrows():
#                 tmpfloat = np.float(res_intersection.area.iloc[nindex] / gdf.loc[[index], 'geometry'].area)
#                 results.append([np.int(precise_matches.index[count]), np.int(row['hru_id_nat']), tmpfloat])
# #                 writer.writerow([np.int(precise_matches.index[count]), np.int(row['hru_id_nat']), tmpfloat])
#                 count += 1
#     return results


# def compute_weights2(index, row, possible_matches, shp_gdf):
#     results = None
#     count = 0

#     precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]

#     if len(precise_matches) != 0:
#         res_intersection = gpd.overlay(shp_gdf.loc[[index]], precise_matches, how='intersection')
#         results = []

#         for nindex, row in res_intersection.iterrows():
#             tmpfloat = np.float(res_intersection.area.iloc[nindex] / shp_gdf.loc[[index], 'geometry'].area)
#             results.append([np.int(precise_matches.index[count]), np.int(row['hru_id_nat']), tmpfloat])
#             count += 1
#     return results

def outer_loop(shp_chunk, src_netcdf):
    results = []
    
    # Have to set spatial_index here instead of passing it in as an argument,
    # otherwise the intersections don't work (all return empty sets)
    # see: https://github.com/geopandas/geopandas/issues/1259
    spatial_index = src_netcdf.sindex
    
    for index, row in shp_chunk.iterrows():
        count = 0

        possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))
        
        if len(possible_matches_index) != 0:
            possible_matches = src_netcdf.iloc[possible_matches_index]
            precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]

            if len(precise_matches) != 0:
                res_intersection = gpd.overlay(shp_chunk.loc[[index]], precise_matches, how='intersection')

                for nindex, row in res_intersection.iterrows():
                    tmpfloat = np.float(res_intersection.area.iloc[nindex] / shp_chunk.loc[[index], 'geometry'].area)
                    results.append([np.int(precise_matches.index[count]), np.int(row['hru_id_nat']), tmpfloat])
                    count += 1
    return results
    
def parallelize():
    cpus = 4
#     cpus = mp.cpu_count()
    isect_dist = []
    
    spatial_index = ncfcells.sindex
    
    shp_chunk = np.array_split(gdf, cpus)
    
    pool = mp.Pool(processes=cpus)
    
    chunk_processes = [pool.apply_async(outer_loop, args=(chunk, ncfcells)) for chunk in shp_chunk]

    intersect_results = [chunk.get() for chunk in chunk_processes]
    isect_dist.append(intersect_results)
    
    return isect_dist


# %%
# %%time
final_res = parallelize()

# %%

with open(f'{workdir}/tmp_weights2.csv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(['grid_ids', 'hru_id_nat', 'w'])
    
    for pp in final_res:
        for xx in pp:
            for dd in xx:
                writer.writerow([dd[0], dd[1], dd[2]])


# %%
# %%time
# 2020-02-24 PAN: Experimenting with parallel processing

# gdf = dask.delayed(gpd.read_file(geodatabase))

# ncfcells = dask.delayed(gpd.GeoDataFrame(df, index=index, geometry=poly))
spatial_index = ncfcells.sindex

final_results = []
tcount = 0

for index, row in gdf.iterrows():
    possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))
    
    if len(possible_matches_index) != 0:
        possible_matches = ncfcells.iloc[possible_matches_index]
        
        # Iterate through all rows in the NHM HRU-shapefile
        out_weights = compute_weights2(index, row, possible_matches, gdf)
        final_results.append(out_weights)
    #     print(out_weights)
        tcount += 1
        if tcount == 100:
            break

# %%

# %% [markdown]
# ## End parallel experiment

# %%
# How many groups do we have?
print("Number of groups:", len(spatial_index.leaves()), '\n')

# Print some basic info for few of them
n_iterations = 10
for i, group in enumerate(spatial_index.leaves()):
    group_idx, indices, bbox = group
    print("Group", group_idx, "contains ", len(indices), "geometries, bounding box:", bbox)
    i+=1
    if i == n_iterations:
        break

# %%
# gdf.loc[[index]]
row

# %%
possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))
print(possible_matches_index)

# %%
possible_matches = ncfcells.iloc[possible_matches_index]
precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]
print(possible_matches)
print(precise_matches)

# %%

# %%
for index, row in gdf.iterrows():
    count = 0

    possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))
    print(possible_matches_index)

    if len(possible_matches_index) != 0:
        possible_matches = ncfcells.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]

        if len(precise_matches) != 0:
            res_intersection = gpd.overlay(gdf.loc[[index]], precise_matches, how='intersection')

            for nindex, row in res_intersection.iterrows():
                tmpfloat = np.float(res_intersection.area.iloc[nindex] / gdf.loc[[index], 'geometry'].area)
#                 writer.writerow([np.int(precise_matches.index[count]), np.int(row['nhm_id']), np.int(row['hru_id_nat']), tmpfloat])
                count += 1
            tcount += 1

#             if tcount%100 == 0:
#                 print(tcount, index)
    else:
        print('no intersection: ', index, np.int(row['hru_id_nat']))

# %%
spatial_index

# %%
pd.DataFrame({'grid': datahandle.isel(time=0).values.flatten()})


# %% [markdown]
# ### Create weighted averages by HRU

# %%
def np_get_wval(ndata, wghts, **kwargs):
    """
    Returns weighted average of ndata with weights = grp
    1) mdata = the subset of values associated with the gridmet id's that are mapped to hru_id.
    2) Some of these values may have nans if the gridmet id is outside of conus so only return values
    that are inside of conus
    3) this means that hru's that are entirely outside of conus will return nans which will ultimately,
    outside of this function get assigned zero's.
    4) the value is assigned the weighted average
    :param ndata: float array of data values
    :param wghts: float array of weights
    :param hru_id hru id number
    :return: numpy weighted averaged - masked to deal with nans associated with
            ndata that is outside of the conus.
    """
    nan_default = kwargs.get('nan', 0)
    cgrid_ids = wghts['grid_ids'].tolist()
    cwgts = wghts['w'].tolist()
    mdata = np.ma.masked_array(ndata[cgrid_ids], np.isnan(ndata[cgrid_ids]))

    tmp = np.ma.average(mdata, weights=cwgts)
    if tmp is masked:
#         print('returning masked value', cgrid_ids, mdata, cwgts)
#         print(tmp, tmp.mask, tmp.fill_value, tmp.data)
        return nan_default

    return tmp


# %%
wgt_df = pd.read_csv(f'{workdir}/era5_weights.csv')
wgt_df.head()

# %%
wgt_key = wgt_df.columns[1]
print(wgt_key)

hru_ids = wgt_df.groupby(wgt_key)
print(hru_ids.head())

nhru = len(hru_ids)

# %%

# %%
src_data = xr.open_dataset(f'{workdir}/crap1_monthly.nc', mask_and_scale=True)
src_data.head()

# %%

# %%

# %%
allsnow = np.zeros((nhru, 12))
temp_default = 273.15  # Default should match units of src_data

for index, row in gdf.iterrows():
    curr_id = row[wgt_key]
    curr_wgts = hru_ids.get_group(curr_id)
    
    for month in range(12):
        aa = src_data['t2m_snow'][month].values.flatten(order="K")
        allsnow[curr_id-1, month] = np_get_wval(aa, curr_wgts, nan=temp_default)

    if index % 10000 == 0:
        print(index, curr_id)

# %%
# Try changing the loop order
allsnow = np.zeros((nhru, 12))
temp_default = 273.15  # Default should match units of src_data

for month in range(12):
    print(f'Month: {month+1}')
    aa = src_data['t2m_snow'][month].values.flatten(order="K")
    
    for index, row in gdf.iterrows():
        curr_id = row[wgt_key]
        curr_wgts = hru_ids.get_group(curr_id)

        allsnow[curr_id-1, month] = np_get_wval(aa, curr_wgts, nan=temp_default)

        if index % 10000 == 0:
            print(index, curr_id)

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
curr_id

# %%
curr_wgts

# %%
aa = src_data['t2m_snow'][0].values.flatten(order="K")   
aa

# %%
type(month)

# %%
aa.data

# %%
csel = curr_wgts['grid_ids'].tolist()

# %%
cdata = aa[csel]

# %%
mdata = np.ma.masked_array(cdata, np.isnan(cdata))

# %%
mdata

# %%
tmp = np.ma.average(mdata, weights=curr_wgts['w'])

# %%
tmp

# %%
tmp is masked

# %%
curr_wgts['w'].tolist()

# %%
curr_wgts['w']

# %%
tmp -= 273.15
tmp *= 1.8 + 32.0
print(tmp)

# %%
