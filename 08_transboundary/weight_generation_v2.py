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
#     display_name: Python [conda env:gis_38]
#     language: python
#     name: conda-env-gis_38-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import geopandas as gpd
import pandas as pd
import numpy as np
import math
import time

import xarray as xr
from shapely.geometry import Polygon
import multiprocessing as mp
import psutil
import os
import resource
import sys

import matplotlib

# %%
from shapely import speedups
speedups.available
speedups.enabled

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/ERA5'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/SNODAS/unmasked/test_on_subset'
# src_netcdf = f'{workdir}/land_sea_mask_upk.nc'
# src_gridspacing = 0.25  # Gridspacing of src_netcdf in degrees

geodatabase = f'{workdir}/nhruv11_sim30_WGS84_subset1.gpkg'
# geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
layer_name = 'nhruv11_sim30_WGS84_subset1'  # Layer name if geodatabase otherwise None
shape_key = 'nhru_v11'  # Name of HRU id attribute

gpkg_filename = f'/Users/pnorton/Projects/National_Hydrology_Model/datasets/SNODAS/unmasked/tmp/snodas_grid_horiz_subset1.gpkg'
gpkg_layer_name = 'snodas_grid_horiz_subset1'
# gpkg_filename = f'{workdir}/snodas_grid_subset1.gpkg'
# gpkg_layer_name = 'snodas_grid_subset1'


weight_file = f'{workdir}/snodas_weights_horiz1.csv'


# %%
def get_mem_usage():
    mem_size = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    if sys.platform == 'darwin':
        # mac OS returns ru_maxrss in bytes (Linux returns it in kilobytes)
        mem_size *= 1024
    return mem_size


def convert_size(size_kbytes):
    # from https://python-forum.io/Thread-Convert-file-sizes-will-this-produce-accurate-results
    if size_kbytes == 0:
        return "0 KB"

    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_kbytes, 1024)))
    power = math.pow(1024, i)
    size = round(size_kbytes / power, 2)
    return f'{size} {size_name[i-1]}'


def create_polygons_v2(lats, lons, grid_spacing, data):

    def get_poly(row, lats, lons):
        ii = row.name * 4
        return Polygon(zip(lons[ii:ii+4], lats[ii:ii+4]))

    t1 = time.time()
    lon, lat = np.meshgrid(lons, lats)
    print(f'meshgrid: {time.time() - t1}')

    # res is half of the 'resolution' (e.g. gridspacing in degrees)
    res = grid_spacing / 2.0

    # For weight generation the datahandle variable is only used for creating
    # the geodataframe from the source netcdf file.
    t1 = time.time()
    df = pd.DataFrame({'grid': data})
#     ddf = dd.from_pandas(df, npartitions=4)
    print(f'DataFrame: {time.time() - t1}')
    print('    mem usage:', convert_size(get_mem_usage()))

    # Create polygon features of the grid-cell bounding boxes
    t1 = time.time()
    lats1 = lat[:, :] - res
    lats2 = lat[:, :] + res

    yy = np.dstack((lats1, lats2, lats2, lats1)).flatten()

    lons1 = lon[:, :] + res
    lons2 = lon[:, :] - res

    xx = np.dstack((lons1, lons1, lons2, lons2)).flatten()
    print(f'numpy: {time.time() - t1}')
    print('    mem usage:', convert_size(get_mem_usage()))

    t1 = time.time()
    
#     df['geometry'] = get_poly(lats=yy, lons=xx)
#     ddf['geometry'] = ddf.apply(get_poly, lats=yy, lons=xx, axis=1)
#     ddf['geometry'] = ddf.map_partitions(get_poly, lats=yy, lons=xx, meta=(None, 'int64')).compute()
    df['geometry'] = df.apply(get_poly, lats=yy, lons=xx, axis=1)
    print(f'df_geometry: {time.time() - t1}')
    print('    mem usage:', convert_size(get_mem_usage()))

    t1 = time.time()
    ncfcells = gpd.GeoDataFrame(df, crs='epsg:4326')
    print(f'GeoDataFrame: {time.time() - t1}')
    print('    mem usage:', convert_size(get_mem_usage()))
    # print(ncfcells.memory_usage(deep=True))
    print('..done.')

    return ncfcells


# %%
# Open the source netcdf file
ds = xr.open_dataset(src_netcdf, mask_and_scale=False)

# Change the variable names to match the src_netcdf file
lathandle = ds['latitude']
lonhandle = ds['longitude']
datahandle = ds['lsm']

# %%
print('Creating ncfcells from dataframe')
ncfcells = create_polygons_v2(lathandle, lonhandle, src_gridspacing, datahandle.isel(time=0).values.flatten())
print('  ... done.')
print('    mem usage:', convert_size(get_mem_usage()))
print(ncfcells.info())

# %%
# 2020-04-03 PAN: This took about an hour to complete; created 8 GB file
ncfcells.to_file(f'{workdir}/GIS/era5_verify.gpkg', layer='geom', driver='GPKG')

# %%
print('Starting mem usage:', convert_size(get_mem_usage()))
print('Reading HRUs')
gdf = gpd.read_file(geodatabase, layer=layer_name)
print('  ... done.')
print('    mem usage:', convert_size(get_mem_usage()))

# # The geospatial fabric v1.1 is in Albers, we need to re-project it to WGS84.
# print(f'Re-projecting {layer_name} to WGS84 (lat/lon)')
# gdf = gdf.to_crs(epsg=4326)
# print('  ...done.')
# print('    mem usage:', convert_size(get_mem_usage()))

# %%
print('Starting mem usage:', convert_size(get_mem_usage()))
t1 = time.time()
ncfcells = gpd.read_file(gpkg_filename, layer=gpkg_layer_name)
print(f'GeoDataFrame: {time.time() - t1}')
print('    mem usage:', convert_size(get_mem_usage()))

# %%
t1 = time.time()
spatial_index = ncfcells.sindex
print(f'gpd.sindex: {time.time() - t1}')
print('    mem usage:', convert_size(get_mem_usage()))


# %%
def process_weights(shp_chunk, shape_key, src_netcdf, queue):
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

                for nindex, irow in res_intersection.iterrows():
                    tmpfloat = np.float(res_intersection.area.iloc[nindex] / shp_chunk.loc[[index], 'geometry'].area)
                    results.append([np.int(precise_matches.index[count]), np.int(irow[shape_key]), tmpfloat])
                    count += 1
    print(convert_size(get_mem_usage()))
    queue.put(results)
    return results


# %%

# %%
# results = []

outer_t1 = time.time()
t1 = time.time()

with open(weight_file, 'w') as f:
    # gdf HRU shapefile
    # ncfcells is source netcdf grid geometry

    f.write(f'grid_ids,{shape_key},w\n')
        
    for index, row in gdf.iterrows():
        count = 0

        if index % 100 == 0:
            print(f'Index: {index}     ({time.time() - t1} seconds)')
            t1 = time.time()
        
        possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))
        # print(f'possible_matches_index: {time.time() - t1}')

        if len(possible_matches_index) != 0:
            possible_matches = ncfcells.iloc[possible_matches_index]

            precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]

            if len(precise_matches) != 0:
                res_intersection = gpd.overlay(gdf.loc[[index]], precise_matches, how='intersection')

                for nindex, irow in res_intersection.iterrows():
                    tmpfloat = np.float(res_intersection.area.iloc[nindex] / gdf.loc[[index], 'geometry'].area)
                    srcidx = np.int(precise_matches.iloc[count]['grid_horiz'])
                    
                    f.write(f'{srcidx},{np.int(irow[shape_key])},{tmpfloat}\n')
#                     f.write(f'{np.int(precise_matches.index[count])},{np.int(irow[shape_key])},{tmpfloat}\n')
                    count += 1

print(f'outer_time: {time.time() - outer_t1}')

# %%
possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))
possible_matches_index

# %%
possible_matches

# %%

ax = res_intersection.plot(cmap='tab10')

# Plot the current HRU
ax = gdf.loc[[index], 'geometry'].plot(ax=ax, facecolor='none', edgecolor='k')

# Plot the grid cells that intersect the current HRU
precise_matches.plot(ax=ax, facecolor='none', edgecolor='k')

# %%
possible_matches.plot(facecolor='none', edgecolor='k')

# %%
precise_matches.plot(facecolor='none', edgecolor='k')

# %%
np.int(precise_matches.index[count-1])

# %%
np.int(precise_matches.iloc[count-1]['grid_vert'])

# %%
precise_matches

# %%
wgt_df = pd.read_csv(weight_file)
wgt_key = wgt_df.columns[1]
hru_ids = wgt_df.groupby(wgt_key)
nhru = len(hru_ids)

print(f'wgt_key: {wgt_key}')
print(f'num HRUs: {nhru}')

# %%
sorted(list(set(wgt_df['nhru_v11'].tolist())))

# %%
totals = hru_ids.sum()

# %%
pr1 = totals[totals.w < 0.9999999]

# %%
for xx in pr1.iterrows():
    print(xx)
    

# %%
pr1.head()

# %%
pr1.to_csv(f'{workdir}/under_one.csv')

# %%
pr2 = totals[totals.w > 1.0000000001]
pr2.to_csv(f'{workdir}/over_one.csv')

# %%
pr2 = totals[totals.w > 1.0000000001]

# %%
pr2.head()

# %%
ncfcells.head()

# %%
a = np.array([[1,2], [3,4]])
a.flatten(order='F')

# %%
ss = set(crap)
len(ss)

# %%
aa = [[l, crap.count(l)] for l in set(crap)]

# %%
total = 0
for xx in aa:
    if xx[1] > 1:
        total += 1
        print(f'{xx[0]} duplicated {xx[1]} times.')
print(f'Total duplicated HRUs = {total}')
print(f'Total unique HRUs = {len(ss)}')

# %%
# Trying to create a simple test grid with incrementing values
# era5_test_grid1.nc
ds2 = xr.open_dataset(f'{workdir}/snodas_gridtest_v1.nc', mask_and_scale=False)

# %%
ds2

# %%
# lat = ds['latitude'].values
# long = ds['longitude'].values
# elevation_band = ds['elevation_band'].values

# mean_elev = np.array([0.1, 0.5, 0.3, 0.6]).reshape((4, 1, 1))

# me = xr.DataArray(mean_elev, coords={'latitude': lat, 'longitude': long, 
#                                 'elevation_band': elevation_band},
#              dims=['elevation_band', 'latitude', 'longitude'])
# ds['mean_elev'] = me

# grid1 = np.array()
ds2['lsm'].values.shape

# %%
# grid1 = np.arange(1, ds2['lsm'].values.size+1).reshape((ds2['lsm'].values.shape), order='F')
grid1 = np.arange(1, ds2['LANDMASK'].values.size+1).reshape((ds2['LANDMASK'].values.shape), order='F')

# %%
grid1

# %%
# For gridmet
# gridxr = xr.DataArray(grid1, coords={'time': ds2['time'].values,
#                                      'latitude': ds2['latitude'].values,
#                                      'longitude': ds2['longitude'].values},
#                       dims=['time', 'latitude', 'longitude'])
# ds2['lsm'] = gridxr

# For SNODAS
gridxr = xr.DataArray(grid1, coords={'lat': ds2['lat'].values,
                                     'lon': ds2['lon'].values},
                      dims=['lat', 'lon'])
ds2['LANDMASK'] = gridxr

# %%
ds2

# %%
ds2['LANDMASK'].plot()

# %%
# Dataset.to_netcdf(self, path=None, mode: str = 'w', format: str = None, group: str = None, engine: str = None, encoding: Mapping = None, unlimited_dims: Iterable[Hashable] = None, compute: bool = True, invalid_netcdf: bool = False) → Union[bytes, ForwardRef('Delayed'), NoneType]
ds2.to_netcdf(f'{workdir}/gridtest_v2.nc')

# %%
