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

# %%
print(gpd.__version__)

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/ERA5'
geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
# geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_wgs84/nhruNationalIdentifier.shp'

# %%
# Open the source netcdf file
ds = xr.open_dataset(f'{workdir}/land_sea_mask_upk.nc', mask_and_scale=False)
print(ds)

# %%
# Read the geodatabase
layer_name = 'nhru_v11'
gdf = gpd.read_file(geodatabase, layer=layer_name)
gdf.head()

# %%
# gdf.crs = 'EPSG:5070'
gdf.crs

# %%
# Re-project to WGS84 lat/lon
gdf = gdf.to_crs(epsg=4326)
gdf.crs

# %%
# Create some variables from the netcdf file for use later
print('The meta data is: \n', json.dumps(ds.attrs, indent=4))

lathandle = ds['latitude']
lonhandle = ds['longitude']
timehandle = ds['time']
datahandle = ds['lsm']

# %%
datahandle.isel(time=0).plot()

# %%
type(lonhandle)

# %%
lon, lat = np.meshgrid(lonhandle, lathandle)
df = pd.DataFrame({'grid': datahandle.isel(time=0).values.flatten()})

# res is the 'resolution' (e.g. gridspacing in degrees)
gridspacing = 0.25
res = gridspacing / 2.0
numcells = lon.size
poly = []
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
ncfcells = gpd.GeoDataFrame(df, index=index, geometry=poly, crs='epsg:4326')
ncfcells.head()

# %%
lon


# %% [markdown]
# ## Begin parallel experiment

# %%

# %%
# This is base on https://swanlund.space/parallelizing-python

def outer_loop(shp_chunk, src_netcdf):
    id_key = 'nhru_v11'
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
                    results.append([np.int(precise_matches.index[count]), np.int(irow[id_key]), tmpfloat])
                    count += 1
    return results

def parallelize():
    cpus = 6
    # cpus = mp.cpu_count()

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
id_key = 'nhru_v11'

# Write the results to a weights file
with open(f'{workdir}/tmp_weights2.csv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(['grid_ids', id_key, 'w'])
    
    for pp in final_res:
        for xx in pp:
            for dd in xx:
                writer.writerow([dd[0], dd[1], dd[2]])


# %%
