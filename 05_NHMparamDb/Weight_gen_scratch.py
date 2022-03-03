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

import xarray as xr

from shapely.geometry import Polygon
import multiprocessing as mp
import psutil
import os
import resource
import math
import time

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/SNODAS/unmasked'
src_netcdf = f'{workdir}/SWE_subset.nc'
src_gridspacing = 0.04166666


# %%
def new_calc_dx_dy(longitude,latitude,shape,radius=6370997.):
    ''' This definition calculates the distance 
        between grid points that are in
        a latitude/longitude format.
        
        Using pyproj GEOD; different Earth Shapes 
        https://jswhit.github.io/pyproj/pyproj.Geod-class.html
        Common shapes: 'sphere', 'WGS84', 'GRS80'
        
        Accepts, 1D arrays for latitude and longitude
        
        Returns: dx, dy; 2D arrays of distances 
                       between grid points in the x and y direction in meters 
    '''
    # from: https://github.com/Unidata/MetPy/issues/288
    from pyproj import Geod
    
    if (radius != 6370997.):
        g = Geod(a=radius,b=radius)
    else:
        g = Geod(ellps=shape)
    
    dx = np.empty(latitude.shape)
    dy = np.zeros(longitude.shape)
    
    for i in range(latitude.shape[1]):
        for j in range(latitude.shape[0]-1):
            _, _, dx[j,i] = g.inv(longitude[j,i],latitude[j,i],longitude[j+1,i],latitude[j+1,i])
    dx[j+1,:] = dx[j,:]
    
    for i in range(latitude.shape[1]-1):
        for j in range(latitude.shape[0]):
            _, _, dy[j,i] = g.inv(longitude[j,i],latitude[j,i],longitude[j,i+1],latitude[j,i+1])
    dy[:,i+1] = dy[:,i]
    
    return dx, dy



# %%
lons = np.array([-124.730225, -124.72317])
lats = np.array([24.878002, 24.884838])

new_calc_dx_dy(lons, lats, 'WGS84')


# %%
def convert_size(size_bytes):
    # from https://python-forum.io/Thread-Convert-file-sizes-will-this-produce-accurate-results
    if size_bytes == 0:
        return "0 B"
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    power = math.pow(1024, i)
    size = round(size_bytes / power, 2)
    return f'{size} {size_name[i]}'


# %%
def create_polygons(lats, lons, grid_spacing, data):
    t1 = time.time()
    lon, lat = np.meshgrid(lons, lats)
    print(f'meshgrid: {time.time() - t1}')
    
    # For weight generation the datahandle variable is not used for anything
    # but creating the geodataframe of the source netcdf file
    t1 = time.time()
    df = pd.DataFrame({'grid': data})
    print(f'DataFrame: {time.time() - t1}')
    
    # res is half of the 'resolution' (e.g. gridspacing in degrees)
    res = grid_spacing / 2.0
    poly = []
    index = []
    count = 0

    # Create polygon features of the grid-cell bounding boxes
    t1 = time.time()
    for i in range(np.shape(lat)[0]):
        for j in range(np.shape(lon)[1]):
            lat_point_list = [lat[i, j] - res, lat[i, j] + res, lat[i, j] + res, lat[i, j] - res]
            lon_point_list = [lon[i, j] + res, lon[i, j] + res, lon[i, j] - res, lon[i, j] - res]
            poly.append(Polygon(zip(lon_point_list, lat_point_list)))
            index.append(count)
            count += 1
    print(f'poly: {time.time() - t1}')
    
    t1 = time.time()
    ncfcells = gpd.GeoDataFrame(df, index=index, geometry=poly, crs='epsg:4326')
    print(f'GeoDataFrame: {time.time() - t1}')

    return ncfcells


# %%
ds = xr.open_dataset(src_netcdf, mask_and_scale=False)

lathandle = ds['lat']
lonhandle = ds['lon']
datahandle = ds['SWE']

# Print some information on the data
print('\n Data sizes are: \n', datahandle.sizes)
print('\n Data coords are: \n', datahandle.coords)



# %%

# %%
print('Before create ncfcells mem usage:', convert_size(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
ncfcells = create_polygons(lathandle, lonhandle, src_gridspacing, datahandle.isel(time=0).values.flatten())
print('After create ncfcells mem usage:', convert_size(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

# %%
lon, lat = np.meshgrid(lonhandle, lathandle)
print(lon)

# %%
print(lon.shape)

# %%
res = src_gridspacing / 2.0

# %%
# i => lat, j => lon; [i, j]
i = 0
j = 0
lat_point_list = [lat[i, j] - res, lat[i, j] + res, lat[i, j] + res, lat[i, j] - res]
lon_point_list = [lon[i, j] + res, lon[i, j] + res, lon[i, j] - res, lon[i, j] - res]
# poly.append(Polygon(zip(lon_point_list, lat_point_list)))

# %%
print(f'{lat[i, j], lon[i, j]}')
print(lat_point_list)
print(lon_point_list)

# %%
ll_2d = list(zip(lon_point_list, lat_point_list))
print(ll_2d)

# %%
cc = lat[:, :] - res
dd = lat[:, :] + res

# np.array([[2],[3],[4]])
yy = np.dstack((cc, dd, dd, cc)).flatten()
# xx = np.dstack([cc, dd])

print(yy.shape)

# %%
# array([  24.85716797,   24.89883463,   24.89883463,   24.85716797,
#        -119.99688528, -119.99688528, -120.03855194, -120.03855194])

#24.85716797
#24.85716796875

#24.89883463
#24.89883462875

#-119.99688528
#-119.99688527731261

# -120.03855194
# -120.03855193731262

# %%
bb = np.ndarray(lat[:,:] - res, lat[:, :] + res, lat[:, :] + res, lat[:, :] - res)

# %%
ff = lon[:, :] + res
gg = lon[:, :] - res

xx = np.dstack((ff, ff, gg, gg)).flatten()
print(xx.shape)

# %%
cc[i, j]

# %%
# np.array([[2],[3],[4]])
yy = np.dstack((cc, dd, dd, cc)).flatten()
# xx = np.dstack([cc, dd])

print(yy.shape)

# %%

# %%
t1 = time.time()
for i in range(np.shape(lat)[0]):
    for j in range(np.shape(lon)[1]):
        lat_point_list = [lat[i, j] - res, lat[i, j] + res, lat[i, j] + res, lat[i, j] - res]
        lon_point_list = [lon[i, j] + res, lon[i, j] + res, lon[i, j] - res, lon[i, j] - res]
print(f'poly: {time.time() - t1}')

# %%
zz = np.dstack((yy, xx))

# %%
zz[0:2,0,:]

# %%
zz.shape


# %%
def create_polygons_v2(lats, lons, grid_spacing, data):
    
    def get_poly(row, lats, lons):
        ii = row.name * 4
        return Polygon(zip(lats[ii:ii+4], lons[ii:ii+4]))
    
    t1 = time.time()
    lon, lat = np.meshgrid(lons, lats)
    print(f'meshgrid: {time.time() - t1}')

    # res is half of the 'resolution' (e.g. gridspacing in degrees)
    res = grid_spacing / 2.0
    
    # For weight generation the datahandle variable is not used for anything
    # but creating the geodataframe of the source netcdf file
    t1 = time.time()
    df = pd.DataFrame({'grid': data})
    print(f'DataFrame: {time.time() - t1}')

    # Create polygon features of the grid-cell bounding boxes
    t1 = time.time()    
    lats1 = lat[:, :] - res
    lats2 = lat[:, :] + res
    
    yy = np.dstack((lats1, lats2, lats2, lats1)).flatten()
    
    lons1 = lon[:, :] + res
    lons2 = lon[:, :] - res

    xx = np.dstack((lons1, lons1, lons2, lons2)).flatten()
    print(f'numpy: {time.time() - t1}')
    print('    mem usage:', convert_size(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    
    t1 = time.time()
    df['geometry'] = df.apply(get_poly, lats=yy, lons=xx, axis=1)
    print(f'df_geometry: {time.time() - t1}')
    print('    mem usage:', convert_size(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

    t1 = time.time()
    ncfcells = gpd.GeoDataFrame(df, crs='epsg:4326')
    print(f'GeoDataFrame: {time.time() - t1}')
    del df
    return ncfcells


# %%
print('Before create ncfcells mem usage:', convert_size(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
ncfcells = create_polygons_v2(lathandle, lonhandle, src_gridspacing, datahandle.isel(time=0).values.flatten())
print('After create ncfcells mem usage:', convert_size(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

# %%
ncfcells.info()

# %%
ncfcells.head()

# %%
print('After create ncfcells mem usage:', convert_size(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

# %%

# %%
