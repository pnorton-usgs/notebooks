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
import numpy as np
import xarray as xr
import matplotlib

# %%
# ERA5
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/ERA5'

# SNODAS
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/SNODAS/unmasked'
# src_netcdf = f'{workdir}/snodas_landmask_NHM.nc'

# DAYMET v3
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v3'
# src_netcdf = f'{workdir}/daymet_v3_landmask_proj.nc'

# DAYMET v4
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v4'
src_netcdf = f'{workdir}/daymet_v4_landmask.nc'

dst_netcdf = f'{workdir}/daymet_v4_landmask_w_patterns.nc'

# %%
# Trying to create a simple test grid with incrementing values
ds2 = xr.open_dataset(src_netcdf, mask_and_scale=False, chunks={})
ds2

# %%
# Create a grid array the same size and shape as the source netcdf variable.
# grid array is just a monotonic series
# grid1 = np.arange(1, ds2['lsm'].values.size+1).reshape((ds2['lsm'].values.shape), order='F')
grid1 = np.arange(1, ds2['LANDMASK'].values.size+1).reshape((ds2['LANDMASK'].values.shape), order='F')
grid_v = np.arange(1, ds2['LANDMASK'].values.size+1).reshape((ds2['LANDMASK'].values.shape), order='F')
grid_h = np.arange(1, ds2['LANDMASK'].values.size+1).reshape((ds2['LANDMASK'].values.shape), order='C')

# %%
# Convert the grid array to an xarray DataArray with the same dimensions as the source netcdf file
# This is kind of lazy because this is just overwriting another variable

# For gridmet
# gridxr = xr.DataArray(grid1, coords={'time': ds2['time'].values,
#                                      'latitude': ds2['latitude'].values,
#                                      'longitude': ds2['longitude'].values},
#                       dims=['time', 'latitude', 'longitude'])
# ds2['lsm'] = gridxr

# xr.merge([ds, ds.rename({'foo': 'bar'})])

# For SNODAS
# gr_vert = xr.DataArray(grid_v, coords={'lat': ds2['lat'].values,
#                                        'lon': ds2['lon'].values},
#                        dims=['lat', 'lon'])

# gr_horiz = xr.DataArray(grid_h, coords={'lat': ds2['lat'].values,
#                                         'lon': ds2['lon'].values},
#                         dims=['lat', 'lon'])

# For Daymet V3
# gr_vert = xr.DataArray(grid_v, coords={'y': ds2['y'].values,
#                                        'x': ds2['x'].values},
#                        dims=['y', 'x'])

# gr_horiz = xr.DataArray(grid_h, coords={'y': ds2['y'].values,
#                                         'x': ds2['x'].values},
#                         dims=['y', 'x'])

# gridxr = xr.DataArray(grid1, dims=['y', 'x'])

# For Daymet V4
gr_vert = xr.DataArray(grid_v, coords={'y': ds2['y'].values, 'x': ds2['x'].values}, dims=['y', 'x'])
gr_horiz = xr.DataArray(grid_h, coords={'y': ds2['y'].values, 'x': ds2['x'].values}, dims=['y', 'x'])
gridxr = xr.DataArray(grid1, dims=['y', 'x'])

ds3 = xr.merge([ds2, xr.Dataset({'grid_vert': gr_vert})])
ds3 = xr.merge([ds3, xr.Dataset({'grid_horiz': gr_horiz})])
# ds2['grid_vert'] = gridxr
# ds2['LANDMASK'] = gridxr

# %%
ds3

# %%

# %%
ds3['grid_horiz'].plot()

# %%
# Write the updated netcdf to a new file

# For SNODAS
# ds3.to_netcdf(dst_netcdf, encoding={'lat': {'_FillValue': None}, 'lon': {'_FillValue': None},})

# For Daymet V3
# ds3.to_netcdf(dst_netcdf, encoding={'x': {'_FillValue': None}, 'y': {'_FillValue': None},
#                                     'lat': {'_FillValue': None}, 'lon': {'_FillValue': None},})

# For Daymet V4
encoding = {}
encoding['x'] = dict(zlib=True, complevel=2, fletcher32=False, shuffle=True, _FillValue=None)
encoding['y'] = dict(zlib=True, complevel=2, fletcher32=False, shuffle=True, _FillValue=None)
encoding['lat'] = dict(zlib=True, complevel=2, fletcher32=False, shuffle=True, _FillValue=None)
encoding['lon'] = dict(zlib=True, complevel=2, fletcher32=False, shuffle=True, _FillValue=None)
encoding['LANDMASK'] = dict(zlib=True, complevel=2, fletcher32=False, shuffle=True)

ds3.to_netcdf(dst_netcdf, encoding=encoding)

# %%
