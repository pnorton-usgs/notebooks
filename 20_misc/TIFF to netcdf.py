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
#     display_name: Python [conda env:holoview]
#     language: python
#     name: conda-env-holoview-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
# import hvplot.xarray
import netCDF4 as cf
import numpy as np
import rasterio
from rasterio.warp import transform
import warnings
import xarray as xr

warnings.filterwarnings('ignore')

# %% [markdown]
# This script converts TIFF-format monthly climatology (12 files) to netcdf, adding appropriate metadata.
#
# 30-Year (1990-2019) Monthly Averages of DAYMET Precipitation and Temperature for North America
# https://www.sciencebase.gov/catalog/item/620a8120d34ec05caca60bb7

# %%
cvar = 'prcp'   # One of prcp, tavg

# Location of TIFF files
work_dir = '/Users/pnorton/tmp/20220214_mike'

# Path to output netcdf file
output_file = f'daymet_{cvar}_1990-2019.nc'

months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# Build list of files to process
filelist = [f'{work_dir}/{cvar}_{mm}.tif' for mm in months]

# %% [markdown]
# ## Create time and time bounds

# %%
st = xr.cftime_range(start='2004', periods=12, freq='MS')
bnd_st = xr.cftime_range(start='1990', periods=12, freq='MS')
bnd_en = xr.cftime_range(start='2019', periods=12, freq='M').shift(1, "D")

time_units = 'days since 1990-01-01'
time_calendar = 'standard'

tb = np.zeros((12, 2), dtype=int)

tb[:, 0] = cf.date2num(bnd_st, units=time_units, calendar=time_calendar)
tb[:, 1] = cf.date2num(bnd_en, units=time_units, calendar=time_calendar)

time_var = xr.Variable('time', cf.date2num(st, units=time_units, calendar=time_calendar),
                       attrs=dict(units=time_units, calendar=time_calendar))
bnd_var = xr.Variable(['time', 'bnds'], tb)

# %% [markdown]
# ## Create data array of the tiff files

# %%
# %%time
# Create data array of the individual tiff files
da = xr.concat([xr.open_rasterio(ii) for ii in filelist], dim=time_var)

# %% [markdown]
# ## Create 2D latitude and longitude arrays from projection coordinates

# %%
# %%time
# Compute the lon/lat coordinates with rasterio.warp.transform
ny, nx = len(da['y']), len(da['x'])
x, y = np.meshgrid(da['x'], da['y'])

# Rasterio works with 1D arrays
lon, lat = transform(da.crs, {'init': 'EPSG:4326'},
                     x.flatten(), y.flatten())
lon = np.asarray(lon).reshape((ny, nx))
lat = np.asarray(lat).reshape((ny, nx))

# Make sure longitude values are in range: 0 to 360
lon[np.where(lon < 0.0)] += 360.0

da.coords['lon'] = (('y', 'x'), lon)
da.coords['lat'] = (('y', 'x'), lat)

# %% [markdown]
# ## Create dataset from the data array

# %%
# %%time
# Create an xarray dataset from the data array, rename the 'band' variable to 'prcp'
ds = da.to_dataset('band')
ds = ds.rename({1: cvar})

# This really isn't needed
# ds = ds.fillna(9.96921e+36)

ds['climatology_bounds'] = bnd_var

# %% [markdown]
# ## Define the coordinate reference system

# %%
# %%time
mycrs = rasterio.crs.CRS.from_string(da.crs, morph_from_esri_dialect=True)
crs_parts = mycrs.to_dict()

ds['crs'] = 0
ds.crs.attrs = {'grid_mapping_name': 'lambert_conformal_conic',
                'crs_wkt': mycrs.to_wkt(),
                'standard_parallel': f'{crs_parts["lat_1"]}, {crs_parts["lat_2"]}',
                'longitude_of_central_meridian': crs_parts['lon_0'],
                'latitude_of_projection_origin': crs_parts['lat_0'],
                'false_easting': crs_parts['x_0'],
                'false_northing': crs_parts['y_0']}

# %% [markdown]
# ## Set attributes for the variables

# %%
ds.x.attrs = {'units': 'm', 'long_name': 'x coordinate of projection', 'standard_name': 'projection_x_coordinate'}
ds.y.attrs = {'units': 'm', 'long_name': 'y coordinate of projection', 'standard_name': 'projection_y_coordinate'}

ds.lon.attrs = {'units': 'degrees_east', 'long_name': 'longitude coordinate', 'standard_name': 'longitude'}
ds.lat.attrs = {'units': 'degrees_north', 'long_name': 'latitude coordinate', 'standard_name': 'latitude'}

ds.time.attrs = dict(bounds='climatology_bounds', units=time_units, calendar=time_calendar, 
                     standard_name='time')

if cvar == 'prcp':
    ds[cvar].attrs = {'units': 'mm', 'grid_mapping': 'crs', 'coordinates': 'lat lon', 'cell_methods':'time: sum within year time: mean over years'}
    
if cvar == 'tavg':
    ds[cvar].attrs = {'units': 'degC', 'grid_mapping': 'crs', 'coordinates': 'lat lon', 'cell_methods':'time: mean within year time: mean over years'}

# %% [markdown]
# ## Add/remove global attributes

# %%
# Remove the existing global attributes
rmlist = [key for key in ds.attrs.keys()]
for aa in rmlist:
    del ds.attrs[aa]

# %%
# Add global attributes
str_var = ''
if cvar == 'prcp':
    str_var = 'precipitation'
elif cvar == 'tavg':
    str_var = 'temperature'
    
global_attrs = dict(title=f'30-year (1990-2019) Monthly Averages of DAYMET {str_var} for North America',
                    Conventions='CF-1.7',
                    keywords=f'North America, {str_var}',
                    source='Daymet: Annual Climate Summaries on a 1-km Grid for North America, Version 4, DOI: https://doi.org/10.3334/ornldaac/1852',
                    summary=f'These data represent the 30-year average monthly {str_var} for North America for the period 01-01-1990 to 12-31-2019')

for kk, vv in global_attrs.items():
    ds.attrs[kk] = vv
    

# %%
ds

# %% [markdown]
# ## Create a netcdf file from the dataset

# %%
# %%time
encoding = {'x': dict(zlib=True, complevel=2, fletcher32=False, shuffle=True, _FillValue=None, dtype=np.float32),
            'y': dict(zlib=True, complevel=2, fletcher32=False, shuffle=True, _FillValue=None, dtype=np.float32),
            'lon': dict(zlib=True, complevel=2, fletcher32=False, shuffle=True, _FillValue=None, dtype=np.float32),
            'lat': dict(zlib=True, complevel=2, fletcher32=False, shuffle=True, _FillValue=None, dtype=np.float32),
            'time': dict(zlib=True, complevel=2, fletcher32=False, shuffle=True, _FillValue=None),
            'climatology_bounds': dict(zlib=True, complevel=2, fletcher32=False, shuffle=True, _FillValue=None),
            cvar: dict(zlib=True, complevel=2, fletcher32=False, shuffle=True, _FillValue=9.96921e+36, dtype=np.float32)}

ds.to_netcdf(output_file, engine='netcdf4', encoding=encoding)

# %% [markdown]
# ## Sample plot of data for March using lat/lon

# %%

# %%
# %%time
ds[cvar].isel(time=2).hvplot.quadmesh(x='lon', y='lat', rasterize=True, cmap='turbo')  # .redim.nodata(value=9.96921e+36)
# da.hvplot.quadmesh(x='x', y='y', rasterize=True, geo=True, tiles='OSM', alpha=0.7, cmap='turbo')

# %% [markdown]
# ## Sample plot of data for March using projection coordinates

# %%
ds[cvar].isel(time=2).plot()

# %%
