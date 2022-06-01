# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np

import pyproj
from pyproj import CRS
from pyproj import Transformer
from pyproj.enums import WktVersion

# %%
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/FY2021_projects/conus404/data'

filename = f'{work_dir}/landmask_6.nc'

# %%
xdf = xr.open_mfdataset(filename, decode_cf=True, engine='netcdf4')

# %%
xdf

# %%
xdf['LANDMASK'].plot()

# %%
t2 = xdf.LANDMASK
xx = xdf.west_east.values
yy = xdf.south_north.values

globe = ccrs.Globe(ellipse='sphere', semimajor_axis=6370000, semiminor_axis=6370000)
lcc = ccrs.LambertConformal(globe=globe, central_longitude=262.1, central_latitude=39.1,
                            standard_parallels=[30.0, 50.0])

# ax = plt.axes(projection=ccrs.Orthographic(-80, 35))
ax = plt.axes(projection=ccrs.PlateCarree())
# ax = plt.axes(projection=lcc)

# ax.plot(x, y, t2, transform=ccrs.PlateCarree())
# t2.plot(ax=ax, transform=ccrs.PlateCarree())
t2.plot.contourf(ax=ax, transform=lcc)
ax.coastlines()

ax.add_feature(cartopy.feature.BORDERS, linestyle='-');
ax.set_extent([xx.min(), xx.max(), yy.min(), yy.max()], crs=lcc)

# %%
xx.max()

# %%
xx.max()

# %%

# %%

# %%
# cen_lon = 262.1
cen_lon = -97.9
cen_lat = 39.1
dx = 4000
dy = 4000
lon_sz = 1367
lat_sz = 1015

# wrf_crs = CRS.from_proj4(f'+proj=latlong +a=6370000 +b=6370000 +pm=0.0')
# wrf_crs = CRS.from_proj4(f'+proj=lcc +lat_1=30.0 +lat_2=50.0 +lat_0={cen_lat} +lon_0={cen_lon}')
wrf_crs = CRS.from_proj4(f'+proj=lcc +lat_1=30.0 +lat_2=50.0 +lat_0={cen_lat} +lon_0={cen_lon} +a=6370000 +b=6370000')
# wrf_crs = CRS.from_proj4(f'+proj=lcc +lat_1=30.0 +lat_2=50.0 +lat_0={cen_lat} +lon_0={cen_lon} +a=6370000 +b=6370000 +towgs84=0,0,0 +no_defs')
print(wrf_crs.to_wkt(WktVersion.WKT1_GDAL, pretty=True))

wgs_crs = CRS.from_epsg(4326)

tform = Transformer.from_crs(wgs_crs, wrf_crs, always_xy=True)

e, n = tform.transform(cen_lon, cen_lat)

# Down left corner of the domain
x0 = -(lon_sz - 1) / 2.0 * dx + e
y0 = -(lat_sz - 1) / 2.0 * dy + n

# 2D grid
xx, yy = np.meshgrid(np.arange(lon_sz) * dx + x0, np.arange(lat_sz) * dy + y0)

# %%
wrf_crs = CRS.from_proj4(f'+proj=latlong +a=6370000 +b=6370000 +pm=0.0')
tform = Transformer.from_crs(wgs_crs, wrf_crs, always_xy=True)
e, n = tform.transform(cen_lon, cen_lat)
print(e,n)

# %%
wrf_crs = CRS.from_proj4(f'+proj=latlong +a=6370000 +b=6370000 +pm=0.0 +towgs84=0,0,0')
tform = Transformer.from_crs(wgs_crs, wrf_crs, always_xy=True)
e, n = tform.transform(cen_lon, cen_lat)
print(e,n)

# %%
print(x0)
print(y0)


# %%
print(f'e = {e}\nn = {n}')

# %%
tform.transform(cen_lon, cen_lat)

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
# wrf_crs = CRS.from_proj4(f'+proj=lcc +lat_1=30.0 +lat_2=50.0 +lat_0={cen_lat} +lon_0={cen_lon} +a=6370000 +b=6370000')

wrf_to_wgs = Transformer.from_crs(wrf_crs, wgs_crs, always_xy=True)

new_lon, new_lat = wrf_to_wgs.transform(xx, yy)

xdf['diff'] = np.sqrt((new_lon - xdf.XLONG)**2 + (new_lat - xdf.XLAT)**2)
# xdf['diff'] = np.sqrt((new_lon - xdf['XLONG'])**2 + (new_lat - xdf['XLAT'])**2)

print(f'The max diff is: {xdf["diff"].max().values}')

# %%
xdf['diff'].plot(cmap='Reds', figsize=(15,10))

# %%
# Plot the difference across longitudes for a given latitude
xdf['diff'].isel(south_north=[500]).plot()

# %%
print(xdf.XLAT.values[10, 10], new_lat[10, 10])
print(xdf.XLONG.values[10, 10], new_lon[10, 10])

# %%
print(xlat[10, 10], new_lat[10, 10])
print(xlon[10, 10], new_lon[10, 10])

# %%
# Bad projection when GCS is not specified as sphere
wrf_crs = CRS.from_proj4(f'+proj=lcc +lat_1=30.0 +lat_2=50.0 +lat_0={cen_lat} +lon_0={cen_lon}')
wrf_to_wgs = Transformer.from_crs(wrf_crs, wgs_crs, always_xy=True)

new_lon, new_lat = wrf_to_wgs.transform(xx, yy)

xdf['diff'] = np.sqrt((new_lon - xdf.XLONG)**2 + (new_lat - xdf.XLAT)**2)

print(f'The max diff is: {xdf["diff"].max().values}')
xdf['diff'].plot(cmap='Reds', figsize=(15,10))

# %%
# Plot the difference across longitudes for a given latitude
xdf['diff'].isel(south_north=[1014]).plot()

# %%

# %%

# %%
lat = xdf.XLAT.values
lon = xdf.XLONG.values

lcc = ccrs.LambertConformal(globe=globe, central_longitude=262.1, central_latitude=39.1,
                            standard_parallels=[30.0, 50.0])

# First, transform from LCC to Orthographic
transform = lcc.transform_points(ccrs.Orthographic(265,25), lon, lat)
x = transform[..., 0]
y = transform[..., 1]

ax = plt.axes(projection=ccrs.Orthographic(265,25))
ax.pcolormesh(x, y, t2, transform=ccrs.PlateCarree())
ax.add_feature(cf.NaturalEarthFeature(
               category='cultural',
               name='admin_1_states_provinces_lines',
               scale='50m',
               facecolor='none'))
ax.coastlines('50m')
ax.add_feature(cf.BORDERS)
ax.set_global()

# %%
xdf.XLAT.data

# %%

# %%
# from pyproj import Transformer
# transformer = Transformer.from_crs("epsg:4326", "epsg:3857")
# transformer.transform(12, 12)

# %%
from pyproj import CRS
from pyproj import Transformer
from pyproj.enums import WktVersion

src_crs = CRS.from_proj4('+proj=lcc +lat_1=30.0 +lat_2=50.0 +lat_0=39.1 +lon_0=262.1 +a=6370000 +b=6370000 +towgs84=0,0,0 +no_defs')
print(src_crs.to_wkt(WktVersion.WKT1_GDAL, pretty=True))
# tform = Transformer.from_crs()
# tform.transform()

# globe = ccrs.Globe(ellipse='sphere', semimajor_axis=6370000, semiminor_axis=6370000)
# lcc = ccrs.LambertConformal(globe=globe, central_longitude=262.1, central_latitude=39.1,
#                             standard_parallels=[30.0, 50.0])

# %%
from pprint import pprint

pprint(src_crs.to_cf())

# %%
dst_crs = CRS.from_epsg(4326)

tform = Transformer.from_crs(dst_crs, src_crs, always_xy=True)
tform

# %%
tform.transform(-97.1, 39.1)

# %%
import pyproj

# dst_proj = pyproj.Proj(f'+proj={wrf_proj} +lat_1={truelat1} +lat_2={truelat2} +lat_0={cen_lat} +lon_0={cen_lon} +a={spheroid_a} +b={spheroid_b} +towgs84=0,0,0 +no_defs')
wrf_proj = pyproj.Proj('+proj=lcc +lat_1=30.0 +lat_2=50.0 +lat_0=39.1 +lon_0=262.1 +a=6370000 +b=6370000 +towgs84=0,0,0 +no_defs')

wgs_proj = pyproj.Proj(proj='latlong', datum='WGS84')
e, n = pyproj.transform(wgs_proj, wrf_proj, 261.1, 39.1)

print(e,n)

# %%

# %%
