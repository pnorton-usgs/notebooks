# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python [conda env:holoview]
#     language: python
#     name: conda-env-holoview-py
# ---

# %%
import xarray as xr

# import cartopy
# import cartopy.crs as ccrs
# import cartopy.feature as cf
import numpy as np

import pyproj
from pyproj import CRS
from pyproj import Transformer
from pyproj.enums import WktVersion


# %%
base_dir = '/Volumes/USGS_NHM2/datasets/wrf-conus404/originals'
filename = f'{base_dir}/wrfout_d01_2020-09-30_00:00:00'

# %%

# %%
xdf = xr.open_mfdataset(filename, decode_cf=True, engine='netcdf4')

# %%
xdf

# %%
xdf.TRUELAT1

# %%
xdf.LANDMASK.plot()

# %%
# %%time
cen_lon = xdf.CEN_LON   # -97.9
cen_lat = xdf.CEN_LAT   # 39.1
dx = xdf.DX   # 4000
dy = xdf.DY   # 4000
lat1 = xdf.TRUELAT1
lat2 = xdf.TRUELAT2
lon_sz = xdf.dims['west_east']    # 1367
lat_sz = xdf.dims['south_north']   # 1015

wrf_crs = CRS.from_proj4(f'+proj=lcc +lat_1={lat1} +lat_2={lat2} +lat_0={cen_lat} +lon_0={cen_lon} +a=6370000 +b=6370000')
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
xx

# %%

# %%

# %%

# %%

# %%

# %%
