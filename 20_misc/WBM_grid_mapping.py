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
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
# 2022-01-10 from Mike W.
"""Example notebook to work with nldas soil moisture netcdf files."""

import geopandas as gpd
import numpy as np
import time
import pickle
import xarray as xr
import cftime
import datetime
from datetime import timedelta
import rioxarray
import pandas as pd
from pyproj import CRS

# %%
crs = CRS("epsg:4326")
cf_grid_mapping = crs.to_cf()
cf_grid_mapping

# %%

# %%
ds = xr.open_dataset("F:/ClimGrid/ClimGrid_WBM.nc")

# %%
ds

# %%

# %%
ds.rio.write_crs('epsg:4326', inplace=True)

# %%
ds.aet.attrs['grid_mapping'] = 'spatial_ref'
ds.pet.attrs['grid_mapping'] = 'spatial_ref'
ds.prcp.attrs['grid_mapping'] = 'spatial_ref'
ds.rain.attrs['grid_mapping'] = 'spatial_ref'
ds.runoff.attrs['grid_mapping'] = 'spatial_ref'
ds.snow.attrs['grid_mapping'] = 'spatial_ref'
ds.soilstorage.attrs['grid_mapping'] = 'spatial_ref'
ds.swe.attrs['grid_mapping'] = 'spatial_ref'
ds.tmean.attrs['grid_mapping'] = 'spatial_ref'

# %%
ds

# %%
ds.to_netcdf("F:/ClimGrid/ClimGrid_WBModel.nc", mode='a')

# %%
