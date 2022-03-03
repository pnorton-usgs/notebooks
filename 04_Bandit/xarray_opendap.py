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
#     display_name: Python [conda env:bandit_38]
#     language: python
#     name: conda-env-bandit_38-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import datetime
import os
import xarray as xr
import pandas as pd

# %%
base_url = 'http://gdp-netcdfdev.cr.usgs.gov:8080/thredds/dodsC/NHM_NWIS'
# base_url = 'http://gdp-netcdfdev.cr.usgs.gov:8080/thredds/fileServer/NHM_NWIS/files'
# base_url = 'http://gdp-netcdfdev.cr.usgs.gov:8080/thredds/catalog/NHM_NWIS/files'

# %%
# files = [f'{base_url}/NWIS_pois.nc', f'{base_url}/HYDAT_pois.nc']

# %%
ds = xr.open_dataset(base_url, decode_cf=True)

# %%
ds

# %%
ds['discharge'].coords

# %%
st_dt = datetime.datetime(1980,10,1)
en_dt = datetime.datetime(1990,9,30)

ds['discharge'].loc[st_dt:en_dt, '01018000'].plot()

# %%

# %%

# %%

# %%

# %%
base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/poi_data'
filename = f'{base_dir}/NWIS_pois.nc'

ds_local = xr.open_dataset(filename, decode_cf=True)

# %%
ds_local

# %%
ds_local['discharge'].loc['01018000', st_dt:en_dt].plot()

# %%
