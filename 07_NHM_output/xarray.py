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
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/byHRU_musk/output/netcdf'
filename = f'{workdir}/crap_125_1044.nc'

# %%
ds = xr.open_dataset(filename)

# %%
# df = ds['hru_actet'].to_dataframe()

# nc3 [-, -]: 17.8G high-water; 13.7G after conversion
# [177, 1476]: 13.7G after read; compressed file size is 4.2G; hru_actet uncompressed is 5.4G
# [1, 109951]: 19G high-water; 13.7G after conversion

# I think this problem is related to: https://github.com/pydata/xarray/issues/2534

# %%
# df.info()

# nc3 [-, -]: killed at memory of 44G
# nc4 [177, 1476]: Memory usage over 55G; killed python process

# %%
# %%time

# Single HRU, all timesteps
ds['hru_actet'][:, 0].plot()

# %%
# %%time

# All HRUs, single timestep
ds['hru_actet'][30, :].plot()

# %%
# %%time

# Compute monthly mean values
mean_month = ds['hru_actet'].resample(time='1M').mean()

# mean_month = albedo.Albedo.resample(freq='M',dim='time')

# %%
mean_month

# %%
# %%time

# single timestep, all HRUs
mean_month[:, 101].plot()

# %%
# %%time

# single HRU, all timesteps
mean_month[105, :].plot()

# %%
