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

# %%

# %%
import dask
import datetime
import fsspec
import pandas as pd
import time
import xarray as xr

from dask.distributed import Client

# %%
work_dir = '/Volumes/USGS_NHM2/datasets/wrf-conus404/originals'
filename = f'{work_dir}/wrfout_d01_2020-09-30_00:00:00'

# %%
client = Client()

# %%
client

# %%
df = xr.open_dataset(filename, chunks={})

# %%
df

# %%
for vv in df.variables:
    cvar = vv
    try:
        minval = df[cvar].min().values
        maxval = df[cvar].max().values
        meanval = df[cvar].mean().values
        print(f'{cvar}: {minval}, {maxval}, {meanval}')
    except TypeError:
        print(f'{cvar}: TypeError')

# %%
df.variables.keys()

# %%

# %%
thresh = 1e9
acswdnb_0 = 8.402427e+08
i_acswdnb_0 = 290.0

acswdnb_1 = 8.420878e+08
i_acswdnb_1 = 290.0

swdnb = (acswdnb_1 + thresh * i_acswdnb_1) - (acswdnb_0 + thresh * i_acswdnb_0)
print(swdnb)

# %%
print((acswdnb_1 + thresh * i_acswdnb_1) / 3600.)

# %%
swdnb / 3600.

# %%

# %%

# %%

# %%

# %%
