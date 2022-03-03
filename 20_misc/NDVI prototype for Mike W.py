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
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import datetime
import xarray as xr
import pandas as pd
import numpy as np

# %%
work_dir = '/Users/pnorton/tmp/NDVI'

# %%
from dask.distributed import Client
client = Client()
client

# %%
xdf = xr.open_mfdataset(f'{work_dir}/MCD13.A2007.unaccum.nc4', engine='netcdf4')

# %%
xdf

# %%
st_date = datetime.datetime(2007, 1, 1)
en_date = datetime.datetime(2007, 1, 31)

# %%
aa = xdf['NDVI'].loc[st_date:en_date, :, :].resample(time='1M').mean()

# %%
aa

# %%
aa.plot()

# %%
aa.max()

# %%
