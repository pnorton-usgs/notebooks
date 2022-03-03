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
#     display_name: Python [conda env:bandit]
#     language: python
#     name: conda-env-bandit-py
# ---

# %%

# %%
import matplotlib.pyplot as plt
import matplotlib as mpl
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

import xarray
from dask.distributed import Client

# import netCDF4 as cf
import pandas as pd
from datetime import datetime
from calendar import monthrange
# from collections import namedtuple
import numpy as np
# import os

# %%
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/conus/output'
filename = '{}/crap.nc'.format(work_dir)

# Open the netcdf file of NHM output variables
xdf = xarray.open_mfdataset(filename, chunks={'nhru': 1040}, combine='by_coords')
xdf = xdf.set_index(nhru='nhm_id')

# %%
my_ids = [46465]
ds_pktemp = xdf['pk_temp'].sel(nhru=my_ids)
print(f'{ds_pktemp.min().values:f}; {ds_pktemp.max().values:f}')

# %%
print(ds_pktemp)

# %%
ds_pktemp.plot()

# %%
ds_pktemp.max().values

# %%
