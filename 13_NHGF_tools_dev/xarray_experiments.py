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
from dask.distributed import Client
# import geopandas
import numpy as np
# import pandas as pd
# import os
import xarray

# %%
client = Client(n_workers=2, threads_per_worker=2, memory_limit='1GB')
# client = Client(processes=False)
client

# %%
client.close()

# %%
work_dir = '/Volumes/USGS_NHM2/datasets/wrf-conus404/originals/WY2020'

# %%
xdf = xarray.open_mfdataset(f'{work_dir}/wrf2d*', concat_dim='Time', combine='nested',
                         parallel=True, coords="minimal", data_vars="minimal", 
                         compat='override', chunks={}) #  , combine='by_coords')

# xdf = xdf.set_index(hruid='hruid')
xdf

# %%
for kk,vv in xdf.dims.items():
    print(kk,vv)

# %%
xdf.dims

# %%
xdf['XTIME']

# %%

# %%

# %%

# %%

# %%
cattrs = xdf['T2'].attrs.keys()
cattrs

# %%
var_attrs_overrides = {'time': {'axis': 'T', 'standard_name': 'time', 'calendar': 'standard'},
                       'T2': {'description': 'Temperature at 2 meters', 'standard_name': 'air_temperature'},
                       'PSFC': {'description': 'Surface air pressure', 'standard_name': 'surface_air_pressure'}}
skip_attrs = ['FieldType', 'MemoryOrder', 'stagger', 'cell_methods']

# %%
cv_or = var_attrs_overrides['T2'].keys()

# %%
# Remove any attributes that should be skipped
wr_attrs = set(cattrs).difference(set(skip_attrs))

# %%
wr_attrs = wr_attrs.union(set(cv_or))

# %%
xdf.T2.values.shape

# %%
xdf.T2.values.ndim

# %%
list(xdf.keys())

# %%

# %%

# %%

# %%
