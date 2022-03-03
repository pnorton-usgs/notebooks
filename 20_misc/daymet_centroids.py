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

# %%

# %%
import numpy as np
import pandas as pd
import netCDF4 as cf
import kdtree_fast as kd

# %%
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v4'

daymet_file = f'{work_dir}/daymet_v4_landmask.nc'

# %%
fhdl = cf.Dataset(daymet_file, 'r')

# %%
lats = fhdl.variables['lat'][...]
lons = fhdl.variables['lon'][...]

# %%
lats

# %%
# %%time
# Load lats/lons into a kdtree
ns = kd.Kdtree_fast(fhdl, 'lat', 'lon')

# %%

# %%

# %%

# %%

# %%

# %%
