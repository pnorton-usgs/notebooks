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
#     display_name: Python [conda env:bandit_py3]
#     language: python
#     name: conda-env-bandit_py3-py
# ---

# %%
import glob
import os

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v3_1980-01-01_2016-12-31/netcdf/*.nc'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v3_1980-01-01_2016-12-31'

it = glob.iglob(workdir)

# %%
cfile = next(it)

# %%
os.path.basename(cfile).split('_')

# %%
