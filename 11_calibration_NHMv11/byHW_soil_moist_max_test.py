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
# import glob
import numpy as np
# import os
import pandas as pd
# import xarray as xr

# %%
headwater = '0259'
hw_suffix = ''
workdir = f'/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHW_sample/HW{headwater}{hw_suffix}'
prcp_file = f'{workdir}/prcp.cbh'

st_date = datetime.datetime(2000, 1, 1)
en_date = datetime.datetime(2010, 12, 31)

# %%
df = pd.read_csv(prcp_file, sep='\s+', skiprows=3, header=None, parse_dates={'time': [0, 1, 2]}, index_col='time',
                 usecols=[0, 1, 2, 6, 7, 8, 9])

# %%
df.head()

# %%
ppt_max = df.max().tolist()

# %%
soil_moist_max = [2.955417, 2.754325, 2.955417, 3.97341]

# %%
for xx, yy in zip(ppt_max, soil_moist_max):
    smidx_max = (1.1 * yy) + (0.5 * xx)
    print(smidx_max)

# %%
# smidx_max(chru) = (1.1 * soil_moist_max(chru)) + (0.5 * ppt_max(chru))

