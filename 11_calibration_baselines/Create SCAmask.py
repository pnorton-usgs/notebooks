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
import pandas as pd
import xarray as xr
import numpy as np

# %%
baseline_var = 'SCA'
sca_var = f'{baseline_var}'
ci_var = f'{baseline_var}ci'

baseline_dir = f'/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/baselines/conus/{baseline_var}'

baseline_file = f'{baseline_dir}/baselines_{baseline_var}.nc'

# %%
st_date = datetime.datetime(2000,1,1)
en_date = datetime.datetime(2010,12,31)

# %%
# Baseline recharge has daily timestep
baseline_df = xr.open_dataset(baseline_file)
baseline_df = baseline_df[[sca_var, ci_var, 'hru']].sel(time=slice(st_date, en_date))
# baseline_df

# Remove July and August from the dataset
baseline_restr = baseline_df.sel(time=baseline_df.time.dt.month.isin([1, 2, 3, 4, 5, 6, 9, 10, 11, 12]))
baseline_restr

# %%
# Create the SCAmask

# Compute lower and upper SCA values based on confidence interval
threshold = 70.0
ci_pct = baseline_restr['SCAci'].where(baseline_restr['SCAci'] >= threshold)
ci_pct /= 100.0

sca_obs = baseline_restr['SCA'].where(~np.isnan(ci_pct))
msk_SCAmax = sca_obs.max(axis=0)
msk_num_obs = (sca_obs > 0.0).sum(axis=0)
msk_num_ann = sca_obs.resample(time='1AS').mean()
msk_num_ann = (msk_num_ann > 0).sum(axis=0)

SCAmask = (msk_num_ann > 0) & (msk_SCAmax > 0.5) & (msk_num_obs > 9)

# %%
# maxSCA[mhru]
# maxSCAyr[nyr, nhru]
# num_sca
# SCAmth

# %%
msk_num_ann[msk_num_ann <= 1]

# %%
msk_num_obs[msk_num_obs < 10]

# %%
msk_SCAmax[msk_SCAmax <= 0.50]

# %%
SCAmask

# %%
