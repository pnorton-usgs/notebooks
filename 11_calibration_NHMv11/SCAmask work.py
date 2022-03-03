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
import glob
import numpy as np
import os
import pandas as pd
import xarray as xr

# %%
headwater = '0259'
hw_suffix = ''
workdir = f'/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHW_sample/HW{headwater}{hw_suffix}/SCA'
ofs_file = f'{workdir}/objfun_{headwater}'

st_date = datetime.datetime(2000, 1, 1)
en_date = datetime.datetime(2010, 12, 31)

# %%
# SCA fields
# "year","month","day","sca","clearidx"
# 2000,1,1,-9.99,-999
# 2000,1,2,-9.99,-999

# %% [markdown]
# ## Check SCA for a headwater to see which HRUs meet criteria

# %%
sca_files = glob.glob(f'{workdir}/HRU_*')
sca_files.sort()

# %%
sca_files

# %%
has_sca_hrus = 0

for ff in sca_files:
    print(os.path.basename(ff))
    
    df = pd.read_csv(ff, sep=',', parse_dates={'time': [0, 1, 2]}, index_col='time', 
                     na_values=[-9.99, -999.0])
    
    df = df[st_date:en_date]
    df = df[df['sca'].notna()]
    df = df[df['clearidx'] >= 70.0]
    df = df[df['sca'] > 0.0]
    
    if len(df.index) >= 10:
        has_sca_hrus += 1
        
#     print(f'\tNumber of records: {len(df.index)}')
print(f'Number HRUs meeting SCA criteria: {has_sca_hrus}')


# %%

# %% [markdown]
# ## Set masks for each HRU in SCA baseline netcdf

# %%
# From fit_baseline_common.py
def get_dataset(filename, f_vars, start_date, end_date):
    # This routine assumes dimension nhru exists and variable nhm_id exists
    df = xr.open_dataset(filename)

    # NOTE: Next line needed if nhm_id variable exists in netcdf file
    df = df.assign_coords(nhru=df.nhm_id)

    if isinstance(f_vars, list):
        df = df[f_vars].sel(time=slice(start_date, end_date))
    else:
        df = df[[f_vars]].sel(time=slice(start_date, end_date))
    return df


# %%
baseline_file = '/Volumes/USGS_NHM1/calibrations/NHMv11/baselines/baseline_SCA_v11.nc'
sca_var = 'snow_cover_extent'
ci_var = 'sca_clear_index'

baseline_df = get_dataset(baseline_file, [sca_var, ci_var, 'nhru'], st_date, en_date)

# Remove July and August from the dataset
baseline_restr = baseline_df.sel(time=baseline_df.time.dt.month.isin([1, 2, 3, 4, 5, 6, 9, 10, 11, 12]))

# Create the SCAmask

# Compute lower and upper SCA values based on confidence interval
threshold = 70.0
ci_pct = baseline_restr[ci_var].where(baseline_restr[ci_var] >= threshold)
ci_pct /= 100.0

# Mask SCA values where CI is masked
sca_obs = baseline_restr[sca_var].where(~np.isnan(ci_pct))

# Maximum SCA value by HRU
msk_SCAmax = sca_obs.max(axis=0)

# Number of daily values > 0.0 by HRU
msk_num_obs = (sca_obs > 0.0).sum(axis=0)

# Number of years of values by HRU
msk_num_ann = sca_obs.resample(time='1AS').mean()
msk_num_ann = (msk_num_ann > 0).sum(axis=0)

# Create SCA mask based on number of years, SCAmax > 0.5, and total number of observations by HRU
SCAmask = (msk_num_ann > 1) & (msk_SCAmax > 0.5) & (msk_num_obs > 9)

# Lower bound of SCA by HRU
baseline_SCAmin = (ci_pct * sca_obs).where(SCAmask)

# Upper bound of SCA by HRU
baseline_SCAmax = (baseline_SCAmin + (1.0 - ci_pct)).where(SCAmask)


# %%
msk_num_obs

# %%
SCAmask

# %%
outfile = open(f'{workdir}/SCAmask_new', 'w')
outfile.write('sca_mask\n')

bool_vals = np.multiply(SCAmask.data, 1)

for xx in bool_vals:
    outfile.write(f'{xx}\n')
    
outfile.close()

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
col_names = ['year', 'month', 'day', 'sca', 'clearidx']
col_types = [int, int, int, float, float]
cols = dict(zip(col_names, col_types))

df = pd.read_csv(f'{workdir}/HRU_151', sep=',', parse_dates={'time': [0, 1, 2]}, index_col='time', 
                 na_values=[-9.99, -999.0])

# %%
df.info()

# %%
# Restrict to calibration date range
df = df[st_date:en_date]
print(f'Number of records: {len(df.index)}')

# %%
# poi_info = poi_info[poi_info['obs_da_actual'].notna()]
df = df[df['sca'].notna()]
print(f'Number of records: {len(df.index)}')

# %%
# Remove records where clearidx is less than 70%
df = df[df['clearidx'] >= 70.0]
print(f'Number of records: {len(df.index)}')

# %%
# df_odd.loc[:, df_odd.notna().sum(axis=0) <= 365].columns.tolist()

# %%
# Remove records where sca == 0
df = df[df['sca'] > 0.0]
print(f'Number of records: {len(df.index)}')

# %%
# Number of rows in dataframe
len(df.index)

# %%

# %%

# %%
