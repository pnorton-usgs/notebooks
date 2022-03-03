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
import numpy as np
import pandas as pd

from pyPRMS.prms_helpers import dparse

# %%
headwater = '0259'
hw_suffix = ''
workdir = f'/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHW_sample/HW{headwater}{hw_suffix}'
statvar_file = f'{workdir}/statvar_{headwater}'
obs_file = f'{workdir}/CAL_DATA.HW{headwater}'

col_names = ['idx', 'year', 'month', 'day', 'hour', 'minute', 'second', 'qval']
col_types = [int, int, int, int, int, int, int, float]
cols = dict(zip(col_names, col_types))

st_date = datetime.datetime(1983, 1, 1)
en_date = datetime.datetime(2009, 12, 31)

df = pd.read_csv(statvar_file, sep='\s+', skiprows=2, skipinitialspace=True, 
                 names=col_names, 
                 parse_dates={'time': [1, 2, 3]}, index_col='time', header=None)

# usecols=obs_col_names, dtype=obs_cols,
#     obs_col_names = ['site_no', 'ts_id', 'year_nu', 'mean_va']
#     obs_col_types = [np.str_, np.int16, np.int16, np.float64]
#     obs_cols = dict(zip(obs_col_names, obs_col_types))

# %%
df.drop(columns=['idx', 'hour', 'minute', 'second'], inplace=True)

# %%
df.head()

# %%
# sim_df.resample(time='1MS').mean()
df.resample('1MS').mean()

# %%
df_odd = df.copy()
df_odd['o_year'] = df_odd.index.year % 2
df_odd = df_odd[df_odd['o_year'] == 1]
df_odd = df_odd[st_date:en_date]

# %%
df_odd.head()

# %% [markdown]
# ### Monthly mean streamflow

# %%
df_mth_mean = df_odd.resample('1MS').mean()

# %%
df_mth_mean.drop(columns=['o_year'], inplace=True)
df_mth_mean.head()

# %% [markdown]
# ### Mean monthly streamflow

# %%
# monthly.resample('M').mean().groupby(monthly.index.month).mean()
df_mth_mean.resample('M').mean().groupby(df_mth_mean.index.month).mean()

# %% [markdown]
# ## Observed streamflow

# %%
col_names = ['CAL', 'year', 'month', 'day', 'efc', 'median_q', 'min_q', 'max_q']
col_types = [int, int, int, int, int, float, float, float]
cols = dict(zip(col_names, col_types))

st_date = datetime.datetime(1983, 1, 1)
en_date = datetime.datetime(2009, 12, 31)

df = pd.read_csv(obs_file, sep='\s+', skiprows=1, skipinitialspace=True, 
                 names=col_names, 
                 parse_dates={'time': [1, 2, 3]}, index_col='time', header=None)
df.drop(columns=['CAL', 'efc'], inplace=True)


# %%
df.head()

# %%
# Select only odd years
df['o_year'] = df.index.year % 2
df = df[df['o_year'] == 1]
df.drop(columns=['o_year'], inplace=True)
df = df[st_date:en_date]

# %%
df.head()

# %%
# Monthly mean streamflow
df_mth_mean = df.resample('1MS').mean()
df_mth_mean.head()

# %%
# Mean monthly streamflow
df_mth_mean.resample('M').mean().groupby(df_mth_mean.index.month).mean()

# %%

# %%
