# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
# ---

# %%
working_dir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t1/06267400/runs/TST1/-0001'
sim_file = '%s/statvar' % working_dir
obs_file = '%s/input/streamflow/06267400.data' % working_dir

# %%
import prms_lib_v3 as obs           # Library for reading the streamflow observation file
import prms_statvar_lib as statvar  # Library for reading the PRMS statvar file
import datetime
import numpy as np
import pandas as pd

# %%
reload(obs)
obs_data = obs.prms(obs_file).data

# %%
obs_data.head()

# %% [markdown]
# ### Subset the observation data to a given date range

# %%
st = datetime.datetime(1985,10,1)
en = datetime.datetime(1990,9,30)

obs_subset = obs_data[st:en]

# %% [markdown]
# #### Compute the monthly mean

# %%
obs_mon_mn = obs_subset.resample('M', how='mean')
obs_mon_mn.head()

# %% [markdown]
# #### Compute the annual mean

# %%
obs_ann = obs_subset.resample('A-SEP', how='mean')
obs_ann.head()

# %% [markdown]
# #### Load the simulation data

# %%
reload(statvar)
sim_data = statvar.statvar(sim_file).data
sim_data.drop(['rec'], axis=1, inplace=True)

# Subset to the given date range
sim_subset = sim_data[st:en]


# %% [markdown]
# #### Compute the monthly mean for the simulation data

# %%
sim_mon_mn = sim_subset.resample('M', how='mean')
sim_mon_mn.head()

# %% [markdown]
# #### Compute the annual mean for the simulation data

# %%
sim_ann = sim_subset.resample('A-SEP', how='mean')
sim_ann.describe()

# %%
a = pd.Series(sim_ann.iloc[:,0])
b = pd.Series(obs_ann.iloc[:,0])

# %%
b - a

# %%

# %% [markdown]
# #### Daily NRMSE

# %%
daily_diff = obs_subset.iloc[:,0] - sim_subset.iloc[:,0]

# Mean daily 
mn_daily = obs_subset.mean(axis=0)

dsq = daily_diff**2
ndsq = (obs_subset.iloc[:,0] - mn_daily.iloc[0])**2

nrmse_daily = np.sqrt(dsq.sum() / ndsq.sum())
print nrmse_daily

# %% [markdown]
# #### Monthly NRMSE

# %%
monthly_mn_diff = obs_mon_mn.iloc[:,0] - sim_mon_mn.iloc[:,0]
monthly_mn_diff.head()

lt_month_mn = obs_mon_mn.mean(axis=0)

nrmse_monthly = np.sqrt((monthly_mn_diff**2).sum() / ((obs_mon_mn.iloc[:,0] - lt_month_mn.iloc[0])**2).sum())
print nrmse_monthly

# %% [markdown]
# #### Annual NRMSE

# %%
annual_diff = obs_ann.iloc[:,0] - sim_ann.iloc[:,0]
annual_diff.head()

mn_annual = obs_ann.mean(axis=0)

nrmse_annual = np.sqrt((annual_diff**2).sum() / ((obs_ann.iloc[:,0] - mn_annual.iloc[0])**2).sum())
print nrmse_annual

# %% [markdown]
# #### Compute mean monthly

# %%
# Compute the mean monthly values for observations
obs_mn_monthly = obs_mon_mn.copy()
obs_mn_mon = obs_mn_monthly.groupby(obs_mn_monthly.index.month).mean().iloc[:,0]

# %%
# Compute the mean monthly values for simulated values
sim_mn_monthly = sim_mon_mn.copy()
sim_mn_mon = sim_mn_monthly.groupby(sim_mn_monthly.index.month).mean().iloc[:,0]


# %% [markdown]
# #### Compute NRMSE for mean monthly

# %%
mn_monthly_diff = obs_mn_mon.iloc[:] - sim_mn_mon.iloc[:]
#print mn_monthly_diff.head()

lt_mn_month = obs_mn_mon.mean(axis=0)
#print lt_mn_month

nrmse_mn_monthly = np.sqrt((mn_monthly_diff**2).sum() / ((obs_mn_mon.iloc[:] - lt_mn_month)**2).sum())
print nrmse_mn_monthly

# %%
