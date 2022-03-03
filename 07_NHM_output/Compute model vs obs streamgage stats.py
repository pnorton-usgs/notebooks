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
# import dask
import datetime
# import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# from collections import OrderedDict

# from pyPRMS.prms_helpers import dparse
# from pyPRMS.ParameterFile import ParameterFile

# %%
# NHMv1.0
# calib_name = 'byHRU_noroute_obs'  # one of: PRECAL, byHRU, ...
# workdir = f'/Volumes/USGS_NHM1/calibrations/NHMv10/DAYMET_releases/{calib_name}'
# param_filename = f'{workdir}/NHM-PRMS.param'
# filename = f'{workdir}/output/NHM-PRMS_data_release.nc'

# NHMv1.1
calib_name = 'byHWobs'  # one of: PRECAL, byHRU, ...
# workdir = '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHRU/20210624_calib'
# workdir = '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHW/20211027_gm_byHW'
workdir = '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHWobs/20211129_gm_byHWobs'
param_filename = f'{workdir}/myparam.param'
filename = f'{workdir}/output_variables/stats.nc'

falcone_dir = '/Users/pnorton/GIS/gagesII_additiona_data/basinchar_and_report_sept_2011'

st_date = datetime.datetime(1979, 10, 1)
en_date = datetime.datetime(2019, 9, 30)

# %%
### Read the model/obs netcdf file
# Open the streamflow observation files
poi_xdf = xr.open_mfdataset(filename, chunks={'poi_id': 2000}, decode_cf=True, engine='netcdf4')
type(poi_xdf)

# %%

# %% [markdown]
# ## Create a POI information dataframe
# Add drainage area related information

# %%
poi_df = poi_xdf[['poi_name', 'latitude', 'longitude', 'drainage_area', 'drainage_area_contrib', 'drainage_area_model']].to_pandas()
# poi_df.rename(columns={'drainage_area': 'da_obs', 'drainage_area_contrib': 'da_contrib_obs'}, inplace=True)

# Compute the correction factor for obs values
poi_df['da_ratio_obs'] = poi_df['drainage_area_contrib'] / poi_df['drainage_area']

# Any NaN values should default to a correction factor of 1.0
poi_df['da_ratio_obs'].fillna(value=1.0, inplace=True)

# Sometimes the full da and contributing da are swapped
poi_df['da_ratio_obs'] = np.where(poi_df['da_ratio_obs'] > 1.0, 1.0 / poi_df['da_ratio_obs'], poi_df['da_ratio_obs'])

poi_df['da_actual_obs'] = poi_df[['drainage_area', 'drainage_area_contrib']].min(axis=1)

# Add drainage area ratio between nwis/hydat and NHM
poi_df['da_ratio'] = poi_df['da_actual_obs'] / poi_df['drainage_area_model']
poi_df['da_ratio'] = np.where(poi_df['da_ratio'] > 1.0, 1.0 / poi_df['da_ratio'], poi_df['da_ratio'])

# Add a drainage area correction factor
poi_df['da_correction'] = poi_df['da_actual_obs'] / poi_df['drainage_area_model']

poi_df.head(10)

# %%

# %%

# %% [markdown]
# ## Read the Falcone information

# %%
col_names = ['STAID', 'CLASS', 'HYDRO_DISTURB_INDX']
col_types = [str, str, int]
cols = dict(zip(col_names, col_types))

falcone_df = pd.read_excel(open(f'{falcone_dir}/gagesII_sept30_2011_conterm.xlsx', 'rb'), sheet_name='Bas_Classif', 
                           usecols=[0, 1], dtype=cols)

falcone_df.rename(columns={'STAID': 'poi_id', 'CLASS': 'falcone_class'}, inplace=True)
falcone_df.set_index('poi_id', inplace=True)
falcone_df.info()

falcone_ids = falcone_df.index.tolist()
falcone_ref = falcone_df[falcone_df['falcone_class'] == 'Ref']
falcone_ref_ids = falcone_ref.index.tolist()

# %%
falcone_df.head()

# %%
poi_df = pd.merge(poi_df, falcone_df, how='left', left_index=True, right_index=True)
poi_df.head(20)

# %%
poi_df.loc['02217500']

# %%

# %% [markdown]
# ## Remove POIs not meeting criteria

# %%
# Remove POIs that lack a DA 
df_reduce_1 = poi_df[poi_df['da_actual_obs'].notna()]

# Remove POIs where the drainage area ratio is less than 0.9
df_reduce_2 = df_reduce_1[df_reduce_1['da_ratio'] >= 0.9]

# List of POIs meeting all criteria
poi_list = df_reduce_2.index.tolist()

# %%
df_reduce_1.info()

# %%
df_reduce_2.info()


# %% [markdown]
# ## Define objective functions

# %%
def nashsutcliffe(obs, sim):
    numerator = np.sum(np.power((sim - obs), 2.0), axis=1)
    mean_obs = obs.mean(skipna=True, axis=1)
    denominator = np.sum(np.power(obs.subtract(mean_obs, axis=0), 2.0), axis=1)
    ns = 1.0 - numerator / denominator
    
    return ns

def bias(obs, sim):
    # Bias
    x1 = sim + 1.0
    x2 = obs + 1.0
    bias = np.sum(np.abs((x1 - x2) / x2), axis=1) / obs.count(axis=1)

    return bias


# %% [markdown]
# ## Get observed streamflow

# %%
obs_df = poi_xdf['discharge'].loc[poi_list, :].to_pandas()
obs_df.head()

# %% [markdown]
# ## Process the model streamflow

# %%
raw_sim_df = poi_xdf['discharge_model'].loc[poi_list, :].to_pandas()

# %% [markdown]
# ### Adjust model streamflow based on drainage area correction between obs and model

# %%
adj_sim_df = raw_sim_df.mul(df_reduce_2.loc[:, 'da_correction'], axis=0)
adj_sim_df.head()

# %% [markdown]
# ## Compute statistics for daily streamflow

# %%
# Nashsutcliffe of the adjusted daily sim values compared to obs
ns = nashsutcliffe(obs_df, adj_sim_df)
ns.name = 'ns'
ns.head()

# %%
# Nashsutcliffe of the log of adjusted daily sim values compared to obs
zero_adj = 0.0001
log_obs_df = np.log(obs_df + zero_adj)

log_adj_sim_df = np.log(adj_sim_df + zero_adj)

nslog = nashsutcliffe(log_obs_df, log_adj_sim_df)
nslog.name = 'nslog'
nslog.head(20)

# %%
# Bias of the daily values
daily_bias = bias(obs_df, adj_sim_df)
daily_bias.name = 'bias'
daily_bias.head(20)

# %% [markdown]
# ## Compute monthly mean for observations

# %%
# 1MS - resample to monthly means from the start of each month
mon_obs_df = obs_df.resample('1MS', axis=1).mean()
mon_obs_df.head()

# %% [markdown]
# ## Compute monthly mean for simulated

# %%
# 1MS - resample to monthly means from the start of each month
mon_adj_sim_df = adj_sim_df.resample('1MS', axis=1).mean()
mon_adj_sim_df.head()

# %%
mon_ns = nashsutcliffe(mon_obs_df, mon_adj_sim_df)
mon_ns.name = 'mon_ns'
mon_ns.head(20)

# %% [markdown]
# ## Compute NS of log of simulated vs. obs

# %%
zero_adj = 0.0001
log_mon_obs_df = np.log(mon_obs_df + zero_adj)
# log_mon_obs_df.head()

log_mon_adj_sim_df = np.log(mon_adj_sim_df + zero_adj)
# log_mon_adj_sim_df.head()

# %%
mon_nslog = nashsutcliffe(log_mon_obs_df, log_mon_adj_sim_df)
mon_nslog.name = 'mon_nslog'
mon_nslog.head(20)

# %%
# Bias of the monthly values
mon_bias = bias(mon_obs_df, mon_adj_sim_df)
mon_bias.name = 'mon_bias'
mon_bias.head(20)

# %% [markdown]
# ## Write the stats to a file

# %%
fields = ['latitude', 'longitude', 'da_actual_obs', 'drainage_area_model', 'falcone_class']
crap = df_reduce_2.loc[:, fields]

crap = pd.merge(crap, ns, left_index=True, right_index=True)
crap = pd.merge(crap, nslog, left_index=True, right_index=True)
crap = pd.merge(crap, mon_ns, left_index=True, right_index=True)
crap = pd.merge(crap, mon_nslog, left_index=True, right_index=True)
crap = pd.merge(crap, daily_bias, left_index=True, right_index=True)
crap = pd.merge(crap, mon_bias, left_index=True, right_index=True)

# %%
col_order = ['ns', 'nslog', 'mon_ns', 'mon_nslog', 'bias', 'mon_bias', 
             'da_actual_obs', 'drainage_area_model', 'falcone_class', 'latitude', 'longitude']
crap.to_csv(f'{workdir}/gage_stats_{calib_name}.csv', sep=',', index=True, columns=col_order)


# %%

# %% [markdown]
# ## Create exceedance plot data

# %%
def get_exceedance_curve(ts, fclass):
    stat_name = f'{ts.name}_{fclass}'
    
    ranked_stat = sorted(ts[ts.notnull()])
    prob = np.arange(len(ranked_stat), dtype=float) + 1.0
    prob = (prob / (len(ranked_stat) + 1.0))
    
    # Return dataframe of exceedence curve values
    return pd.DataFrame({'exceedance': prob, stat_name: ranked_stat}, columns=['exceedance', stat_name])


# %%

# %%

# %%

# %%

# %%

# %% [markdown]
# ## Plot stuff

# %%
cmap = plt.cm.get_cmap('RdYlBu', 5)
norm = Normalize(vmin=crap['ns'].min(), vmax=crap['ns'].max())

# %%
# df.plot(x="longitude", y="latitude", kind="scatter", c="brightness", colormap="YlOrRd")
crap.plot(x='longitude', y='latitude', kind='scatter', c='ns', figsize=(15,10), colormap='YlOrRd')

# %%
crap.min()

# %%

# %%
