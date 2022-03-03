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
import os
import pandas as pd
import xarray as xr

from collections import OrderedDict

from pyPRMS.ControlFile import ControlFile
from pyPRMS.ParameterFile import ParameterFile


# %%
# Taken from nhm_utilities/nhm_streamflow_to_netcdf.py
def read_csv_header(filename):
    fhdl = open(filename, 'r')

    # First and second rows are headers
    hdr1 = fhdl.readline().strip()

    fhdl.close()

    tmp_flds = hdr1.split(' ')
    tmp_flds.remove('Date')

    flds = {nn+3: hh for nn, hh in enumerate(tmp_flds)}

    # poi_flds maps column index to POI and is used to rename the dataframe columns from indices to station IDs
    poi_flds = OrderedDict()

    # poi_seg_flds maps POI to the related segment ID
    poi_seg_flds = OrderedDict()

    for xx, yy in flds.items():
        tfld = yy.split('_')
        segid = int(tfld[2]) - 1  # Change to zero-based indices
        poiid = tfld[4]

        poi_flds[xx] = poiid
        poi_seg_flds[poiid] = segid

    return poi_flds, poi_seg_flds


def read_streamflow(filename, field_names):
    df = pd.read_csv(filename, sep='\s+', header=None, skiprows=2, parse_dates={'time': [0, 1, 2]},
                     index_col='time')

    df.rename(columns=field_names, inplace=True)

    return df


# %%

# %%

# %%

# %%

# %%

# %%

# %%
hw = '0833'
workdir = f'/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHW_obs_sample/hw_{hw}'

poi_data = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/poi_data/*_pois.nc'
ctl_filename = f'{workdir}/control.default.bandit'
# ctl_filename = f'{workdir}/control_HW{hw}'

# %%
## Read observed streamflow
poi_xdf = xr.open_mfdataset(poi_data, chunks={'poi_id': 2000}, combine='nested', concat_dim='poi_id', 
                            decode_cf=True, engine='netcdf4')
poi_xdf

# %%

# %%

# %%
ctl = ControlFile(ctl_filename)

sim_filename = f'{workdir}/{ctl.get("csv_output_file").values}'
print(f'Simulated streamflow filename: {sim_filename}')

# %%
# Read the simulated streamflow
poi_flds, poi_seg_flds = read_csv_header(sim_filename)
model_pois_segments = list(poi_seg_flds.values())

sim_df = read_streamflow(sim_filename, field_names=poi_flds).T
sim_df.index.name = 'poi_id'
sim_df.head()

# %%
model_pois = sim_df.index.tolist()

st_date = min(sim_df.columns.tolist())
en_date = max(sim_df.columns.tolist())

print(f'Number of POIs: {len(model_pois)}')
print(f'Start date: {st_date}')
print(f'End date: {en_date}')

# %%

# %%
# Read the model parameter file
param_file = os.path.normpath(os.path.join(f'{workdir}/{ctl.get("param_file").values}'))
pfile = ParameterFile(param_file, verbose=True, verify=True)

# Get the list of POI IDs in the parameter file
poi_gage_id = pfile.parameters['poi_gage_id']
poi_gage_id_list = poi_gage_id.data.tolist()

# Get the cumulative drainage area for each of the POI segments
seg_cum_area = pfile.parameters.get_dataframe('seg_cum_area')

# Get model area for each POI in square miles
# using iloc because model_pois_segments are 0-based local indices
model_pois_area = seg_cum_area.iloc[model_pois_segments].to_numpy(dtype=float) / 640.0

# %%

# %%

# %%

# %%

# %%
# Get the POI information
poi_info = poi_xdf[['poi_name', 'latitude', 'longitude', 'drainage_area', 'drainage_area_contrib']].to_pandas()

# Restrict to the POIs in the model
poi_info = poi_info.loc[model_pois]

# Add the model drainage area
poi_info['drainage_area_model'] = model_pois_area

# Compute the correction factor for obs values; NaN obs_da_ratio defaults to 1.0
poi_info['obs_da_ratio'] = poi_info['drainage_area_contrib'] / poi_info['drainage_area']
poi_info['obs_da_ratio'].fillna(value=1.0, inplace=True)

# Sometimes the full da and contributing da are swapped in NWIS
poi_info['obs_da_ratio'] = np.where(poi_info['obs_da_ratio'] > 1.0, 1.0 / poi_info['obs_da_ratio'], poi_info['obs_da_ratio'])

poi_info['obs_da_actual'] = poi_info[['drainage_area', 'drainage_area_contrib']].min(axis=1)

# Add drainage area ratio between nwis/hydat and NHM
poi_info['sim_da_ratio'] = poi_info['obs_da_actual'] / poi_info['drainage_area_model']
poi_info['sim_da_ratio'] = np.where(poi_info['sim_da_ratio'] > 1.0, 1.0 / poi_info['sim_da_ratio'], poi_info['sim_da_ratio'])

# Add a drainage area correction factor
poi_info['sim_da_correction'] = poi_info['obs_da_actual'] / poi_info['drainage_area_model']

poi_info.info()

# %%
# Remove POIs which are missing more than a given number of obs
# df_reduce_1 = poi_info[poi_info['missing_obs'] < 3650]

# Remove POIs that lack a DA 
df_reduce_2 = poi_info[poi_info['obs_da_actual'].notna()]

# Remove POIs where obs-to-sim da ratio < 0.9
df_reduce_3 = df_reduce_2[df_reduce_2['sim_da_ratio'] >= 0.9]

# %%
df_reduce_2

# %%
df_reduce_3

# %%

# %%

# %%

# %%

# %%

# %%
obs_df = poi_xdf['discharge'].loc[model_pois, st_date:en_date].to_pandas()

# Drop POIs where there are no valid values
obs_df.dropna(how='all', axis=0, inplace=True)

# Get the number of valid values for each streamgage
obs_cnts = obs_df.notna().sum(axis=1)
obs_cnts = obs_cnts[obs_cnts > 365]
obs_cnts.name = 'obs_cnts'
print(obs_cnts)

included_pois = obs_cnts.index.tolist()
print(f'\nPOIs meeting calibration criteria: {included_pois}\n')

print(obs_df.head())

# %%

# %% [markdown]
# ## Reduce poi_info, sim_df and obs_df dataframes to just the included_pois

# %%
poi_info = poi_info.loc[included_pois]

obs_df = obs_df.loc[included_pois]

sim_df = sim_df.loc[included_pois]

# %%
# Get the number of valid values for each streamgage
sim_cnts = sim_df.notna().sum(axis=1)
sim_cnts.name = 'sim_cnts'
print(sim_cnts)

## Adjust the simulated values by the ratio of obs_da vs. sim_da for each POI
adj_sim_df = sim_df.mul(poi_info.loc[:, 'sim_da_correction'], axis=0)
adj_sim_df.head()

# %%
# Ratio of valid obs values vs valid sim values
val_ratio = obs_cnts / sim_cnts
val_ratio.name = 'val_ratio'

print(val_ratio)

# %%
# Compute the drainage area weights for each POI relative to the total drainage area of the POIs
obs_total_da = poi_info['obs_da_actual'].sum()
obs_total_da

poi_info['obs_weight'] = poi_info['obs_da_actual'] / obs_total_da
poi_info


# %%

# %%
def nashsutcliffe(obs, sim, weight=None):
    numerator = np.sum(np.power((sim - obs), 2.0), axis=1)
    mean_obs = obs.mean(skipna=True, axis=1)
    denominator = np.sum(np.power(obs.subtract(mean_obs, axis=0), 2.0), axis=1)
    ns = 1.0 - numerator / denominator
    
    if weight is not None:
        # Adjust where NS value where NS < 0.0 when weight is specified
        # This adjustment was in the original NS code for the byHW_obs calibration
        ns.where(ns >= 0.0, ns * weight, inplace=True)
    return ns


# %% [markdown]
# ## Compute Nash-Sutcliffe on daily and log(daily) sim/obs values

# %%
ns = nashsutcliffe(obs_df, adj_sim_df, weight=0.1)
ns.name = 'ns'
ns.head()

# %%
# Nashsutcliffe of the log of adjusted daily sim values compared to obs
zero_adj = 0.0001
log_obs_df = np.log(obs_df + zero_adj)

log_adj_sim_df = np.log(adj_sim_df + zero_adj)

nslog = nashsutcliffe(log_obs_df, log_adj_sim_df, weight=0.1)
nslog.name = 'nslog'
nslog.head()


# %% [markdown]
# ## Compute the objective function

# %%
# compute_OF(use_gage, wght, wghtNUM, NS, NSlog, DAcor, NSsum, prmsOF)
#       NSsum = 0.0
#       num_gages = size(use_gage)

#       do ii = 1, num_gages
#         if (use_gage(ii) == 0) then
#           w = 3.0
#           wlog = 1.0
#           x = (w * (1.0 - NS(ii))) + (wlog * (1.0 - NSlog(ii)))

# #           ! wght is the ratio of each gage area to the total gage area
# #           ! wghtNUM is ratio of number of obs values to number of sim values
#           NSsum = NSsum + (DAcor(ii) * wghtNUM(ii) * wght(ii) * x)
#         end if
#       end do

#       prmsOF = NSsum / (float(num_gages))

# %%
def compute_objfcn(obs_weights, numvals_ratio, ns, nslog, da_ratio):
    # This routine does not implement the use_gage filtering criteria 
    # as written in byHW_obs code
    w = 3.0
    wlog = 1.0
    x = (w * (1.0 - ns)) + (wlog * (1.0 - nslog))
    
    ns_sum = (da_ratio * numvals_ratio * obs_weights * x).sum()
    objfcn = ns_sum / len(ns)
    print(f'ns_sum = {ns_sum}')
    print(f'objfcn = {objfcn}')


# %%
compute_objfcn(poi_info['obs_weight'], val_ratio, ns, nslog, poi_info['sim_da_ratio'])

# %%

# %%

# %%

# %%

# %%

# %%
