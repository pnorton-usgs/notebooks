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
nhm_version = 1.0
baseline_info = {'AET': {'max': 'aet_max', 'min': 'aet_min', 'prms_var': 'hru_actet'},
                 'RCH': {'max': 'recharge_max_norm', 'min': 'recharge_min_norm', 'prms_var': 'recharge'},
                 'RUN': {'max': 'runoff_max', 'min': 'runoff_min', 'prms_var': 'hru_outflow'},
                 'SCA': {'sca': 'snow_cover_extent', 'scaci': 'sca_clear_index', 'prms_var': 'snowcov_area'},
                 'SOMann': {'max': 'soil_moist_max_norm', 'min': 'soil_moist_min_norm', 'prms_var': 'soil_rechr'},
                 'SOMmth': {'max': 'soil_moist_max_norm', 'min': 'soil_moist_min_norm', 'prms_var': 'soil_rechr'}}

baseline_var = 'SCA'
sca_var = baseline_info[baseline_var]['sca']
ci_var = baseline_info[baseline_var]['scaci']

if nhm_version == 1.0:
    baseline_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/baselines'

    sim_dir = '/Volumes/USGS_NHM1/NHM_output/NHMv10/byHRU/netcdf'
    # sim_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/Maurer/byHRU_musk'

    baseline_file = f'{baseline_dir}/baseline_{baseline_var}_v10.nc'
elif nhm_version == 1.1:
    baseline_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/baselines'

    sim_dir = '/Volumes/USGS_NHM1/NHM_output/NHMv11/byHRU'
    # sim_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/Maurer/byHRU_musk'

    baseline_file = f'{baseline_dir}/baseline_{baseline_var}_v11.nc'

sim_var = baseline_info[baseline_var]['prms_var']
sim_file = f'{sim_dir}/{sim_var}.nc'

# %%
st_date = datetime.datetime(2000,1,1)
en_date = datetime.datetime(2010,12,31)

# Mask July and August baseline values or just zero them out
remove_ja = False

# %% [markdown]
# ### Read the baseline information

# %%
# Baseline recharge has daily timestep
baseline_df = xr.open_dataset(baseline_file)
baseline_df = baseline_df.assign_coords(nhru=(baseline_df.nhm_id))
baseline_df = baseline_df[[sca_var, ci_var, 'nhru']].sel(time=slice(st_date, en_date))

# Remove July and August from the dataset
if remove_ja:
    baseline_restr = baseline_df.sel(time=baseline_df.time.dt.month.isin([1, 2, 3, 4, 5, 6, 9, 10, 11, 12]))
else:
    baseline_restr = baseline_df

# %%
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


# %%
# Lower bound of SCA by HRU
baseline_SCAmin = (ci_pct * sca_obs).where(SCAmask)

# Upper bound of SCA by HRU
baseline_SCAmax = (baseline_SCAmin + (1.0 - ci_pct)).where(SCAmask)

# %%
print(f'SCAci: min={np.nanmin(ci_pct)}, max={np.nanmax(ci_pct)}')
print(f'SCAmin: min={np.nanmin(baseline_SCAmin)}, max={np.nanmax(baseline_SCAmin)}')
print(f'SCAmax: min={np.nanmin(baseline_SCAmax)}, max={np.nanmax(baseline_SCAmax)}')

# %% [markdown]
# ### Read the model output from the simulation

# %%
# Model output has daily timestep
sim_df = xr.open_dataset(sim_file)

# NOTE: Next line needed if nhm_id variable exists in netcdf file
sim_df = sim_df.assign_coords(nhru=(sim_df['nhm_id']))

sim_df = sim_df.sel(time=slice(st_date, en_date))

# List ordered 1..nhru
hru_order = list(range(1, sim_df['nhru'].size + 1))

if remove_ja:
    sim_restr = sim_df.sel(time=sim_df.time.dt.month.isin([1, 2, 3, 4, 5, 6, 9, 10, 11, 12]))
else:
    sim_restr = sim_df

# %% [markdown]
# ### Compute the by-HRU counts of simulation values that are within baseline bounds

# %%
ss = sim_restr[sim_var].loc[:, hru_order].data

# Baseline upper and lower values ordered by nhm_id
bb_min = baseline_SCAmin.loc[:, hru_order].data
bb_max = baseline_SCAmax.loc[:, hru_order].data

daily_within_count = ((ss >= bb_min) & (ss <= bb_max)).sum(axis=0)

# daily_count = ss.shape[0]
daily_count = (~np.isnan(bb_min)).sum(axis=0)
ss_count = (~np.isnan(ss)).sum(axis=0)

print(f'daily_count: min={np.nanmin(daily_count)}, max={np.nanmax(daily_count)}')
print(f'daily_within_count: min={np.nanmin(daily_within_count)}, max={np.nanmax(daily_within_count)}')

print(f'ss_count: min={np.nanmin(ss_count)}, max={np.nanmax(ss_count)}')

# %% [markdown]
# # Plot the results

# %%
import matplotlib.pyplot as plt

from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ParameterNetCDF import ParameterNetCDF

# %%
if nhm_version == 1.0:
    workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/paramdb_daymet_CONUS'
    # filename = f'{workdir}/something.nc'

    hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp'
    hru_layer_name = None
    hru_shape_key='hru_id_nat'

    seg_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/GF_nat_reg_lcc/nsegmentNationalIdentifier.shp'
    seg_layer_name = None
    seg_shape_key = 'seg_id_nat'
elif nhm_version == 1.1:
    workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS'
    
    hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
    hru_layer_name = 'nhruv11_sim30'
    hru_shape_key='nhru_v11'

    # Segment lines
    seg_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
    seg_layer_name = 'nsegment_v11'
    seg_shape_key = 'nsegment_v11'

# %%
# pdb = ParameterFile(filename, verbose=True, verify=True)
# pdb = ParameterNetCDF(filename=filename, verbose=False, verify=True)
pdb = ParamDb(paramdb_dir=workdir, verbose=True, verify=False)

pdb.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)
# pdb.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %% [markdown]
# ### Add the monthly runoff climatology as percent within bounds

# %%
new_param = 'fit_snowcov_area_daily'
# pdb.parameters.remove(new_param)

pdb.parameters.add(name=new_param, datatype=2, minimum=0.0, maximum=100.0, default=0.0,
                   description='Percent of daily SCA values within baseline bounds')

tparam = pdb.parameters.get(new_param)

tparam.dimensions.add('nhru', size=daily_within_count.shape[0])

pct = (daily_within_count / daily_count) * 100.0
tparam.data = pct

# %% [markdown]
# ### Plot the baseline results

# %%
cparam = 'fit_snowcov_area_daily'
use_orig_colors = True

if use_orig_colors:
    # Colormap using Lauren's original colors for the T&M
    from matplotlib.colors import LinearSegmentedColormap

    colors = [(0, 0, 1), (0.13, 1, 1), (0.13, 1, 0.02), (0.99, 0.58, 0.04), (0.98, 0, 1)]
    n_bin = 5
    cmap_name = 'lauren_col'
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
else:
    cmap = plt.cm.get_cmap('RdYlBu', 5)

pdb.parameters.plot(cparam, limits='valid', linewidth=0.0, edgecolor='whitesmoke', mask_defaults='darkgrey', cmap=cmap)

# %%
chru = 0

print(f'daily_within_count: {daily_with_count[chru]}, total_count: {daily_count[chru]}')

# %%

# %%
