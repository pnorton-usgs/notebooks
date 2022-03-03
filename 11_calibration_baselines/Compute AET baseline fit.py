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

baseline_var = 'AET'
calib = 'byHRU_noroute_obs'
max_var = baseline_info[baseline_var]['max']
min_var = baseline_info[baseline_var]['min']

if nhm_version == 1.0:
    baseline_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/baselines'

    sim_dir = f'/Volumes/USGS_NHM1/calibrations/NHMv10/DAYMET_releases/{calib}'
    # sim_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/Maurer/byHRU_musk'

    baseline_file = f'{baseline_dir}/baseline_{baseline_var}_v10.nc'
elif nhm_version == 1.1:
    baseline_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/baselines'

    sim_dir = '/Volumes/USGS_NHM1/NHM_output/NHMv11/byHRU'
    # sim_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/Maurer/byHRU_musk'

    baseline_file = f'{baseline_dir}/baseline_{baseline_var}_v11.nc'

sim_var = baseline_info[baseline_var]['prms_var']
sim_file = f'{sim_dir}/{sim_var}.nc'

out_csv = f'/Volumes/USGS_NHM1/calibrations/NHMv10/DAYMET_releases/baseline_fit_files/{baseline_var}_fitBASELINE_{calib}'

# %%
st_date = datetime.datetime(2000,1,1)
en_date = datetime.datetime(2009,12,31)

# %% [markdown]
# ### Read the baseline information

# %%
baseline_df = xr.open_dataset(baseline_file)
baseline_df = baseline_df.assign_coords(nhru=(baseline_df.nhm_id))
baseline_df = baseline_df[[max_var, min_var, 'nhru']].sel(time=slice(st_date, en_date))

# %% [markdown]
# #### Compute the monthly climatology for the baseline monthly mean

# %%
# Compute the monthly climatology for the baseline
baseline_mean_month = baseline_df.groupby('time.month').mean('time')
# baseline_mean_month
# mean_month_base = base_df.groupby('time.month').mean('time')

# %% [markdown]
# ### Read the model output from the simulation

# %%
# Read the simulation results
sim_df = xr.open_dataset(sim_file)

# NOTE: Next line needed if nhm_id variable exists in netcdf file
sim_df = sim_df.assign_coords(nhru=(sim_df.nhm_id))

sim_df = sim_df[[sim_var]].sel(time=slice(st_date, en_date))

# List ordered 1..nhru
hru_order = list(range(1, sim_df['nhru'].size + 1))

# %% [markdown]
# #### Compute monthly mean from daily model output

# %%
# Compute the monthly mean for the simulation
# 1MS - resample to monthly means from the start of each month
sim_month_mean = sim_df.resample(time='1MS').mean()
# sim_month_mean

# %% [markdown]
# #### Compute the monthly climatology from the simulation monthly mean

# %%
# Compute the monthly climatology for the simulation
sim_mean_month = sim_month_mean.groupby('time.month').mean('time')
# sim_mean_month

# %% [markdown]
# ### Compute the by-HRU counts of simulation monthly climatology values that are within baseline bounds

# %%
ss = sim_mean_month[sim_var].loc[:, hru_order].data
bb_min = baseline_mean_month[min_var].loc[:, hru_order].data
bb_max = baseline_mean_month[max_var].loc[:, hru_order].data

clim_within_count = np.logical_and(ss >= bb_min, ss <= bb_max).sum(axis=0)
clim_month_count = ss.shape[0]

# %% [markdown]
# ### Compute the by-HRU counts of simulation monthly mean values that are within baseline bounds

# %%
ss = sim_month_mean[sim_var].loc[:, hru_order].data
bb_min = baseline_df[min_var].loc[:, hru_order].data
bb_max = baseline_df[max_var].loc[:, hru_order].data

monthly_within_count = np.logical_and(ss >= bb_min, ss <= bb_max).sum(axis=0)
monthly_month_count = ss.shape[0]

# %% [markdown]
# # Plot the results
#

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
new_param = 'fit_aet_monclim'
# pdb.parameters.remove(new_param)

pdb.parameters.add(name=new_param, datatype=2, minimum=0.0, maximum=100.0, default=0.0,
                   description='Percent of monthly climatology AET values within baseline bounds')

tparam = pdb.parameters.get(new_param)

tparam.dimensions.add('nhru', size=clim_within_count.shape[0])

pct = (clim_within_count / clim_month_count) * 100.0
tparam.data = pct

# %% [markdown]
# ### Add the monthly runoff climatology as percent within bounds

# %%
new_param = 'fit_aet_monthly'
# pdb.parameters.remove(new_param)

pdb.parameters.add(name=new_param, datatype=2, minimum=0.0, maximum=100.0, default=0.0,
                   description='Percent of monthly mean AET values within baseline bounds')

tparam = pdb.parameters.get(new_param)

tparam.dimensions.add('nhru', size=monthly_within_count.shape[0])

pct = (monthly_within_count / monthly_month_count) * 100.0
tparam.data = pct

# %%

# %%
cparam = 'fit_aet_monthly'
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

pdb.parameters.plot(cparam, limits='valid', linewidth=0.0, edgecolor='whitesmoke', cmap=cmap)    
# pdb.parameters.plot(cparam, limits='valid', linewidth=0.0, edgecolor='whitesmoke', mask_defaults='darkgrey', cmap=cmap)

# %%
# Write the legacy fit output files (but as csv)
fhdl = open(out_csv, 'w')
fhdl.write('hru,lat,lon,monthly_within,monthly_total,mn_monthly_within,mn_monthly_total\n')
lats = pdb.parameters['hru_lat'].tolist()
lons = pdb.parameters['hru_lon'].tolist()

# monthly_within_count
# monthly_month_count
# clim_within_count
# clim_month_count

for chru in hru_order:
    idx = chru - 1
    fhdl.write(f'{chru},{lats[idx]},{lons[idx]},{monthly_within_count[idx]},{monthly_month_count},{clim_within_count[idx]},{clim_month_count}\n')
    
fhdl.close()

# %%
out_csv

# %%
