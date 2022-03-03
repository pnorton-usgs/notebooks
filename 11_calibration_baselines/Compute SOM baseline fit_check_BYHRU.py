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
baseline_var = 'SOM'
max_var_ann = f'soil_moist_max_norm'
min_var_ann = f'soil_moist_min_norm'
max_var_mth = f'soil_moist_max_norm'
min_var_mth = f'soil_moist_min_norm'

baseline_dir = f'/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/baselines/byHRU/{baseline_var}'

sim_dir = '/Volumes/parker_rocks/NHM_output/netcdf'
# sim_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/Maurer/byHRU_musk'

baseline_file_ann = f'{baseline_dir}/baseline_{baseline_var}ann.nc'
baseline_file_mth = f'{baseline_dir}/baseline_{baseline_var}mth.nc'

sim_var = 'soil_rechr'
sim_file = f'{sim_dir}/{sim_var}.nc'

# %%
st_date = datetime.datetime(1982,1,1)
en_date = datetime.datetime(2010,12,31)

num_yrs = (en_date.year - st_date.year) + 1

# %% [markdown]
# ### Read the baseline information

# %%
baseline_ann_df = xr.open_dataset(baseline_file_ann)
baseline_ann_df = baseline_ann_df[[max_var_ann, min_var_ann, 'nhru']].sel(time=slice(st_date, en_date))

baseline_mth_df = xr.open_dataset(baseline_file_mth)
baseline_mth_df = baseline_mth_df[[max_var_mth, min_var_mth, 'nhru']].sel(time=slice(st_date, en_date))

# %% [markdown]
# ### Read the model output from the simulation

# %%
# Model output has daily timestep
sim_df = xr.open_dataset(sim_file)

# NOTE: Next line needed if nhm_id variable exists in netcdf file
sim_df = sim_df.assign_coords(nhru=(sim_df.nhm_id))

sim_df = sim_df[[sim_var]].sel(time=slice(st_date, en_date))
# sim_df

# %% [markdown]
# #### Compute monthly total from daily model output

# %%
# Compute the monthly total from the simulation daily output
# 1MS - resample to monthly means from the start of each month
sim_month_total = sim_df.resample(time='1MS').sum()

# %% [markdown]
# #### Compute the annual mean from daily model output

# %%
# sim_ann_total = sim_df.resample(time='1AS').sum()
sim_ann_total = sim_month_total.resample(time='1AS').sum()

# %%
sim_month_total[sim_var].data.max(axis=0).shape

# %%
# Make sure the array is ordered 1..nhru for the nhm_id
hru_order = list(range(1, 109952))

ss = sim_ann_total[sim_var].loc[:, hru_order].data
ss_min = ss.min(axis=0)
ss_max = ss.max(axis=0)

# Where ss_max == ss_min increment ss_max by a small amount
# This is a bit different from how Lauren handled this
#       if (max == min) max = max + (0.01 * max)
#       if (max == min) max = max + 0.0001
ss_max = np.where(ss_max == ss_min, ss_max + (0.01 * ss_max), ss_max)
ss_max = np.where(ss_max == ss_min, ss_max + 0.0001, ss_max)

# Normalize the simulation values 0.0 to 1.0 for comparison to baseline
ss = (ss - ss_min) / (ss_max - ss_min)

bb_min_ann = baseline_ann_df[min_var_ann].loc[:, hru_order].data
bb_max_ann = baseline_ann_df[max_var_ann].loc[:, hru_order].data

ann_within_count = ((ss >= bb_min_ann) & (ss <= bb_max_ann)).sum(axis=0)
# ann_within_count = np.logical_and(ss >= bb_min_ann, ss <= bb_max_ann).sum(axis=0)
# ann_count = ss.shape[0]
ann_count = (~np.isnan(ss)).sum(axis=0)

# %%
ann_within_count

# %%
# Make sure the array is ordered 1..nhru for the nhm_id
hru_order = list(range(1, 109952))

ss = sim_month_total[sim_var].loc[:, hru_order].data

# Compute the min and max for each month for 1982-2010 by HRU
ss_min = sim_month_total[sim_var].loc[:, hru_order].groupby('time.month').min('time').data
ss_max = sim_month_total[sim_var].loc[:, hru_order].groupby('time.month').max('time').data

# Reshape from (12, nhru) to (nobs, nhru)
ss_min = np.tile(ss_min, (29, 1))
ss_max = np.tile(ss_max, (29, 1))

# Where ss_max == ss_min increment ss_max by a small amount
ss_max = np.where(ss_max == ss_min, ss_max + (0.01 * ss_max), ss_max)
ss_max = np.where(ss_max == ss_min, ss_max + 0.0001, ss_max)

# Normalize the simulation values 0.0 to 1.0 for comparison to baseline
ss = (ss - ss_min) / (ss_max - ss_min)

bb_min_mth = baseline_mth_df[min_var_mth].loc[:, hru_order].data
# bb_max_mth = np.repeat(bb_max_ann, int(bb_min_mth.shape[0] / num_yrs), axis=0)
bb_max_mth = baseline_mth_df[max_var_mth].loc[:, hru_order].data

month_within_count = ((ss >= bb_min_mth) & (ss <= bb_max_mth)).sum(axis=0)
# month_within_count = np.logical_and(ss >= bb_min, ss <= bb_max).sum(axis=0)
# month_count = ss.shape[0]
month_count = (~np.isnan(ss)).sum(axis=0)

# %%

# %%
ss.shape

# %%
ss_min.shape

# %%
month_within_count[0]

# %%
print(bb_max_ann.shape)
print(bb_max_ann[:, 0])

# %%
bb_max_mth.shape

# %%
bb_max_mth[0:14, 0]

# %%

# %%

# %%

# %% [markdown]
# ### Plot the results

# %%
import matplotlib.pyplot as plt

from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ParameterNetCDF import ParameterNetCDF

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/paramdb_daymet_CONUS'
# filename = f'{workdir}/something.nc'

hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp'
hru_layer_name = None
hru_shape_key='hru_id_nat'

seg_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/GF_nat_reg_lcc/nsegmentNationalIdentifier.shp'
seg_layer_name = None
seg_shape_key = 'seg_id_nat'

# %%
# pdb = ParameterFile(filename, verbose=True, verify=True)
# pdb = ParameterNetCDF(filename=filename, verbose=False, verify=True)
pdb = ParamDb(paramdb_dir=workdir, verbose=True, verify=False)

pdb.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)
# pdb.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %% [markdown]
# ### Add the annual total as percent within bounds

# %%
new_param = 'fit_SOM_anntot'
# pdb.parameters.remove(new_param)

pdb.parameters.add(name=new_param, datatype=2, minimum=0.0, maximum=100.0, default=0.0,
                   description='Percent of annual total values within baseline bounds')

tparam = pdb.parameters.get(new_param)

tparam.dimensions.add('nhru', size=ann_within_count.shape[0])

pct = (ann_within_count / ann_count) * 100.0
tparam.data = pct

# %% [markdown]
# ### Add the monthly total as percent within bounds

# %%
new_param = 'fit_SOM_montot'

pdb.parameters.add(name=new_param, datatype=2, minimum=0.0, maximum=100.0, default=0.0,
                   description='Percent of monthly total values within baseline bounds')

tparam = pdb.parameters.get(new_param)

tparam.dimensions.add('nhru', size=month_within_count.shape[0])

pct = (month_within_count / month_count) * 100.0
tparam.data = pct

# %%

# %%
cparam = 'fit_SOM_montot'
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

# %%

# %%

# %%
chru = 19408

print(f'month_within_count: {month_within_count[chru]}, total_count: {month_count[chru]}')
print(f'ann_within_count: {ann_within_count[chru]}, total_count: {ann_count[chru]}')

# %%
sim_month_total[sim_var].data.shape

# %%
sim_ann_total[sim_var].data.shape

# %%
