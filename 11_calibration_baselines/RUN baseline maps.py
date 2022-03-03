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
baseline_var = 'RUN'

max_var = 'runoff_max'
min_var = 'runoff_min'

# baseline_dir = f'/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/baselines/conus/{baseline_var}'
baseline_dir = f'/Volumes/parker_rocks/calibrations/NHMv11/baselines/{baseline_var}'

baseline_file = f'{baseline_dir}/baseline_{baseline_var}.nc'


# %%
st_date = datetime.datetime(1982,1,1)
en_date = datetime.datetime(2010,12,31)

# %% [markdown]
# ### Read the baseline information

# %%
baseline_df = xr.open_dataset(baseline_file)
baseline_df = baseline_df.assign_coords(nhru=(baseline_df.nhm_id))
baseline_df = baseline_df[[max_var, min_var, 'nhru']].sel(time=slice(st_date, en_date))

# %%
baseline_df

# %% [markdown]
# #### Compute the monthly climatology for the baseline monthly mean

# %%
baseline_month_clim = baseline_df.groupby('time.month').mean('time')
# baseline_month_clim

# %% [markdown]
# #### Compute monthly mean from daily model output

# %%
# Compute the monthly mean from the simulation daily output
# 1MS - resample to monthly means from the start of each month
sim_month_mean = sim_df.resample(time='1MS').mean()
# sim_month_mean

# %% [markdown]
# #### Compute the monthly climatology from the simulation monthly mean

# %%
# Compute the monthly climatology for the simulation
sim_month_clim = sim_month_mean.groupby('time.month').mean('time')
# sim_month_clim

# Number of observations of each month for period of record
# sim_month_clim.groupby('time.month').count('time')['hru_outflow']

# Total number of observations for monthly climatology
# sim_month_clim.count('time')

# %% [markdown]
# ### Compute the by-HRU counts of simulation monthly climatology values that are within baseline bounds

# %%
# Make sure the array is ordered 1..nhru for the nhm_id
hru_order = list(range(1, 109952))

ss = sim_month_clim[sim_var].loc[:, hru_order].data
bb_min = baseline_month_clim[min_var].loc[:, hru_order].data
bb_max = baseline_month_clim[max_var].loc[:, hru_order].data

clim_within_count = np.logical_and(ss >= bb_min, ss <= bb_max).sum(axis=0)
clim_month_count = ss.shape[0]

# %% [markdown]
# ### Compute the by-HRU counts of simulation monthly mean values that are within baseline bounds

# %%
# Make sure the array is ordered 1..nhru for the nhm_id
hru_order = list(range(1, 109952))

ss = sim_month_mean[sim_var].loc[:, hru_order].data
bb_min = baseline_df[min_var].loc[:, hru_order].data
bb_max = baseline_df[max_var].loc[:, hru_order].data

monthly_within_count = np.logical_and(ss >= bb_min, ss <= bb_max).sum(axis=0)
monthly_month_count = ss.shape[0]

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
# ### Add the monthly runoff climatology as percent within bounds

# %%
new_param = 'fit_runoff_monclim'
# pdb.parameters.remove(new_param)

pdb.parameters.add(name=new_param, datatype=2, minimum=0.0, maximum=100.0, default=0.0,
                   description='Percent of monthly climatology values within baseline bounds')

tparam = pdb.parameters.get(new_param)

tparam.dimensions.add('nhru', size=clim_within_count.shape[0])

pct = (clim_within_count / clim_month_count) * 100.0
tparam.data = pct

# %% [markdown]
# ### Add the monthly mean as percent within bounds

# %%
new_param = 'fit_runoff_monthly'

pdb.parameters.add(name=new_param, datatype=2, minimum=0.0, maximum=100.0, default=0.0,
                   description='Percent of monthly mean values within baseline bounds')

tparam = pdb.parameters.get(new_param)

tparam.dimensions.add('nhru', size=monthly_within_count.shape[0])

pct = (monthly_within_count / monthly_month_count) * 100.0
tparam.data = pct

# %%

# %%
cparam = 'fit_runoff_monclim'
use_orig_colors = False

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
