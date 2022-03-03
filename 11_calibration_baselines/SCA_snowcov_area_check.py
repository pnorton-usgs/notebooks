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

sim_dir = '/Volumes/parker_rocks/NHM_output/netcdf'
# sim_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/Maurer/byHRU_musk'

baseline_file = f'{baseline_dir}/baselines_{baseline_var}.nc'

sim_var = 'snowcov_area'
sim_file = f'{sim_dir}/{sim_var}.nc'

# %%
st_date = datetime.datetime(2000,1,1)
en_date = datetime.datetime(2010,12,31)

# %% [markdown]
# ### Read the baseline information

# %%
# Baseline recharge has daily timestep
baseline_df = xr.open_dataset(baseline_file)
baseline_df = baseline_df[[sca_var, ci_var, 'hru']].sel(time=slice(st_date, en_date))
# baseline_df

# Remove July and August from the dataset
# baseline_restr = baseline_df.sel(time=baseline_df.time.dt.month.isin([1, 2, 3, 4, 5, 6, 9, 10, 11, 12]))
baseline_restr = baseline_df

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

SCAmask = (msk_num_ann > 1) & (msk_SCAmax > 0.5) & (msk_num_obs > 9)


# %%
# msk_num_yr = 
msk_ann = sca_obs.resample(time='1AS').mean()
kk = (msk_ann > 0).sum(axis=0)
np.count_nonzero(kk > 0)

# %%
#         if (CI >= 70.0) then
#           if (all([7, 8] /= nm)) then
# #             ! Compute obsSCA for all months except July and August
#             obsSCA(1, n) = (CI / 100.0) * obs
#             obsSCA(2, n) = obsSCA(1, n) + (100.0 - CI) / 100.0
#           end if
#         else
#           obsSCA(1, n) = -888.0
#           obsSCA(2, n) = -888.0
#         end if

# Compute lower and upper SCA values based on confidence interval
# threshold = 70.0
# ci_pct = baseline_restr['SCAci'].where(baseline_restr['SCAci'] >= threshold)
# ci_pct /= 100.0

# sca_obs = baseline_restr['SCA'].where(baseline_restr['SCA'] > 0.0)
# sca_obs = baseline_restr['SCA']
baseline_SCAmin = (ci_pct * sca_obs).where(SCAmask)
# baseline_SCAmin = baseline_SCAmin.where(~np.isnan(ci_pct))

baseline_SCAmax = (baseline_SCAmin + (1.0 - ci_pct)).where(SCAmask)
# baseline_SCAmax = baseline_SCAmax.where(~np.isnan(ci_pct))

# baseline_SCAmax = baseline_SCAmax.where(baseline_SCAmax > 0.005)

# %%
print(f'SCAci: min={np.nanmin(ci_pct)}, max={np.nanmax(ci_pct)}')
print(f'SCAmin: min={np.nanmin(baseline_SCAmin)}, max={np.nanmax(baseline_SCAmin)}')
print(f'SCAmax: min={np.nanmin(baseline_SCAmax)}, max={np.nanmax(baseline_SCAmax)}')

# %%
ci_pct[0:364, 82974]

# %%

# %% [markdown]
# ### Read the model output from the simulation

# %%
# Model output has daily timestep
sim_df = xr.open_dataset(sim_file)

# NOTE: Next line needed if nhm_id variable exists in netcdf file
sim_df = sim_df.assign_coords(nhru=(sim_df['nhm_id']))

# sim_df = sim_df[[sim_var]].sel(time=slice(st_date, en_date))
sim_df = sim_df.sel(time=slice(st_date, en_date))

sim_restr = sim_df.sel(time=sim_df.time.dt.month.isin([8]))
# sim_restr = sim_df.sel(time=sim_df.time.dt.month.isin([1, 2, 3, 4, 5, 6, 9, 10, 11, 12]))
# sim_df

# %%
sim_restr

# %%
mon_tot = sim_restr.resample(time='1MS').sum()
# mon_tot = sim_restr['snowcov_area'].resample(time='1MS').sum()
mon_tot

# %%
mon_clim = mon_tot.groupby('time.month').sum('time')
mon_clim

# %%
hru_order = list(range(1, 109952))

ss = mon_clim[sim_var].loc[:, hru_order].data
np.nanmax(ss[7,:])

# %% [markdown]
# # Plot the results

# %%
import matplotlib.pyplot as plt

from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ParameterNetCDF import ParameterNetCDF

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/paramdb_daymet_CONUS'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/daymet_fix_work'
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
new_param = 'snowcov_area'
pdb.parameters.remove(new_param)

pdb.parameters.add(name=new_param, datatype=2, minimum=0.0, maximum=500.0, default=0.0,
                   description='Snow covered area')

tparam = pdb.parameters.get(new_param)

tparam.dimensions.add('nhru', size=ss.shape[1])

tparam.data = ss[7, :]

# %%
cparam = 'snowcov_area'

cmap = plt.cm.get_cmap('RdYlBu', 5)

# Colormap using Lauren's original colors for the T&M
# from matplotlib.colors import LinearSegmentedColormap

# colors = [(0, 0, 1), (0.13, 1, 1), (0.13, 1, 0.02), (0.99, 0.58, 0.04), (0.98, 0, 1)]
# n_bin = 5
# cmap_name = 'lauren_col'
# cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)

# Plot
# cmap.set_bad('black', 1.)
pdb.parameters.plot(cparam, limits='valid', linewidth=0.0, edgecolor='whitesmoke', mask_defaults='white', cmap=cmap)

# %%

# %%
