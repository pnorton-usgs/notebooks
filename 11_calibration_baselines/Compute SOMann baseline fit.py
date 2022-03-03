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

baseline_var = 'SOMann'
max_var = baseline_info[baseline_var]['max']
min_var = baseline_info[baseline_var]['min']

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
st_date = datetime.datetime(1982,1,1)
en_date = datetime.datetime(2010,12,31)

num_yrs = (en_date.year - st_date.year) + 1

# %% [markdown]
# ### Read the baseline information

# %%
baseline_df = xr.open_dataset(baseline_file)
baseline_df = baseline_df.assign_coords(nhru=(baseline_df.nhm_id))
baseline_df = baseline_df[[max_var, min_var, 'nhru']].sel(time=slice(st_date, en_date))

# %% [markdown]
# ### Read the model output from the simulation

# %%
# Model output has daily timestep
sim_df = xr.open_dataset(sim_file)

# NOTE: Next line needed if nhm_id variable exists in netcdf file
sim_df = sim_df.assign_coords(nhru=(sim_df.nhm_id))

# List ordered 1..nhru
hru_order = list(range(1, sim_df['nhru'].size + 1))

sim_df = sim_df[[sim_var]].sel(time=slice(st_date, en_date))
# sim_df

# %% [markdown]
# #### Compute monthly total from daily model output

# %%
# Compute the monthly total from the simulation daily output
# 1MS - resample to monthly means from the start of each month
# sim_month_total = sim_df.resample(time='1MS').sum()

# %% [markdown]
# #### Compute the annual total from daily model output

# %%
# sim_ann_total = sim_df.resample(time='1AS').sum()
# sim_ann_total = sim_month_total.resample(time='1AS').sum()
sim_ann_total = sim_df.resample(time='1AS').sum()

# %%
sim_month_total[sim_var].data.max(axis=0).shape

# %%
ss = sim_ann_total[sim_var].loc[:, hru_order].data
ss_min = ss.min(axis=0)
ss_max = ss.max(axis=0)

# Where ss_max == ss_min increment ss_max by a small amount
ss_max = np.where(ss_max == ss_min, ss_max + (0.01 * ss_max), ss_max)
ss_max = np.where(ss_max == ss_min, ss_max + 0.0001, ss_max)

# Normalize the simulation values 0.0 to 1.0 for comparison to baseline
ss = (ss - ss_min) / (ss_max - ss_min)

bb_min = baseline_df[min_var].loc[:, hru_order].data
bb_max = baseline_df[max_var].loc[:, hru_order].data

ann_within_count = ((ss >= bb_min) & (ss <= bb_max)).sum(axis=0)
ann_count = (~np.isnan(ss)).sum(axis=0)

# %%
ann_count.max()

# %% [markdown]
# ### Plot the results

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

# %%
cparam = 'fit_SOM_anntot'
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
chru = 1

print(f'month_within_count: {month_within_count[chru]}, total_count: {month_count[chru]}')
print(f'ann_within_count: {ann_within_count[chru]}, total_count: {ann_count[chru]}')

# %%
sim_month_total[sim_var].data.shape

# %%
sim_ann_total[sim_var].data.shape

# %%
