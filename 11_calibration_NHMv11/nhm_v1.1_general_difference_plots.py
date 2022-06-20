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
import matplotlib.pyplot as plt

from pyPRMS.ParameterNetCDF import ParameterNetCDF

# %%
import warnings
warnings.filterwarnings('ignore')

# %%
# v1.1
#           A_minus_B_20210331.nc

base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/FY2022_projects/IWAAs/20220418_CONUS_check_calib_param_files/byHRU'

# filename = f'{base_dir}/byHW_minus_byHRU.nc'
filename = f'{base_dir}/new_minus_old.nc'

plot_dir = f'{base_dir}/diff_plots'

hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
hru_layer_name = 'nhruv11_sim30'
hru_shape_key='nhru_v11'

# Segment lines
seg_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
seg_layer_name = 'nsegment_v11'
seg_shape_key = 'nsegment_v11'

# %%
pnc = ParameterNetCDF(filename=filename, verbose=True, verify=True)

pnc.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)
pnc.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %%
# Calibration parameters not included: K_coef, poi_gage_id, snarea_curve
cal_params = ['carea_max', 'fastcoef_lin', 'freeh2o_cap', 'gwflow_coef'] #, 
              # 'jh_coef', 'rad_trncf', 'rain_cbh_adj', 'slowcoef_sq', 
              # 'smidx_coef', 'smidx_exp','snarea_thresh', 'snow_cbh_adj', 
              # 'soil2gw_max', 'soil_moist_max', 'soil_rechr_max_frac', 
              # 'tmax_allrain_offset', 'tmax_allsnow', 
              # 'tmax_cbh_adj', 'tmin_cbh_adj']

for cparam in cal_params:
    print(f'Plotting {cparam}')
#     pnc.parameters.plot(cparam, output_dir=plot_dir, limits='absolute', cmap='bwr')
    pnc.parameters.plot(cparam, output_dir=plot_dir, limits='centered', cmap='bwr')

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
import xarray as xr

df = xr.open_dataset(filename, mask_and_scale=False, decode_timedelta=False)

# %%
aa = df['gw_tau']

# %%
aa.T.dims

# %%
len(aa.T.dims)

# %%
aa.values

# %%
df.info()

# %%

# %%
