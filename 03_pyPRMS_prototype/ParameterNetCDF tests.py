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
from pyPRMS.ParameterNetCDF import ParameterNetCDF

# %%


# V1.1 gridmet byHRU calibration results
# /Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/gridmet_byHRU
# basedir = '/Users/pnorton/Projects/National_Hydrology_Model'
# workdir = f'{basedir}/calibrations/NHMv11/gridmet_byHRU'
# filename = f'{workdir}/A_minus_B.nc'

# plot_dir = f'{workdir}/diff_plots'

# HRU polygons
# hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
# hru_layer_name = 'nhruv11_sim30'
# hru_shape_key = 'nhru_v11'


# V1.0
cal_name = 'C_minus_D_new_byHW_musk_obs'
# # cal_name = 'B_minus_C'

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/parameters'
filename = f'{workdir}/{cal_name}.nc'

plot_dir = f'{workdir}/cal_plots_{cal_name}'

hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp'
hru_layer_name = None
hru_shape_key='hru_id_nat'

seg_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/GF_nat_reg_lcc/nsegmentNationalIdentifier.shp'
seg_layer_name = None
seg_shape_key = 'seg_id_nat'

# %%
pnc = ParameterNetCDF(filename=filename, verbose=False, verify=True)

# %%
pnc.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)
pnc.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %%
# pnc.parameters.plot('tmax_cbh_adj', use_drange=True, cmap='bwr')

# %%
# Calibration parameters not included: K_coef, poi_gage_id, snarea_curve
cal_params = ['adjmix_rain', 'carea_max', 'cecn_coef', 'emis_noppt', 'fastcoef_lin', 
              'freeh2o_cap', 'gwflow_coef', 'jh_coef', 'mann_n', 'potet_sublim', 
              'rad_trncf', 'radmax', 'rain_cbh_adj', 'slowcoef_sq', 'smidx_coef', 'smidx_exp', 
              'snow_cbh_adj', 'snowinfil_max', 'soil2gw_max', 'soil_moist_max', 
              'soil_rechr_max_frac', 'ssr2gw_exp', 'ssr2gw_rate', 'tmax_allrain_offset', 
              'tmax_allsnow', 'tmax_cbh_adj', 'tmin_cbh_adj']

# %%
for cparam in cal_params:
    print(f'Plotting {cparam}')
#     pnc.parameters.plot(cparam, output_dir=plot_dir, use_drange=False, cmap='bwr')
    pnc.parameters.plot(cparam, output_dir=plot_dir, use_drange=True, cmap='bwr')

# %%
pnc.parameters.plot('mann_n', output_dir=plot_dir, use_drange=True, cmap='bwr')

# %%
aa = 'HW1476'

# %%
aa[-4:]

# %%
