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

# %%
# import IPython

# IPython.Application.instance().kernel.do_shutdown(True)

# # %%javascript
# IPython.notebook.kernel.restart()

# %%
import matplotlib.pyplot as plt

from pyPRMS.ParameterNetCDF import ParameterNetCDF

# %%
# v1.1
cal_name = 'B_20211027_minus_A'

workdir = '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHW/20211027_gm_byHW/calibration_diff_work'
filename = f'{workdir}/{cal_name}.nc'

plot_dir = f'{workdir}/cal_plots_{cal_name}'

hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
hru_layer_name = 'nhruv11_sim30'
hru_shape_key='nhru_v11'

# Segment lines
seg_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
seg_layer_name = 'nsegment_v11'
seg_shape_key = 'nsegment_v11'

# %%
pnc = ParameterNetCDF(filename=filename, verbose=False, verify=True)

pnc.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)
pnc.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %%
# Calibration parameters not included: K_coef, poi_gage_id, snarea_curve
cal_params = ['adjmix_rain', 'carea_max', 'cecn_coef', 'emis_noppt', 
              'fastcoef_lin', 'freeh2o_cap', 'gwflow_coef', 'mann_n',
              'potet_sublim', 'rad_trncf', 'radmax', 'rain_cbh_adj', 'slowcoef_sq', 
              'smidx_coef', 'smidx_exp','snowinfil_max', 'snow_cbh_adj', 
              'soil2gw_max', 'soil_moist_max', 'soil_rechr_max_frac', 
              'ssr2gw_exp', 'ssr2gw_rate', 
              'tmax_allrain_offset', 'tmax_allsnow', 
              'tmax_cbh_adj', 'tmin_cbh_adj']

for cparam in cal_params:
    print(f'Plotting {cparam}')
    if cparam == 'mann_n':
        pnc.parameters.plot(cparam, output_dir=plot_dir, limits='centered', linewidth=1.0, 
                            facecolor='snow', edgecolor='whitesmoke', var_width=True,
                            vary_color=True, cmap='bwr')
    else:
        # pnc.parameters.plot(cparam, output_dir=plot_dir, limits='absolute', cmap='bwr')
        pnc.parameters.plot(cparam, output_dir=plot_dir, limits='centered', cmap='bwr')

# %%
