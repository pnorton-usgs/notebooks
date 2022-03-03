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
from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ParameterNetCDF import ParameterNetCDF
# from pyPRMS.plot_helpers import get_projection

# %%
# V1.0
cal_name = 'B_byHRU_param'
# # cal_name = 'B_minus_C'

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/maurer_fix_work'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/daymet_fix_work'
filename = f'{workdir}/{cal_name}.nc'

plot_dir = f'{workdir}/cal_plots_{cal_name}'

hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp'
hru_layer_name = None
hru_shape_key='hru_id_nat'

seg_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/GF_nat_reg_lcc/nsegmentNationalIdentifier.shp'
seg_layer_name = None
seg_shape_key = 'seg_id_nat'

# %%
# pdb = ParameterFile(filename, verbose=True, verify=True)
pdb = ParameterNetCDF(filename=filename, verbose=False, verify=True)
# pdb = ParamDb(paramdb_dir=work_dir, verbose=True, verify=True)

pdb.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)
pdb.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %%

# %%
# pdb.parameters.plot('hru_area', output_dir=None, linewidth=6.0, 
#                     facecolor='snow', edgecolor='whitesmoke', 
#                     vary_color=True, vary_width=True, cmap='viridis')

# pdb.parameters.plot('radmax', output_dir='/Users/pnorton/tmp', limits='valid', linewidth=0.0, edgecolor='whitesmoke', 
#                     mask_defaults='darkgrey', cmap='Greens')

# pdb.parameters.plot('mann_n', limits='absolute', linewidth=1.0, 
#                     facecolor='snow', edgecolor='whitesmoke', var_width=True, 
#                     vary_color=True, cmap='plasma_r')

# %%
# byHRU calibration parameters
calib_params = {'carea_max': 'terrain_r', 
                'fastcoef_lin': 'Blues', 
                'freeh2o_cap': 'Blues', 
                'gwflow_coef': 'YlGnBu', 
                'jh_coef': 'YlOrRd', 
                'rad_trncf': 'Oranges', 
                'rain_cbh_adj': 'YlGnBu', 
                'slowcoef_sq': 'GnBu', 
                'smidx_coef': 'gist_heat_r', 
                'smidx_exp': 'gist_heat_r', 
                'snarea_thresh': 'YlGnBu', 
                'snow_cbh_adj': 'YlGnBu', 
                'soil2gw_max': 'gist_earth_r', 
                'soil_moist_max': 'gist_earth_r', 
                'soil_rechr_max_frac': 'gist_earth_r', 
                'tmax_allrain_offset': 'YlOrRd', 
                'tmax_allsnow': 'plasma', 
                'tmax_cbh_adj': 'RdBu_r', 
                'tmin_cbh_adj': 'RdBu_r'}

# byHW calibration parameters
# calib_params = {'adjmix_rain': 'Greens', 
#                 'carea_max': 'terrain_r', 
#                 'cecn_coef': 'gnuplot', 
#                 'emis_noppt': 'Reds', 
#                 'fastcoef_lin': 'Blues', 
#                 'freeh2o_cap': 'Blues', 
#                 'gwflow_coef': 'YlGnBu', 
#                 'jh_coef': 'YlOrRd', 
#                 'mann_n': 'plasma_r', 
#                 'potet_sublim': 'viridis', 
#                 'rad_trncf': 'Oranges', 
#                 'radmax': 'YlOrRd', 
#                 'rain_cbh_adj': 'YlGnBu', 
#                 'slowcoef_sq': 'GnBu', 
#                 'smidx_coef': 'gist_heat_r', 
#                 'smidx_exp': 'gist_heat_r', 
#                 'snowinfil_max': 'YlGnBu', 
#                 'snow_cbh_adj': 'YlGnBu', 
#                 'soil2gw_max': 'gist_earth_r', 
#                 'soil_moist_max': 'gist_earth_r', 
#                 'soil_rechr_max_frac': 'gist_earth_r', 
#                 'ssr2gw_exp': 'gist_earth_r', 
#                 'ssr2gw_rate': 'gist_earth_r', 
#                 'tmax_allrain_offset': 'YlOrRd', 
#                 'tmax_allsnow': 'plasma', 
#                 'tmax_cbh_adj': 'RdBu_r', 
#                 'tmin_cbh_adj': 'RdBu_r'}

for cparam, col in calib_params.items():
    print(f'Plotting {cparam}')
    
    if cparam == 'mann_n':
        pdb.parameters.plot(cparam, output_dir=plot_dir, limits='absolute', linewidth=1.0, 
                            facecolor='snow', edgecolor='whitesmoke', var_width=True, 
                            vary_color=True, cmap=col)
    else:
        pdb.parameters.plot(cparam, output_dir=plot_dir, limits='valid', 
                            linewidth=0.0, edgecolor='whitesmoke', mask_defaults='darkgrey', cmap=col)
    

# %%
# pdb.parameters.plot('seg_width', output_dir=None, linewidth=6.0, 
#                     facecolor='snow', edgecolor='whitesmoke', 
#                     vary_color=True, vary_width=True, cmap='tab20')

# %%
pdb.parameters['jh_coef'].stats()

# %%

# %%
import copy

# %%
