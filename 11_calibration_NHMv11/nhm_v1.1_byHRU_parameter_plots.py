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

# %%
# v1.1
cal_name = 'A_master_20210707'

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/gridmet_byHRU'
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
pdb = ParameterNetCDF(filename=filename, verbose=False, verify=True)

pdb.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)
pdb.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %%
# Calibration parameters not included: K_coef, poi_gage_id, snarea_curve
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
aa = pdb.parameters['segment_type']

# %%
aa.stats()

# %%
cparam = 'segment_type'
col = 'gist_earth_r'
pdb.parameters.plot(cparam, limits='absolute', linewidth=1.0, 
                    facecolor='snow', edgecolor='whitesmoke', var_width=True, 
                    vary_color=True, cmap=col)

# %%
import numpy as np

# unique, counts = np.unique(x, return_counts=True)
unique, counts = np.unique(aa.data, return_counts=True)

# %%
unique

# %%
counts

# %%
print(np.asarray((unique, counts)).T)

# %%
