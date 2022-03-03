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
from pyPRMS.plot_helpers import get_projection

# %%
calib = 'byHRU'
base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10'
work_dir = f'{base_dir}/nhmparamdb_lhay_maurer_{calib}'
out_dir = f'{base_dir}/maurer_calib_plots/maurer_{calib}'
# ------------------------------------------

base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10'
work_dir = f'{base_dir}/paramdb_v10_daymet_CONUS'
out_dir = f'{base_dir}/paramdb_daymet_CONUS_plots'

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20200415_red_river_v2'
# filename = f'{workdir}/myparam.param'

# shpfile = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
# layer_name = 'nsegment_v11'

# HRU polygons
# hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
# hru_layer_name = 'nhruv11_sim30'
# hru_shape_key = 'nhru_v11'

hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp'
hru_layer_name = None
hru_shape_key='hru_id_nat'

# Segment lines
# seg_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
# seg_layer_name = 'nsegment_v11'
# seg_shape_key = 'nsegment_v11'

# %%
# pdb = ParameterFile(filename, verbose=True, verify=True)
pdb = ParamDb(paramdb_dir=work_dir, verbose=True, verify=True)

# %%
# Load the HRU shapes
# pdb.parameters.shapefile_hrus(f'{workdir}/GIS/HRU_subset.shp', layer_name=None, shape_key='nhru_v11')
pdb.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)

# seg_shape_key = 'nsegment_v'
# pdb.parameters.shapefile_segments(f'{workdir}/GIS/Segments_subset.shp', layer_name=None, shape_key=seg_shape_key)

# pdb.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %%
# pdb.parameters.plot('hru_area', output_dir=None, linewidth=6.0, 
#                     facecolor='snow', edgecolor='whitesmoke', 
#                     vary_color=True, vary_width=True, cmap='viridis')

pdb.parameters.plot('smidx_exp', cmap='tab20')

# %%
# calib_params = ['adjmix_rain', 'carea_max', 'cecn_coef', 'dprst_frac', 'emis_noppt', 
#                 'fastcoef_lin', 'freeh2o_cap', 'gwflow_coef', 'jh_coef', 'mann_n', 
#                 'potet_sublim', 'rad_trncf', 'radmax', 'rain_cbh_adj', 'slowcoef_sq', 
#                 'smidx_coef', 'smidx_exp', 'snarea_curve', 'snowinfil_max', 'snow_cbh_adj', 
#                 'snow_intcp', 'soil2gw_max', 'soil_moist_max', 'soil_rechr_max_frac', 
#                 'srain_intcp', 'ssr2gw_exp', 'ssr2gw_rate', 'tmax_allrain_offset', 
#                 'tmax_allsnow', 'tmax_cbh_adj', 'tmin_cbh_adj', 'wrain_intcp']
calib_params = {'adjmix_rain': None, 
                'carea_max': 'gist_stern', 
                'cecn_coef': None, 
                'dprst_frac': 'tab20c', 
                'emis_noppt': 'nipy_spectral', 
                'fastcoef_lin': 'nipy_spectral', 
                'freeh2o_cap': 'rainbow', 
                'gwflow_coef': 'gist_ncar', 
                'jh_coef': 'tab20', 
                'mann_n': None, 
                'potet_sublim': None, 
                'rad_trncf': None, 
                'radmax': None, 
                'rain_cbh_adj': None, 
                'slowcoef_sq': 'rainbow', 
                'smidx_coef': 'tab20c', 
                'smidx_exp': 'tab20c', 
                'snarea_curve': None, 
                'snowinfil_max': 'terrain', 
                'snow_cbh_adj': None, 
                'snow_intcp': 'terrain', 
                'soil2gw_max': 'terrain_r', 
                'soil_moist_max': 'terrain_r', 
                'soil_rechr_max_frac': 'terrain_r', 
                'srain_intcp': 'terrain', 
                'ssr2gw_exp': 'terrain_r', 
                'ssr2gw_rate': 'terrain_r', 
                'tmax_allrain_offset': None, 
                'tmax_allsnow': None, 
                'tmax_cbh_adj': None, 
                'tmin_cbh_adj': None, 
                'wrain_intcp': 'terrain'}
for cparam, col in calib_params.items():
    print(f'Plotting {cparam}')
    pdb.parameters.plot(cparam, output_dir=out_dir, cmap=col)
    

# %%
# pdb.parameters.plot('seg_width', output_dir=None, linewidth=6.0, 
#                     facecolor='snow', edgecolor='whitesmoke', 
#                     vary_color=True, vary_width=True, cmap='tab20')

# %%
pdb.parameters['jh_coef'].stats()

# %%
pdb.parameters['jh_coef'].maximum

# %%
pdb.parameters['jh_coef'].check()

# %%
pdb.parameters['jh_coef'].stats()

# %%
pdb.parameters['jh_coef'].minimum

# %%
pdb.parameters['jh_coef'].check_values()

# %%
