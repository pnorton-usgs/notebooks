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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import matplotlib.pyplot as plt

from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.plot_helpers import get_projection

# %%

base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10'
work_dir = f'{base_dir}/paramdb_v10_daymet_CONUS'

# base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/byHRU_test1'
# work_dir = f'{base_dir}/testDb'
out_dir = f'{base_dir}/plots'
# ------------------------------------------

# HRU polygons
hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp'
hru_layer_name = None
hru_shape_key='hru_id_nat'



# %%
pdb = ParamDb(paramdb_dir=work_dir, verbose=True, verify=True)

# %%
# Load the HRU shapes
pdb.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)

# %%
# pdb.parameters.plot('hru_area', output_dir=None, linewidth=6.0, 
#                     facecolor='snow', edgecolor='whitesmoke', 
#                     vary_color=True, vary_width=True, cmap='viridis')

# pdb.parameters.plot('dday_intcp', output_dir='/Users/pnorton/tmp/crap', cmap='tab20')
pdb.parameters.plot('carea_max', cmap='viridis')

# %%
# calib_params = {'adjmix_rain': None, 
#                 'carea_max': 'gist_stern', 
#                 'cecn_coef': None, 
#                 'dprst_frac': 'tab20c', 
#                 'emis_noppt': 'nipy_spectral', 
#                 'fastcoef_lin': 'nipy_spectral', 
#                 'freeh2o_cap': 'rainbow', 
#                 'gwflow_coef': 'gist_ncar', 
#                 'jh_coef': 'tab20', 
#                 'mann_n': None, 
#                 'potet_sublim': None, 
#                 'rad_trncf': None, 
#                 'radmax': None, 
#                 'rain_cbh_adj': None, 
#                 'slowcoef_sq': 'rainbow', 
#                 'smidx_coef': 'tab20c', 
#                 'smidx_exp': 'tab20c', 
#                 'snarea_curve': None, 
#                 'snowinfil_max': 'terrain', 
#                 'snow_cbh_adj': None, 
#                 'snow_intcp': 'terrain', 
#                 'soil2gw_max': 'terrain_r', 
#                 'soil_moist_max': 'terrain_r', 
#                 'soil_rechr_max_frac': 'terrain_r', 
#                 'srain_intcp': 'terrain', 
#                 'ssr2gw_exp': 'terrain_r', 
#                 'ssr2gw_rate': 'terrain_r', 
#                 'tmax_allrain_offset': None, 
#                 'tmax_allsnow': None, 
#                 'tmax_cbh_adj': None, 
#                 'tmin_cbh_adj': None, 
#                 'wrain_intcp': 'terrain'}
calib_params = {'carea_max': 'gist_stern', 
                'fastcoef_lin': 'nipy_spectral', 
                'freeh2o_cap': 'rainbow', 
                'gwflow_coef': 'gist_ncar', 
                'jh_coef': 'tab20', 
                'rad_trncf': None, 
                'radmax': None, 
                'rain_cbh_adj': None, 
                'slowcoef_sq': 'rainbow', 
                'smidx_coef': 'tab20c', 
                'smidx_exp': 'tab20c', 
                'snow_cbh_adj': None, 
                'soil2gw_max': 'terrain_r', 
                'soil_moist_max': 'terrain_r', 
                'soil_rechr_max_frac': 'terrain_r', 
                'srain_intcp': 'terrain', 
                'tmax_allrain_offset': None, 
                'tmax_allsnow': None, 
                'tmax_cbh_adj': None, 
                'tmin_cbh_adj': None}
for cparam, col in calib_params.items():
    print(f'Plotting {cparam}')
    pdb.parameters.plot(cparam, output_dir=out_dir, cmap=col)

# %%
sc = pdb.parameters.get_dataframe('snarea_curve')

# %%
sc.head()

# %%
sc

# %%
# Plot the snow depletion curves
crv = pdb.parameters.get_dataframe('snarea_curve')
# crv.T.plot(cmap = plt.cm.get_cmap('tab10', crv.shape[0]))
crv.T.iloc[:, 10].plot(cmap='tab10')

# %%
crv.shape[0]

# %%
crv.T.iloc[:, 1]

# %%
crv.T

# %%
