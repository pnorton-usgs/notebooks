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

# from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterNetCDF import ParameterNetCDF
# from pyPRMS.plot_helpers import get_projection

# %%
# V1.0
cal_name = 'G_new_minus_F_new'
# # cal_name = 'B_minus_C'

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/maurer_fix_work'
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

pnc.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)
pnc.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %%

# %%

# %%
# byHW calibration parameters
cal_params = ['adjmix_rain', 'carea_max', 'cecn_coef', 'emis_noppt', 'fastcoef_lin', 
              'freeh2o_cap', 'gwflow_coef', 'jh_coef', 'mann_n', 'potet_sublim', 
              'rad_trncf', 'radmax', 'rain_cbh_adj', 'slowcoef_sq', 'smidx_coef', 'smidx_exp', 
              'snow_cbh_adj', 'snowinfil_max', 'soil2gw_max', 'soil_moist_max', 
              'soil_rechr_max_frac', 'ssr2gw_exp', 'ssr2gw_rate', 'tmax_allrain_offset', 
              'tmax_allsnow', 'tmax_cbh_adj', 'tmin_cbh_adj']

for cparam in cal_params:
    print(f'Plotting {cparam}')
#     pnc.parameters.plot(cparam, output_dir=plot_dir, limits='absolute', cmap='bwr')
    pnc.parameters.plot(cparam, output_dir=plot_dir, limits='centered', cmap='bwr')

# %%
pnc.parameters['adjmix_rain'].data.shape

# %%
import pandas as pd

pd.set_option('display.max_rows', None)
param_stats = []

for pp in pnc.parameters.values():
    param_stats.append(pp.stats())
    
df = pd.DataFrame.from_records(param_stats, columns=['name', 'min', 'max', 'mean', 'median'])

# %%
df.head(126)

# %%
aa = pnc.parameters.get_dataframe('mann_n')
# aa.head(10)

for xx in [3705, 4226, 3703]:
    print(aa.loc[xx].tolist())

# %%
