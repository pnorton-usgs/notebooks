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
import glob
import numpy as np
import os
import pandas as pd
import sys

from collections import Counter

from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ControlFile import ControlFile
from pyPRMS.ValidParams import ValidParams

# %%
base_dir = '/Users/pnorton/Projects/National_Hydrology_Model'
paramdb_dir = f'{base_dir}/calibrations/NHMv10/byHRU_test1/paramdb_v10_daymet_CONUS'

workdir = f'{base_dir}/src/code_lauren/CALIBRATION_RESULTS_JUNE2019/Daymet_byHRU_musk'
calibdir = f'{workdir}/final_params'
outdir = f'{base_dir}/datasets/paramdb_v10/parameters'
# paramfile = f'{calibdir}/FINALparams_0000'

hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp'
hru_layer_name = None
hru_shape_key='hru_id_nat'

seg_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/GF_nat_reg_lcc/nsegmentNationalIdentifier.shp'
seg_layer_name = None
seg_shape_key = 'seg_id_nat'

# %%
pdb_orig = ParamDb(paramdb_dir, verbose=True)

# %%
byHW_parameters = ['adjmix_rain', 'carea_max', 'cecn_coef', 'emis_noppt', 
                   'fastcoef_lin', 'freeh2o_cap', 'gwflow_coef', 'mann_n', 
                   'potet_sublim', 'rad_trncf', 'radmax', 'rain_cbh_adj', 
                   'slowcoef_sq', 'smidx_coef', 'smidx_exp', 'snow_cbh_adj', 
                   'snowinfil_max', 'soil2gw_max', 'soil_moist_max', 
                   'soil_rechr_max_frac', 'ssr2gw_exp', 'ssr2gw_rate', 
                   'tmax_allrain_offset', 'tmax_allsnow', 'tmax_cbh_adj', 'tmin_cbh_adj']
print(len(byHW_parameters))

# %%
# pf = ParameterFile(filename=f'{calibdir}/FINALparams_0000', verbose=True, verify=True)

for chw in range(7265):
#     print(f'{chw:04d}')
    cfile = f'{calibdir}/FINALparams_{chw:04d}'
    sys.stdout.write(f'\rUpdating HW {chw:04d}:    ')
    
    try:
        pf = ParameterFile(filename=cfile, verbose=False, verify=True)
        nhm_id = pf.parameters.get('nhm_id').data
#         cparam = 'tmax_cbh_adj'

        for cparam in byHW_parameters:
            for cid, vals in zip (nhm_id, pf.parameters.get(cparam).data):
                pdb_orig.parameters.update_element(cparam, cid, vals)
    except FileNotFoundError:
        print(f'\nMissing headwater: FINALparams_{chw:04d}')
    
# pdb_orig.parameters.update_element(cparam, chru, pfile_hru.parameters[cparam].data[0])

# %%
pdb_orig.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)
pdb_orig.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %%
pdb_orig.parameters.plot('mann_n', use_drange=False, cmap='bwr')

# %%



# %%
pdb_orig.write_netcdf(filename=f'{outdir}/C_byHRU_musk_param_alt1.nc')

# %%
svar = pf.parameters.get('tmax_cbh_adj').data
nhmid = pf.parameters.get('nhm_id').data

# %%
for xx, yy in zip(nhmid, svar):
    print(xx, yy)

# %%
svar[0, :]

# %%
