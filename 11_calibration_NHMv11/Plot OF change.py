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
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import numpy as np
import pandas as pd

from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile

# %%

# %%

# %%

# %%

# %%
headwater = '3020'
workdir = '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHW/20210927_gm_byHW'
ofs_file = f'{workdir}/objfun/objfun_{headwater}'

paramdb_dir = f'{workdir}/paramdb_v11_gridmet_CONUS'

hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
hru_layer_name = 'nhruv11_sim30'
hru_shape_key='nhru_v11'


# %%
pdb_orig = ParamDb(paramdb_dir, verbose=True)

nhru = pdb_orig.dimensions.get('nhru').size

# %%
pdb_orig.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)

# %%
cfile = f'{workdir}/FINALparams/FINALparams_{headwater}'

pf = ParameterFile(filename=cfile, verbose=False, verify=True)
nhm_id = pf.parameters.get('nhm_id').data
nhm_seg = pf.parameters.get('nhm_seg').data

# %%

# %%
def_data = np.zeros((nhru))
def_data[:] = -999.0

# %%
pdb_orig.parameters.add(name='of_prms_start', datatype=2, units=None, description='Starting OF', help='Starting OF',
                        default=-999.0)
pdb_orig.parameters.get('of_prms_start').dimensions.add(name='nhru', size=nhru)

tp = pdb_orig.parameters.get('of_prms_start')
tp.data = def_data

# %%

# %%
col_names = ['HW', 'step', 'of_mth_range', 'of_mnmth_range', 'of_mnmth_median', 
             'of_daily_range', 'of_high_daily_median', 'of_low_daily_median',
             'of_prms', 'of_run', 'of_aet', 'of_sca', 'of_rch', 'of_som']

# NOTE: Must use Int64 dtype for nullable-integer fields
col_types = [str, int, float, float, float, float, float, float, float, float, float, float, float, float]

cols = dict(zip(col_names, col_types))

df = pd.read_csv(ofs_file, sep='\s+', skipinitialspace=True, dtype=cols, na_values=['**************'])        

x_vars = df.columns.tolist()[2:]

# %%
first_df = df.iloc[0]
last_df = df.iloc[-1]

# %%
chg_df = df.iloc[-1, 2:] - df.iloc[0, 2:]
chg_df

# %%
chg_df.loc['of_prms']

# %%
last_df

# %%

# %%

# %%
for cid in nhm_id:
    pdb_orig.parameters.update_element('of_prms_start', cid, first_df['of_prms'])

# %%
pdb_orig.parameters.plot('of_prms_start', limits='valid', 
                         linewidth=0.0, edgecolor='whitesmoke', mask_defaults='darkgrey', cmap='YlGnBu')

# %%

# %%
col_names = ['HW', 'step', 'of_mth_range', 'of_mnmth_range', 'of_mnmth_median', 
             'of_daily_range', 'of_high_daily_median', 'of_low_daily_median',
             'of_prms', 'of_run', 'of_aet', 'of_sca', 'of_rch', 'of_som']

# NOTE: Must use Int64 dtype for nullable-integer fields
col_types = [str, int, float, float, float, float, float, float, float, float, float, float, float, float]

cols = dict(zip(col_names, col_types))

df = pd.read_csv(f'{workdir}/final_of.txt', sep='\s+', skipinitialspace=True, dtype=cols, na_values=['**************'])        

x_vars = df.columns.tolist()[2:]

# %%
# df.loc[:, df.notna().sum(axis=0)
of_nan_df = df[df.isna().any(axis=1)]

# %%
of_nan_df.to_csv(f'{workdir}/of_prms_nan.csv', index=False)

# %%

# %%

# %%

# %%

# %%

# %%

# %%
