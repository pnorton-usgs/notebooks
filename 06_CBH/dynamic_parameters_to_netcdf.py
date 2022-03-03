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
#     display_name: Python [conda env:bandit_py3]
#     language: python
#     name: conda-env-bandit_py3-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import datetime
import pandas as pd
import netCDF4
import numpy as np

from pyPRMS.prms_helpers import dparse
from Bandit.dynamic_parameters import DynamicParameters

# %%
DATETIME_INDEX_COLS = [0, 1, 2]

param = 'wrain_intcp'
param_name = 'dyn_{}'.format(param)

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/dynamic_params/2019-08_dynamic_parameters'
cbh_file = '{}/{}.param'.format(workdir, param_name)

dynparam_datatype = {'cov_type': 'i4',
                     'hru_percent_imperv': 'f4',
                     'snow_intcp': 'f4',
                     'srain_intcp': 'f4',
                     'wrain_intcp': 'f4'}
dynparam_units = {'cov_type': 'none',
                  'hru_percent_imperv': 'fraction',
                  'snow_intcp': 'inches',
                  'srain_intcp': 'inches',
                  'wrain_intcp': 'inches'}

# st_date = datetime.datetime(1931, 1, 1)
# en_date = datetime.datetime(2100, 12, 31)

# idx_retrieve = {15: 15, 16: 16, 17: 17}

# load_cols = list(CBH_INDEX_COLS)
# load_cols.extend([xx+5 for xx in idx_retrieve.keys()])

# %% [markdown]
# ### Read ASCII dynamic parameter file and save to netCDF

# %%
df = pd.read_csv(cbh_file, sep=' ', skipinitialspace=True, comment='#',
                 skiprows=1, engine='c', memory_map=True, dtype=dynparam_datatype[param],
                 date_parser=dparse, parse_dates={'time': DATETIME_INDEX_COLS},
                 index_col='time', na_values=[-99.0, -999.0, 'NaN', 'inf'])

# Create a netCDF file for the dynamic parameter data
nco = netCDF4.Dataset('{}/{}.nc'.format(workdir, param_name), 'w', clobber=True)
nco.createDimension('hru', len(df.columns))
nco.createDimension('time', None)

timeo = nco.createVariable('time', 'f4', ('time'))
hruo = nco.createVariable('hru', 'i4', ('hru'))

varo = nco.createVariable(param, dynparam_datatype[param], ('time', 'hru'), zlib=True)

# NOTE: Do not use _FillValue for the dynamic parameter files. The xarray library that is used to read these
#       netcdf files will automatically create float32 arrays for integers with the _FillValue attribute set.
#                           fill_value=netCDF4.default_fillvals[dynparam_datatype[param]], zlib=True)

nco.setncattr('Description', 'Dynamic {} parameter by HRU'.format(param))
# nco.setncattr('Bandit_version', __version__)
# nco.setncattr('NHM_version', nhmparamdb_revision)

timeo.calendar = 'standard'
# timeo.bounds = 'time_bnds'
timeo.units = 'days since 1900-01-01 00:00:00'

hruo.long_name = 'Hydrologic Response Unit ID (HRU)'

varo.long_name = 'Dynamic {}'.format(param)
varo.units = dynparam_units[param]

# Write the HRU ids
hruo[:] = df.columns.values

timeo[:] = netCDF4.date2num(df.index.tolist(), units='days since 1900-01-01 00:00:00',
                            calendar='standard')

# Write the CBH values
varo[:, :] = df.values

nco.close()

# %% [markdown]
# ### Read the netCDF dynamic parameter file

# %%
st_date = datetime.datetime(1980, 1, 1)
en_date = datetime.datetime(1981, 12, 31)
nhm_hrus = [1,5,2,3]

mydyn = DynamicParameters('{}/{}.nc'.format(workdir, param_name), param, st_date, en_date, nhm_hrus)

# %%
mydyn.read()

# %%

# %%


mydyn.read_netcdf()

# %%
mydyn.data

# %%
out_order = [kk for kk in nhm_hrus]
for cc in ['day', 'month', 'year']:
    out_order.insert(0, cc)

header = ' '.join(map(str, out_order))

# Output ASCII files
out_ascii = open('crap.param', 'w')
out_ascii.write('{}\n'.format(param))
out_ascii.write('{}\n'.format(header))
out_ascii.write('####\n')
mydyn.data.to_csv(out_ascii, columns=out_order, na_rep='-999', 
                    sep=' ', index=False, header=False, encoding=None, chunksize=50)
out_ascii.close()

# %%
mydyn.ds

# %%
