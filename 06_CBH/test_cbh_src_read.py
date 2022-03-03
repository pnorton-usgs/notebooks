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
#     display_name: Python [conda env:bandit]
#     language: python
#     name: conda-env-bandit-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import datetime
import pandas as pd
import numpy as np

from pyPRMS.prms_helpers import dparse

# %%
CBH_INDEX_COLS = [0, 1, 2, 3, 4, 5]

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/tmp'
# cbh_file = f'{workdir}/prcpHRU.cbh'
cbh_file = f'{workdir}/prcp.cbh'

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/fortran/prms6/tests/conus'
# cbh_file = '{}/prcp.cbh'.format(workdir)

st_date = datetime.datetime(1980, 1, 1)
en_date = datetime.datetime(2010, 12, 31)

# idx_retrieve = {15: 15, 16: 16, 17: 17}

# load_cols = list(CBH_INDEX_COLS)
# load_cols.extend([xx+5 for xx in idx_retrieve.keys()])

# %%
# df = pd.read_csv(cbh_file, sep=' ', skipinitialspace=True, usecols=load_cols,
#                  skiprows=3, engine='c', memory_map=True,
#                  date_parser=dparse, parse_dates={'time': CBH_INDEX_COLS},
#                  index_col='time', header=None, na_values=[-99.0, -999.0, 'NaN', 'inf'])

df = pd.read_csv(cbh_file, sep=' ', skipinitialspace=True,
                 skiprows=3, engine='c', memory_map=True,
                 date_parser=dparse, parse_dates={'time': CBH_INDEX_COLS},
                 index_col='time', header=None, na_values=[-99.0, -999.0, 'NaN', 'inf'])

# %%
df.info()

# %%
ppt_max = np.max(df[6].to_numpy())

# %%
soil_moist_max = 2.4547219

smidx_max = (1.1 * soil_moist_max) + (0.5 * ppt_max)

print(smidx_max)

# %%
out_cbh = 'test_cbh.cbh'
df.to_csv(out_cbh, na_rep='-999', float_format='%0.3f', sep=' ', index=False, header=False, encoding=None, chunksize=50)

# %%
# Create a netCDF file for the CBH data
nco = netCDF4.Dataset('{}/{}.nc'.format(outdir, vv), 'w', clobber=True)
nco.createDimension('hru', len(df2.columns))
nco.createDimension('time', None)

timeo = nco.createVariable('time', 'f4', ('time'))
hruo = nco.createVariable('hru', 'i4', ('hru'))

varo = nco.createVariable(vv, 'f4', ('time', 'hru'), fill_value=NC_FILL_FLOAT, zlib=True)

nco.setncattr('Description', 'Climate by HRU')
nco.setncattr('Bandit_version', __version__)
nco.setncattr('NHM_version', nhmparamdb_revision)

timeo.calendar = 'standard'
# timeo.bounds = 'time_bnds'
timeo.units = 'days since 1980-01-01 00:00:00'

hruo.long_name = 'Hydrologic Response Unit ID (HRU)'

varo.long_name = var_desc[vv]
varo.units = var_units[vv]

# Write the HRU ids
hruo[:] = df2.columns.values

timeo[:] = netCDF4.date2num(df2.index.tolist(), units='days since 1980-01-01 00:00:00',
                            calendar='standard')

# Write the CBH values
varo[:, :] = df2.values

nco.close()

