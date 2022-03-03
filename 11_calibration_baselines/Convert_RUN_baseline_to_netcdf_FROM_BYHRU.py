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
import datetime
import numpy as np
import netCDF4 as nc
import pandas as pd
import sys

from pyPRMS.prms_helpers import dparse

# %%
# ===============
# NHM v1.1 format
# <no header>
# Year Month MWBMq min max

# nhru = 114958
# Space delimited

# filename pattern
# MWBM_*.csv

# ===============
# NHM v1.0 format
# <no header>
# Year Month MWBMq min max

# nhru = 109951
# Space delimited

# filename pattern
# HRU_*

# %%
baseline = 'RUN'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/baselines/conus'
workdir = '/Volumes/parker_rocks/calibrations/NHMv11/baselines'
file_prefix = 'MWBM_'
file_suffix = '.csv'
delimiter = ' '
date_fields = [0, 1]

st_date = datetime.datetime(1980,1,1)
en_date = datetime.datetime(2010,12,31)

nhru = 114958

hrus = np.arange(1, nhru+1)

# %%
the_vars = {'runoff_mwbm': {'standard_name': 'runoff_mwbm',
                            'description': 'Surface runoff from Monthly Water Balance Model (MWBM)',
                            'datatype': 'f4',
                            'units': 'ft3 s-1',
                            'fill_value': nc.default_fillvals['f4']},
            'runoff_max': {'standard_name': 'runoff_max',
                           'description': 'Monthly maximum surface runoff',
                           'datatype': 'f4', 
                           'units': 'ft3 s-1',
                           'fill_value': nc.default_fillvals['f4']},            
            'runoff_min': {'standard_name': 'runoff_min',
                           'description': 'Monthly minimum surface runoff',
                           'datatype': 'f4', 
                           'units': 'ft3 s-1',
                           'fill_value': nc.default_fillvals['f4']}}

# %%
first = True

for chru in hrus:
    sys.stdout.write(f'\rProcessing {chru:06d}')
    sys.stdout.flush()
    
    filename = f'{workdir}/{baseline}/{file_prefix}{chru}{file_suffix}'
    
    # NOTE: RUN baseline has no column header
    df = pd.read_csv(filename, sep=delimiter, skipinitialspace=True,
                     header=None,
                     date_parser=dparse, parse_dates={'date': date_fields},
                     engine='c', memory_map=True, index_col='date',
                     na_values=[-999.0])
    df.rename(columns={2: 'mwbmq', 3: 'min', 4: 'max'}, inplace=True)
    
    # NOTE: 2021-04-26 - there is a change in date range at 26709 (1949-) so we have to specify 
    #       the date range we want (1980-2010)
    df = df[st_date:en_date]
    
    if first:
        date_vals = df.index.tolist()
        st_date = df.index[0]
        max_vals = np.zeros((len(df.index), len(list(hrus))))
        min_vals = np.zeros((len(df.index), len(list(hrus))))
        mwbm_vals = np.zeros((len(df.index), len(list(hrus))))
        
        # Start with NaNs so it's easier to figure out if HRUs are missing
        max_vals.fill(np.nan)
        min_vals.fill(np.nan)
        mwbm_vals.fill(np.nan)
        first = False
        
    max_vals[:, chru-1] = df['max'].to_numpy()
    min_vals[:, chru-1] = df['min'].to_numpy()
    mwbm_vals[:, chru-1] = df['mwbmq'].to_numpy()

# %%
print(f'max_vals: min={np.nanmin(max_vals)}, max={np.nanmax(max_vals)}')
print(f'min_vals: min={np.nanmin(min_vals)}, max={np.nanmax(min_vals)}')
print(f'mwbm_vals: min={np.nanmin(mwbm_vals)}, max={np.nanmax(mwbm_vals)}')

# %%
max_vals[np.isnan(max_vals)] = nc.default_fillvals['f4']
min_vals[np.isnan(min_vals)] = nc.default_fillvals['f4']
mwbm_vals[np.isnan(mwbm_vals)] = nc.default_fillvals['f4']

# %%

# %%

# %%
outfile = f'{workdir}/{baseline}/baseline_{baseline}.nc'

# Create a netCDF file for the CBH data
nco = nc.Dataset(outfile, 'w', clobber=True)
nco.createDimension('nhru', nhru)
nco.createDimension('time', None)

timeo = nco.createVariable('time', 'f4', ('time'))
timeo.calendar = 'standard'
timeo.units = f'days since {st_date.year}-{st_date.month:02}-01 00:00:00'

nhruo = nco.createVariable('nhru', 'i4', ('nhru'))
nhruo.long_name = 'Local HRU id'
nhruo.description = 'Local Hydrologic Response Unit ID (HRU)'

nhmido = nco.createVariable('nhm_id', 'i4', ('nhru'))
nhmido.long_name = 'NHM HRU id'
nhmido.description = 'National Hydrologic Response Unit ID (HRU)'

for var, info in the_vars.items():
    varo = nco.createVariable(var, info['datatype'], ('time', 'nhru'), fill_value=info['fill_value'], zlib=True)
    varo.standard_name = info['standard_name']
    varo.description = info['description']
    varo.units = info['units']

nco.setncattr('Description', 'Baseline surface runoff')

# Write the HRU ids
nhruo[:] = hrus
nhmido[:] = hrus

timeo[:] = nc.date2num(date_vals, units=timeo.units, calendar=timeo.calendar)

nco.variables['runoff_max'][:, :] = max_vals
nco.variables['runoff_min'][:, :] = min_vals
nco.variables['runoff_mwbm'][:, :] = mwbm_vals

nco.close()

# %%
# nco.close()
outfile

# %%
# Columns (82897,82942,82976) have mixed types (RUNmax)

for ii in range(nhru):
    ttype = df.iloc[:, ii].dtype
    
    if ttype != 'float64':
        print(f'{ii}: {ttype}')

# %%
df.iloc[:, 82974].values

# %%
np.argwhere(np.isnan(max_vals))
# np.argwhere(np.isnan(df.to_numpy()))

# %%
df.iloc[188, 82895]

# %%
df.iloc[:,1]

# %%
aa = np.argwhere(max_vals < 0.0)

# %%
aa.shape

# %%
bb = list(set(aa[:, 1]))

# %%
bb.sort()

# %%
bb

# %%
