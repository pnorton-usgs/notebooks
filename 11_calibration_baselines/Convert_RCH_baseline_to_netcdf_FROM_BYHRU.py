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
import numpy as np
import netCDF4 as nc
import pandas as pd
import sys

from pyPRMS.prms_helpers import dparse

# %%
# ===============
# NHM v1.1 format
# "Year","min","max"

# nhru=114958

# Comma delimited

# filename pattern
# RCH_*.csv

# Values are normalized 0 to 1


# ===============
# NHM v1.0 format
# CAL_YR  MIN  MAX   REITZ WATERGAP

# space-delimited

# filename pattern
# HRU_*

# nhru = 109951

# %%
baseline = 'RCH'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/Maurer/baselines/byHRU/RUN/byHRU'
workdir = '/Volumes/parker_rocks/calibrations/NHMv11/baselines'
file_prefix = 'RCH_'
file_suffix = '.csv'
delimiter = ','
date_fields = [0]

nhru = 114958

hrus = np.arange(1, nhru+1)

# %%
the_vars = {'recharge_min_norm': {'standard_name': 'recharge_min_norm',
                                    'description': 'Normalized annual minimum recharge',
                                    'datatype': 'f4',
                                    'units': '1',
                                    'fill_value': nc.default_fillvals['f4']},
            'recharge_max_norm': {'standard_name': 'recharge_max_norm',
                                    'description': 'Normalized annual maximum recharge',
                                    'datatype': 'f4',
                                    'units': '1',
                                    'fill_value': nc.default_fillvals['f4']}}

# var_info = {'RCHmax': ['Maximum recharge', 'inch'], 
#             'RCHmin': ['Minimum recharge', 'inch'], 
#             'RCHmean': ['Mean recharge', 'inch'],
#             'RCHrng' : ['Range of recharge', 'inch'],
#             'RCHreitz': ['Recharge from Reitz', 'inch'],
#             'RCHwatergap': ['Recharge from Watergap model', 'inch']}

# %%
first = True

for chru in hrus:
    sys.stdout.write(f'\rProcessing {chru:06d}')
    sys.stdout.flush()
    
    filename = f'{workdir}/{baseline}/{file_prefix}{chru}{file_suffix}'
    
    df = pd.read_csv(filename, sep=delimiter, skipinitialspace=True,
                     date_parser=dparse, parse_dates={'date': date_fields},
                     engine='c', memory_map=True, index_col='date',
                     na_values=[-999.0])
    
    if first:
        date_vals = df.index.tolist()
        st_date = df.index[0]
        max_vals = np.zeros((len(df.index), len(list(hrus))))
        min_vals = np.zeros((len(df.index), len(list(hrus))))
        
        # Start with NaNs so it's easier to figure out if HRUs are missing
        max_vals.fill(np.nan)
        min_vals.fill(np.nan)
        first = False
        
    min_vals[:, chru-1] = df['min'].to_numpy()
    max_vals[:, chru-1] = df['max'].to_numpy()

# %%
print(f'min_vals: min={np.nanmin(min_vals)}, max={np.nanmax(min_vals)}')
print(f'max_vals: min={np.nanmin(max_vals)}, max={np.nanmax(max_vals)}')

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

nco.setncattr('Description', 'Baseline recharge')

# Write the HRU ids
nhruo[:] = hrus
nhmido[:] = hrus

timeo[:] = nc.date2num(date_vals, units=timeo.units, calendar=timeo.calendar)

nco.variables['recharge_max_norm'][:, :] = max_vals
nco.variables['recharge_min_norm'][:, :] = min_vals

nco.close()

# %%

# %%
nco.close()

# %%
# Columns (82897,82942,82976) have mixed types (RUNmax)

for ii in range(nhru):
    ttype = df.iloc[:, ii].dtype
    
    if ttype != 'float64':
        print(f'{ii}: {ttype}')

# %%
df.iloc[:, 82974].values

# %%
