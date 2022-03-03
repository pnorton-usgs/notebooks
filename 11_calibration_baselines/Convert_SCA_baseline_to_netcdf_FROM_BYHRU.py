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
# "year","month","day","sca","clearidx"

# filename pattern
# SCA_*.csv

# Comma-delimited


# ===============
# NHM v1.0 format
# <NO HEADER>
# year month day sca ci

# filename pattern
# HRU_*

# space delimited

# %%
baseline = 'SCA'
workdir = '/Volumes/parker_rocks/calibrations/NHMv11/baselines'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/baselines/conus'
file_prefix = 'SCA_'
file_suffix = '.csv'
delimiter = ','
date_fields = [0, 1, 2]

nhru = 114958

hrus = np.arange(1, nhru+1)

# %%
the_vars = {'snow_cover_extent': {'standard_name': 'snow_cover_extent',
                                  'description': 'Daily snow extent',
                                  'datatype': 'f4',
                                  'units': '1',
                                  'fill_value': nc.default_fillvals['f4']},
            'sca_clear_index': {'standard_name': 'sca_clear_index',
                                'description': 'Clear index for the daily snow map',
                                'datatype': 'f4',
                                'units': 'percent',
                                'fill_value': nc.default_fillvals['f4']}}

# %%
first = True

for chru in hrus:
    sys.stdout.write(f'\rProcessing {chru:06d}')
    sys.stdout.flush()
    
    filename = f'{workdir}/{baseline}/{file_prefix}{chru}{file_suffix}'
    
    df = pd.read_csv(filename, sep=delimiter, skipinitialspace=True,
                     date_parser=dparse, parse_dates={'date': date_fields},
                     engine='c', memory_map=True, index_col='date',
                     na_values=[-9.0, -9.99, -99.0, -99.99, -999.0, -9999.0])
    
    if first:
        date_vals = df.index.tolist()
        st_date = df.index[0]
        sca_vals = np.zeros((len(df.index), len(list(hrus))))
        scaci_vals = np.zeros((len(df.index), len(list(hrus))))
        
        # Start with NaNs so it's easier to figure out if HRUs are missing
        sca_vals.fill(np.nan)
        scaci_vals.fill(np.nan)
        first = False
        
    sca_vals[:, chru-1] = df['sca'].to_numpy()
    scaci_vals[:, chru-1] = df['clearidx'].to_numpy()

# %%
print(f'sca_vals: min={np.nanmin(sca_vals)}, max={np.nanmax(sca_vals)}')
print(f'scaci_vals: min={np.nanmin(scaci_vals)}, max={np.nanmax(scaci_vals)}')

# %%
sca_vals[np.isnan(sca_vals)] = nc.default_fillvals['f4']
scaci_vals[np.isnan(scaci_vals)] = nc.default_fillvals['f4']

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

nco.setncattr('Description', 'Baseline snow covered extent')

# Write the HRU ids
nhruo[:] = hrus
nhmido[:] = hrus

timeo[:] = nc.date2num(date_vals, units=timeo.units, calendar=timeo.calendar)

nco.variables['snow_cover_extent'][:, :] = sca_vals
nco.variables['sca_clear_index'][:, :] = scaci_vals

nco.close()

# %%

# %%

# %%
df.info()

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
df.max()

# %%
