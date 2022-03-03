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
# nhru = 114958

# location on yeti
# /cxfs/projects/usgs/water/mows/NHM/calibration/gf_v11/byHRU/baselines/AET

# filename pattern
# AET_*.csv

# comma-delimited

# header
# Year, Month, MWBM(in), MOD16(in), SSEBop(in)


# ===============
# NHM v1.0 format
# nhru = 109951

# location on yeti
# /cxfs/projects/usgs/water/mows/NHM/NHM_CONUS/NHM_CONUS_CALIBRATION/OLD/RUN_CALbyHRU/AET/byHRU

# filename pattern
# HRU_*

# space-delimited

# header
# Year Month MWBM(in) MOD16(in) SSEBop(in)

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/Maurer/baselines/byHRU/RUN/byHRU'
baseline = 'AET'
workdir = '/Volumes/parker_rocks/calibrations/NHMv11/baselines'

outfile = f'{workdir}/{baseline}/baseline_{baseline}.nc'

file_prefix = 'AET_'
file_suffix = '.csv'
delimiter = ','
date_fields = [0, 1]
file_nan_values = [-999.0]

global_description = 'Baseline actual evapotranspiration'

nhru = 114958

hrus = np.arange(1, nhru+1)

# %%
the_vars = {'aet_mwbm': {'standard_name': 'aet_mwbm',
                         'description': 'Actual evapotranspiration from Monthly Water Balance Model (MWBM)',
                         'datatype': 'f4',
                         'units': 'inch',
                         'fill_value': nc.default_fillvals['f4']},
            'aet_mod16': {'standard_name': 'aet_mod16',
                          'description': 'Actual evapotranspiration from MODIS-16 global evapotranspiration product',
                          'datatype': 'f4',
                          'units': 'inch',
                          'fill_value': nc.default_fillvals['f4']},
            'aet_ssebop': {'standard_name': 'aet_ssebop',
                           'description': 'Actual evapotranspiration from Simplified Surface Energy Balance (SSEBop)',
                           'datatype': 'f4',
                           'units': 'inch',
                           'fill_value': nc.default_fillvals['f4']},
            'aet_max': {'standard_name': 'aet_max',
                        'description': 'Monthly maximum actual evapotranspiration from MWBM, MOD16, and SSEBop',
                        'datatype': 'f4', 
                        'units': 'inch',
                        'fill_value': nc.default_fillvals['f4']},            
            'aet_min': {'standard_name': 'aet_min',
                        'description': 'Monthly minimum actual evapotranspiration from MWBM, MOD16, and SSEBop',
                        'datatype': 'f4', 
                        'units': 'inch',
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
                     na_values=[-999.0])
    
    if first:
        date_vals = df.index.tolist()
        st_date = df.index[0]
        ntime = len(df.index)
        
        # Start with NaNs so it's easier to figure out if HRUs are missing
        max_vals = np.zeros((ntime, nhru)).fill(np.nan)
        min_vals = np.zeros((ntime, nhru)).fill(np.nan)
        mwbm_vals = np.zeros((ntime, nhru)).fill(np.nan)
        mod16_vals = np.zeros((ntime, nhru)).fill(np.nan)
        ssebop_vals = np.zeros((ntime, nhru)).fill(np.nan)
        
        first = False
        
    mwbm_vals[:, chru-1] = df['MWBM(in)'].to_numpy()
    mod16_vals[:, chru-1] = df['MOD16(in)'].to_numpy()
    ssebop_vals[:, chru-1] = df['SSEBop(in)'].to_numpy()
    max_vals[:, chru-1] = np.fmax.reduce([mwbm_vals[:, chru-1], mod16_vals[:, chru-1], ssebop_vals[:, chru-1]])
    min_vals[:, chru-1] = np.fmin.reduce([mwbm_vals[:, chru-1], mod16_vals[:, chru-1], ssebop_vals[:, chru-1]])

# %%
print(f'max_vals: min={np.nanmin(max_vals)}, max={np.nanmax(max_vals)}')
print(f'min_vals: min={np.nanmin(min_vals)}, max={np.nanmax(min_vals)}')
print(f'mwbm_vals: min={np.nanmin(mwbm_vals)}, max={np.nanmax(mwbm_vals)}')
print(f'mod16_vals: min={np.nanmin(mod16_vals)}, max={np.nanmax(mod16_vals)}')
print(f'ssebop_vals: min={np.nanmin(ssebop_vals)}, max={np.nanmax(ssebop_vals)}')

# %%
# Convert the NaN values to the fill value for each variable
max_vals[np.isnan(max_vals)] = the_vars['aet_max']['fill_value']
min_vals[np.isnan(min_vals)] = the_vars['aet_min']['fill_value']
mwbm_vals[np.isnan(mwbm_vals)] = the_vars['aet_mwbm']['fill_value']
mod16_vals[np.isnan(mod16_vals)] = the_vars['aet_mod16']['fill_value']
ssebop_vals[np.isnan(ssebop_vals)] = the_vars['aet_ssebop']['fill_value']

# %%
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

nco.setncattr('Description', global_description)

# Write the HRU ids
nhruo[:] = hrus
nhmido[:] = hrus

timeo[:] = nc.date2num(date_vals, units=timeo.units, calendar=timeo.calendar)

nco.variables['aet_max'][:, :] = max_vals
nco.variables['aet_min'][:, :] = min_vals
nco.variables['aet_mwbm'][:, :] = mwbm_vals
nco.variables['aet_mod16'][:, :] = mod16_vals
nco.variables['aet_ssebop'][:, :] = ssebop_vals

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
adate = [2001, 10, 1]
pd.to_datetime(adate)

# %%
st_date = datetime.datetime(1980,1,1)
pd.to_datetime(st_date)

# %%
df = pd.read_csv('/Volumes/parker_rocks/calibrations/NHMv10/AET/byHRU/HRU_1', sep=' ', skipinitialspace=True,
                 parse_dates={'date': date_fields},
#                  date_parser=dparse, parse_dates={'date': date_fields},
                 engine='c', memory_map=True, index_col='date',
                 na_values=[-999.0])

# %%
df.info()

# %%
df = pd.read_csv('/Volumes/parker_rocks/calibrations/NHMv10/RCH/byHRUranges/HRU_1', sep=' ', skipinitialspace=True,
                 parse_dates={'date': [0]},
                 engine='c', memory_map=True, index_col='date',
                 na_values=[-999.0])

# %%
df.info()

# %%
df['MIN']

# %%
