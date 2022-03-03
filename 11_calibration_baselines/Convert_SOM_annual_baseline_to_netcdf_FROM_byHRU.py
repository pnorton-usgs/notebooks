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

# comma delimited

# filename pattern
# SOMann_*.csv


# ===============
# NHM v1.0 format
# CAL_YR   MIN     MAX    MERRA    NCEP-NCAR   MOSAIC     NOAH

# space delimited

# filename pattern
# annHRU_*

# %%
baseline = 'SOMann'
workdir = '/Volumes/parker_rocks/calibrations/NHMv11/baselines'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/baselines/byHRU'
file_prefix = 'SOMann_'
file_suffix = '.csv'
delimiter = ','
date_fields = [0]

nhru = 114958

hrus = np.arange(1, nhru+1)

# %%

# %%
first = True

for chru in hrus:
#     sys.stdout.write('\r                                       ')
    sys.stdout.write(f'\rProcessing {chru:06d}')
    sys.stdout.flush()
    
    filename = f'{workdir}/{baseline}/{file_prefix}{chru}{file_suffix}'
    
    df = pd.read_csv(filename, sep=delimiter, skipinitialspace=True,
                     date_parser=dparse, parse_dates={'date': date_fields},
                     engine='c', memory_map=True, index_col='date')
    
    if first:
        date_vals = df.index.tolist()
        st_date = df.index[0]
        max_vals = np.zeros((len(df.index), len(list(hrus))))
        min_vals = np.zeros((len(df.index), len(list(hrus))))
        
        # Start with NaNs so it's easier to figure out if HRUs are missing
        max_vals.fill(np.nan)
        min_vals.fill(np.nan)
        first = False
        
    max_vals[:, chru-1] = df['max'].to_numpy()
    min_vals[:, chru-1] = df['min'].to_numpy()

# %%
print(f'max_vals: min={np.nanmin(max_vals)}, max={np.nanmax(max_vals)}')
print(f'min_vals: min={np.nanmin(min_vals)}, max={np.nanmax(min_vals)}')

# %%
#  MIN     MAX    MERRA    NCEP-NCAR   MOSAIC     NOAH
# MIN: soil_moist_min_norm
# MAX: soil_moist_max_norm
# MERRA: soil_moist_merra_norm
# NCEP-NCAR: soil_moist_ncep_norm
# MOSAIC: soil_moist_mosaic_norm
# NOAH: soil_moist_noah_norm

the_vars = {'soil_moist_min_norm': {'standard_name': 'soil_moist_min_norm',
                                    'description': 'Normalized annual minimum top-layer soil moisture',
                                    'datatype': 'f4',
                                    'units': 'fraction',
                                    'fill_value': nc.default_fillvals['f4']},
            'soil_moist_max_norm': {'standard_name': 'soil_moist_max_norm',
                                    'description': 'Normalized annual maximum top-layer soil moisture',
                                    'datatype': 'f4',
                                    'units': 'fraction',
                                    'fill_value': nc.default_fillvals['f4']}}
#             'soil_moist_merra_norm': {'standard_name': 'soil_moist_merra_norm',
#                                       'description': 'Normalized monthly top-layer soil moisture derived from MERRA-Land reanalysis',
#                                       'datatype': 'f4',
#                                       'units': 'fraction',
#                                       'fill_value': nc.default_fillvals['f4']},
#             'soil_moist_ncep_norm': {'standard_name': 'soil_moist_ncep_norm',
#                                      'description': 'Normalized monthly top-layer soil moisture derived from NCEP-NCAR reanalysis',
#                                      'datatype': 'f4',
#                                      'units': 'fraction',
#                                      'fill_value': nc.default_fillvals['f4']},
#             'soil_moist_mosaic_norm': {'standard_name': 'soil_moist_mosaic_norm',
#                                        'description': 'Normalized monthly top-layer soil moisture derived from NLDAS-MOSAIC model',
#                                        'datatype': 'f4',
#                                        'units': 'fraction',
#                                        'fill_value': nc.default_fillvals['f4']},
#             'soil_moist_noah_norm': {'standard_name': 'soil_moist_noah_norm',
#                                      'description': 'Normalized monthly top-layer soil moisture derived from NLDAS-NOAH model',
#                                      'datatype': 'f4',
#                                      'units': 'fraction',
#                                      'fill_value': nc.default_fillvals['f4']}}


outfile = f'{workdir}/{baseline}/baseline_SOMann.nc'

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

nco.setncattr('Description', 'Baseline soil moisture information')

# Write the HRU ids
nhruo[:] = hrus
nhmido[:] = hrus

timeo[:] = nc.date2num(date_vals, units=timeo.units, calendar=timeo.calendar)

nco.variables['soil_moist_min_norm'][:, :] = min_vals
nco.variables['soil_moist_max_norm'][:, :] = max_vals

nco.close()


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

# %%
np.argwhere(np.isnan(min_vals))

# %%
