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

from pyPRMS.prms_helpers import dparse

# %%
# SOM files
# baselines.SOMmthMAX
# baselines.SOMmthMEAN
# baselines.SOMmthMERRA
# baselines.SOMmthMIN
# baselines.SOMmthMOSAIC
# baselines.SOMmthMdR
# baselines.SOMmthNCEP
# baselines.SOMmthNOAH
# baselines.SOMmthRNG

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/Maurer/baselines/byHRU/RUN/byHRU'
cvar = 'SOMmthNOAH'
baseline = 'SOM'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/baselines/conus'

nhru = 109951

hrus = np.arange(1, nhru+1)

# %%
var_info = {'SOMmthMAX': ['Maximum soil moisture', 'inch'], 
            'SOMmthMIN': ['Minimum soil moisture', 'inch'], 
            'SOMmthMEAN': ['Mean soil moisture', 'inch'],
            'SOMmthMdR': ['mdr soil moisture', 'inch'],
            'SOMmthRNG' : ['Range of soil moisture', 'inch'],
            'SOMmthMERRA': ['Soil moisture from MERRA', 'inch'],
            'SOMmthMOSAIC': ['Soil moisture from MOSAIC', 'inch'],
            'SOMmthNCEP': ['Soil moisture from NCEP', 'inch'],
            'SOMmthNOAH': ['Soil moisture from NOAH', 'inch']}

delim = {'SOMmthMAX': ' ', 'SOMmthMIN': ' ', 'SOMmthMEAN': ',', 'SOMmthMdR': ',', 'SOMmthRNG': ',',
         'SOMmthMERRA': ' ', 'SOMmthMOSAIC': ' ', 'SOMmthNCEP': ' ', 'SOMmthNOAH': ' '}

filename = f'{workdir}/{baseline}/baselines.{cvar}'

df = pd.read_csv(filename, sep=delim[cvar], skipinitialspace=True,
                 date_parser=dparse, parse_dates={'date': [0, 1]},
                 engine='c', memory_map=True, index_col='date')

# %%
df.info()

# %%
st_date = df.index[0]

outfile = f'{workdir}/{baseline}/{cvar}.nc'

# Create a netCDF file for the CBH data
nco = nc.Dataset(outfile, 'w', clobber=True)
nco.createDimension('hru', nhru)
nco.createDimension('time', None)

timeo = nco.createVariable('time', 'f4', ('time'))
timeo.calendar = 'standard'
# timeo.bounds = 'time_bnds'

# FIXME: Days since needs to be set to the starting date of the model pull
timeo.units = f'days since {st_date.year}-{st_date.month:02}-01 00:00:00'

hruo = nco.createVariable('hru', 'i4', ('hru'))
hruo.long_name = 'Hydrologic Response Unit ID (HRU)'

varo = nco.createVariable(cvar, 'f4', ('time', 'hru'), fill_value=nc.default_fillvals['f4'], zlib=True)
varo.long_name = var_info[cvar][0]
varo.units = var_info[cvar][1]

nco.setncattr('Description', 'Baseline soil moisture information')

# Write the HRU ids
hruo[:] = hrus

timeo[:] = nc.date2num(df.index.tolist(), units=timeo.units, calendar='standard')
nco.variables[cvar][:, :] = df.values

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
np.argwhere(np.isnan(df.to_numpy()))

# %%
