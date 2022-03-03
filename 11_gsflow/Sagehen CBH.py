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
# %matplotlib inline
import matplotlib.pyplot as plt

import datetime
import netCDF4 as nc
import numpy as np
import pandas as pd
from pyPRMS.CbhAscii import CbhAscii
from pyPRMS.prms_helpers import dparse

# %%
pd.__version__

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/sagehen/prms6'
st_date = datetime.datetime(1980, 10, 1)
en_date = datetime.datetime(1996, 9, 30)

CBH_INDEX_COLS = [0, 1, 2, 3, 4, 5]

hrus = list(range(1,129))


# %%

# %%

# %%
def read_cbh(filename, hruids):
    incl_cols = list(CBH_INDEX_COLS)

    for xx in hrus:
        incl_cols.append(xx+5)  # include an offset for the having datetime info
    # print(incl_cols)

    # Columns 0-5 always represent date/time information
    df = pd.read_csv(filename, sep=' ', skipinitialspace=True, usecols=incl_cols,
                     skiprows=3, engine='c', memory_map=True,
                     date_parser=dparse, parse_dates={'time': CBH_INDEX_COLS},
                     index_col='time', header=None, na_values=[-99.0, -999.0])
    
#     if self.__stdate is not None and self.__endate is not None:
#         self.__data = self.__data[self.__stdate:self.__endate]

    # self.__data.reset_index(drop=True, inplace=True)

    # Rename columns with NHM HRU ids
    ren_dict = {v + 5: v for v in hruids}

    # NOTE: The rename is an expensive operation
    df.rename(columns=ren_dict, inplace=True)
    
    return df


# %%
def get_cbh_xarray(workdir, varname, hrus):
    df = read_cbh(f'{workdir}/{varname}.day', hrus)
    # df.head()

    # Unpivot into correct form for transforming to xarray
    df2 = pd.melt(df, ignore_index=False, var_name='hru_id', value_name=varname, value_vars=list(range(1,129)))
    df2.reset_index(inplace=True)
    df2.set_index(['time', 'hru_id'], inplace=True)
    # df2.head()

    # Convert to an xarray
    xdf = df2.to_xarray()
    # xdf.info()
    return xdf


# %%
# varname = 'precip'

tmax = get_cbh_xarray(workdir, 'tmax', hrus)['tmax']
tmin = get_cbh_xarray(workdir, 'tmin', hrus)['tmin']
precip = get_cbh_xarray(workdir, 'precip', hrus)['precip']

# %%

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(17,11))
ax = axes.flatten()

precip.plot(ax=ax[0], cmap='tab20c_r')
tmin.plot(ax=ax[1], cmap='bwr')
tmax.plot(ax=ax[2], cmap='bwr')

# %%
precip[:, 70].plot()

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20, 6))
ax = axes.flatten()

hru0 = 70
precip[:, hru0].plot(ax=ax[0])
tmin[:, hru0].plot(ax=ax[1])
tmax[:, hru0].plot(ax=ax[2])

# %%
