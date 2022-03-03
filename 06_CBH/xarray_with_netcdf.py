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
import os
import xarray as xr
import pandas as pd

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/gridmet/gm_gf_v11_US_with_derived'
workdir = '/Volumes/USGS_NHM2/datasets/gridmet/gm_gf_v11_US_UNITS_filled'
# single_file = '{}/daymet_v3_cbh_all_19800101-19801231.nc'.format(workdir)
# multifiles = '{}/daymet_v3_cbh_all_*.nc'.format(workdir)

# %%
# ds = xr.open_dataset(single_file)
# add either: 
# combine='nested' (original behavior before v0.15)
# combine='by_coords' (future default behavior)
ds = xr.open_mfdataset(f'{workdir}/*_1979.nc', chunks={'hruid': 1040}, combine='by_coords', decode_cf=True,
                       attrs_file=f'{workdir}/gm_climate_1979.nc')

# %%
ds

# %%
ds['tmin'].encoding['_FillValue']

# %%

# %%

# %%
for xx in ds.data_vars:
    print(xx)

# %%
for xx in ds['tmax'].coords:
    print(xx)

# %%
ds['time'].attrs

# %%
ds.encoding

# %%
# ds['hru'].loc[nhm_hrus]
st_date = datetime.datetime(1980,1,1)
en_date = datetime.datetime(1980,6,6)
pd.to_datetime(ds['time'].loc[st_date:en_date].values).tolist()

# %%
# nhm_hrus = [1,5,2,3]
nhm_hrus = [57873, 57874, 57877, 57880, 57867, 57872, 57878, 57879, 57881, 57882, 57868, 57869, 57863, 57864]
cvar = 'tmax'
da = ds[cvar].loc[:, nhm_hrus]

# %%
da

# %%
# Select a date range
# self.__dataset[var].loc[self.__stdate:self.__endate, self.__nhm_hrus].to_pandas()
da = ds[cvar].loc[st_date:en_date].to_pandas()
da.head()

# %%
# Select date range and HRUs as the same time
da = ds[cvar].loc[st_date:en_date, nhm_hrus].to_pandas()
da.head()

# %%
print(type(st_date))

# %%
from pyPRMS.Cbh import Cbh

# %%
st_date = datetime.datetime(1981, 1, 1)
en_date = datetime.datetime(1984, 12, 31)
nhm_hrus = [1,5,2,3]

cbh_obj = Cbh(filename=multifiles, st_date=st_date, en_date=en_date, nhm_hrus=nhm_hrus)

# %%
cbh_obj.write_netcdf('crap_all.nc')

# %%
import os

os.path.exists(multifiles)

# %%
cbh_obj.write_ascii('crap')

# %%
cbh_obj.write_netcdf('crap.nc')

# %%
os.path.basename(multifiles)

# %%
aa = [1,2,3]
isinstance(aa, list)

# %%
type(aa)

# %%
da

# %%
