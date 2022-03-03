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
import os
import xarray as xr
import pandas as pd

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v3_1980-01-01_2016-12-31/by_year'

single_file = '{}/daymet_v3_cbh_all_19800101-19801231.nc'.format(workdir)
multifiles = '{}/daymet_v3_cbh_all_*.nc'.format(workdir)

# %%
# ds = xr.open_dataset(single_file)
ds = xr.open_mfdataset(multifiles, chunks={'hru': 1000})


# ds = xr.open_mfdataset("/glade/p/cesm/community/ASD-HIGH-RES-CESM1/hybrid_v5_rel04_BC5_ne120_t12_pop62/ocn/proc/tseries/monthly/*.nc", 
#                        parallel=True, coords="minimal", data_vars="minimal", compat='override')

# %%
ds.keys()

# %%
for xx in ds.data_vars:
    print(xx)

# %%
# ds['hru'].loc[nhm_hrus]
pd.to_datetime(ds['time'].loc[st_date:en_date].values).tolist()

# %%
st_date = datetime.datetime(1980, 1, 1)
en_date = datetime.datetime(1981, 12, 31)
nhm_hrus = [1,5,2,3]

cvar = 'tmax'
da = ds[cvar].loc[:, nhm_hrus]

# %%
da

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
