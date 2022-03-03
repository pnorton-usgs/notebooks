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
import xarray as xr

from Bandit.model_output import ModelOutput

# %%
workdir = '/Users/pnorton/tmp'
filename = f'{workdir}/seg_gwflow.nc'

filename_rain = f'{workdir}/hru_rain.nc'

st_date = datetime.datetime(1980, 10, 1)
en_date = datetime.datetime(1980, 10, 31)

# the_segs = [49153, 49154, 49155, 49156, 49157]
the_segs = [49153]

# %%
modout = ModelOutput(filename=filename, varname='seg_gwflow', startdate=st_date, enddate=en_date, nhm_segs=the_segs)

# %%
modout.get_var('seg_gwflow').head()

# %%
modout.dataset['seg_gwflow'].to_pandas().head()

# %%
# aa = xr.open_dataset(filename, decode_coords=True, chunks={'nsegment': 1000})
aa = xr.open_dataset(filename_rain, decode_coords=True, chunks={'nhru': 1000})

# %%
aa

# %%
# da.sel(a=da.c.to_index().get_indexer(['x', 'y']))
# aa[seg_gwflow].sel(nsegment=)

a = aa['seg_gwflow']
# a.sel(country=a.currency == 'EUR')
a.sel(nsegment=a.seg_id == 49154).to_pandas()

# %%
# In [63]: da = xr.DataArray(np.random.rand(3,2), dims=list('ab'), coords={'c':(('a',),list('xyz'))})
# In [64]: da.sel(a=(np.isin(da.c, list('xy'))))
a.sel(nsegment=(np.isin(a.seg_id, [49157, 49154]))).to_pandas()

# %%

# %%

# %%
a

# %%
b = aa['seg_id']
# b.sel(nsegment=b.seg_id in [49154, 49157]).to_pandas()
# data = self.__dataset[varname].loc[:, self.__nhm_hrus].to_pandas()
b.loc[[49154, 49157]].to_pandas()

# %%
# Get the indices for the NHM (global) ids
cc = b.to_index().get_indexer([49157, 49154])

# %%
dd = a.loc[:, cc].to_pandas()

# %%
dd

# %%
ee = a.loc[:, cc]

# %%
ee

# %%
# da.assign_coords(lon=(((da.lon + 180) % 360) - 180))
# aa = xr.open_dataset(filename, decode_coords=True, chunks={'nsegment': 1000})

# aa = xr.open_dataset(filename, chunks={'nsegment': 1000})
aa = xr.open_dataset(filename_rain, decode_coords=True, chunks={'nhru': 1000})
aa

# %%
ba = aa.assign_coords(nhru=(aa.nhm_id))

# %%
dd = ba['hru_rain'].loc[:, [101, 102]].to_pandas()

# %%
dd

# %%
ba

# %%
list(range(1,11))

# %%
