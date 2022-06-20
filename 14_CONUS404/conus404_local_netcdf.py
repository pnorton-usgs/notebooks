# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python [conda env:holoview]
#     language: python
#     name: conda-env-holoview-py
# ---

# %%
import fsspec
import xarray as xr
import hvplot.xarray
import intake
import os
import numpy as np
import warnings
from matplotlib import path
import panel as pn  

# %%
warnings.filterwarnings('ignore')

# %%
from dask.distributed import Client, progress
from dask import delayed

# Dask cluster creation option from notebook
from dask_kubernetes import KubeCluster
cluster = KubeCluster()
cluster.scale(20)
client = Client(cluster)
# client

# %%
cluster

# %%
client.close()

# %%
fs = fsspec.filesystem('s3', requester_pays=True, profile='nhgf-s3')


# %%
def nearxy(x,y,xi,yi):
    ind = np.ones(len(xi),dtype=int)
    for i in range(len(xi)):
        dist = np.sqrt((x-xi[i])**2+(y-yi[i])**2)
        ind[i] = dist.argmin()
    return ind

def ind2sub(array_shape, ind):
    rows = int(ind.astype('int') / array_shape[1])
    cols = int(ind.astype('int') % array_shape[1]) # or numpy.mod(ind.astype('int'), array_shape[1])
    return (rows, cols)


# %%

# %%
# %%time
ds = xr.open_dataset(fs.open('s3://nhgf-development/ashalper/conus404_drb/DRB_WY2017.nc'))
#                      chunks={'time':2921, 'y':38, 'x':29})

# %%
ds

# %%

# %%
lat = 41.3588104
lon = -76.2230835

# Get the west_east, south_north indices for the lat/lon
[jj,ii] = ind2sub(ds.lon.shape, nearxy(ds.lon.values, ds.lat.values, [lon], [lat]))

# %%
[jj, ii]

# %% jupyter={"outputs_hidden": true}
# %%time
ss = ds.isel(x=ii, y=jj)
ss
# frozen_precip = (ds.SNOW_ACC_NC + ds.GRAUPEL_ACC_NC).sum(dim='time')

# %% jupyter={"source_hidden": true}
# %%time
frozen_precip = (ss.SNOW_ACC_NC + ss.GRAUPEL_ACC_NC).sum(dim='time')
frozen_precip.compute()

# %%
# %%time
# Why doesn't this work?
chg_snowmelt = ss['ACSNOM'].sel(time='2017-09-30 23:00:00') - ss['ACSNOM'].sel(time='2016-10-01 01:00:00')

# chg_snowmelt = ss['ACSNOM'][-1] - ss['ACSNOM'][0]

# %%
chg_snowmelt

# %%
ss['ACSNOM'].sel(time='2017-09-30 23:00:00')

# %%
ss['ACSNOM'][-1]

# %%
ss['SNOW'].plot()

# %%
ss['SNOW'].sum(dim='time')

# %%

# %%

# %%
