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

# %% [markdown]
# This notebook is from https://nbviewer.jupyter.org/gist/rsignell-usgs/d3dfaf3cd3d8b39894a69b22127dfe38

# %%
import xarray as xr
import fsspec
import numpy as np

import hvplot.pandas
import hvplot.xarray
import geoviews as gv
from holoviews.operation.datashader import rasterize
import cartopy.crs as ccrs

# Also requires: s3fs, zarr

# %%
from dask.distributed import Client
client = Client()
client

# %%
url = 's3://noaa-nwm-retro-v2-zarr-pds'

# %%
# %%time
ds = xr.open_zarr(fsspec.get_mapper(url, anon=True), consolidated=True)

# %%
var = 'streamflow'
ds[var]

# %%
print(f'Variable size: {ds[var].nbytes/1e12:.1f} TB')

# %% [markdown]
# ## Find the site with the largest streamflow on 2017-06-01

# %%
# %%time
imax = ds[var].sel(time='2017-06-01 00:00:00').argmax().values

# %% [markdown]
# ### Plot the hindcast time series at that location

# %%
# %%time
ds[var][:,imax].hvplot(grid=True)

# %% [markdown]
# ## Compute mean discharge during April 2010 on all rivers

# %%
streamflow_April_2010 = ds[var].sel(time=slice('2010-04-01 00:00','2010-04-30 23:00'))

# %%
print(f'Variable size: {streamflow_April_2010.nbytes/1e9:.1f} GB')

# %%
# %%time
var_mean = streamflow_April_2010.mean(dim='time').compute()

# %% [markdown]
# ### Visualize the mean discharge with hvplot

# %%
df = var_mean.to_pandas().to_frame()

# %% [markdown]
# The dataframe just has streamflow, so add longitude and latitude as columns

# %%
df = df.assign(latitude=ds['latitude'])
df = df.assign(longitude=ds['longitude'])
df.rename(columns={0: "transport"}, inplace=True)

# %%
p = df.hvplot.points('longitude', 'latitude', crs=ccrs.PlateCarree(),
                     c='transport', colorbar=True, size=14)

# %% [markdown]
# We don't want to plot all the 2.7M points individually, so aggregate to 0.02 degree resolution and rasterize with datashader. Use a log scale for visualization since there is a large dynamic range in streamflow.

# %%
g = rasterize(p, aggregator='mean', x_sampling=0.02, y_sampling=0.02, width=500).opts(tools=['hover'], 
                aspect='equal', logz=True, cmap='viridis', clim=(1e-2, np.nan))

# %% [markdown]
# Plot the rasterized streamflow data on an OpenStreetMap tile service basemap

# %%
g * gv.tile_sources.OSM

# %%
#client.close(); cluster.close()
