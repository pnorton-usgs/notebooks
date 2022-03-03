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
#     display_name: Python [conda env:dask_1]
#     language: python
#     name: conda-env-dask_1-py
# ---

# %% [markdown]
# # Explore the National Water Model Reanalysis
# Use [Xarray](http://xarray.pydata.org/en/stable/), [Dask](https://dask.org) and [hvPlot](https://hvplot.holoviz.org) from the [HoloViz](https://holoviz.org) tool suite to explore the National Water Modle Reanalysis Version 2.  We read from a cloud-optimized [Zarr](https://zarr.readthedocs.io/en/stable/) dataset that is part of the [AWS Open Data Program](https://aws.amazon.com/opendata/), and we use a Dask cluster to parallelize computation and reading of data chunks.  

# %%
import xarray as xr
import fsspec
import numpy as np

# %%
import hvplot.pandas
import hvplot.xarray
import geoviews as gv
from holoviews.operation.datashader import rasterize
import cartopy.crs as ccrs

from dask.distributed import Client

# %% [markdown]
# ### Start a Dask cluster
# This is not required, but speeds up computations.  
#
# Here we are using [Qhub](https://www.quansight.com/post/announcing-qhub) which allows us to pass a specified conda environment to a Dask cluster running on Kubernetes. 
#
# You could also start a local cluster that just uses the cores available on the computer running the notebook server, and there are [many other ways to set up Dask clusters](https://docs.dask.org/en/latest/setup.html).
#

# %%
# from dask_gateway import Gateway
# gateway = Gateway()

# %%
# options = gateway.cluster_options()
# options

# %%
# options.environment='pangeo'

# %%
# # %%time
# cluster = gateway.new_cluster(cluster_options=options)

# %%
# client = Client(cluster)

# %%
# client

# %%
#cluster.adapt(minimum=4, maximum=20);
# cluster.scale(10)

# %%
# Local dask client
client = Client()
client

# %%

# %% [markdown]
# Open Zarr datasets in Xarray using a mapper from fsspec.  We use `anon=True` for free-access public buckets like the AWS Open Data Program, and `requester_pays=True` for requester-pays public buckets. 

# %%
url = 's3://noaa-nwm-retro-v2-zarr-pds'

# %%
# %%time
ds = xr.open_zarr(fsspec.get_mapper(url, anon=True), consolidated=True)

# %%
var='streamflow'

# %%
ds[var]


# %%
def plot_nwm_field(da, label=None):
    # Convert Xarray to Pandas dataframe so we can use hvplot.points for visualization
    df = da.to_pandas().to_frame()
    #The dataframe just has streamflow, so add longitude and latitude as columns
    df = df.assign(latitude=ds['latitude'])
    df = df.assign(longitude=ds['longitude'])
    df.rename(columns={0: "transport"}, inplace=True)
    p = df.hvplot.points('longitude', 'latitude', geo=True,
                     c='transport', colorbar=True, size=14, label=label)
    # We don't want to plot all the 2.7M points individually, so aggregate 
    # to 0.02 degree resolution and rasterize with datashader. 
    # Use a log scale for visualization since there is a large dynamic range in streamflow.
    g = rasterize(p, aggregator='mean', x_sampling=0.02, y_sampling=0.02, width=500).opts(tools=['hover'], 
                aspect='equal', logz=True, cmap='viridis', clim=(1e-2, np.nan))
    return (g * gv.tile_sources.OSM)


# %% [markdown]
# ### Read and plot data for all the stations at a specific time

# %%
# %%time
select_time = '2017-06-01 00:00:00'
da = ds[var].sel(time=select_time)
plot_nwm_field(da, label=f'{var}:{select_time}')

# %% [markdown]
# ### Read and plot data for entire time series at a specific location 
# Just as an example we pick the location with the largest stream flow from the specific time above

# %%
# %%time
imax = da.argmax().values
ds[var][:,imax].hvplot(grid=True)

# %% [markdown]
# ### Compute mean discharge during May 2017 on all rivers

# %%
da= ds[var].sel(time=slice('2017-04-01 00:00','2017-05-01 00:00'))

# %%
da

# %%
# %%time
var_mean = da.mean(dim='time').compute()

# %%
plot_nwm_field(var_mean, 'Mean Streamflow: May 2017')

# %%
client.close(); cluster.shutdown()

# %%
