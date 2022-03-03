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
#     display_name: Python [conda env:gis_38]
#     language: python
#     name: conda-env-gis_38-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import cartopy.crs as ccrs

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm, PowerNorm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import shapely

from dask.distributed import Client
import geopandas
import numpy as np
import pandas as pd
import pyproj as prj
import os
import xarray

from osgeo import ogr

from pyPRMS.plot_helpers import set_colormap, get_projection, plot_polygon_collection


# %%
client = Client(n_workers=2, threads_per_worker=2, memory_limit='1GB')
# client = Client(processes=False)
client

# %%
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v3/tmp'
filename = '{}/crap.nc'.format(work_dir)


geofile = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
layer_name = 'nhruv11_sim30'
shape_key = 'nhru_v11'
# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v3_1980-01-01_2016-12-31/netcdf'
# filename = '{}/*.nc'.format(work_dir)

# %%

# %%
def shapefile_hrus(filename, layer_name=None):
    '''Read a shapefile or geodatabase that corresponds to HRUs
    '''

    hru_poly = geopandas.read_file(filename, layer=layer_name)

    if hru_poly.crs.name == 'USA_Contiguous_Albers_Equal_Area_Conic_USGS_version':
        print(f'Overriding USGS aea crs with EPSG:5070')
        hru_poly.crs = 'EPSG:5070'
#     self.__hru_shape_key = shape_key

    return hru_poly


def plot_cbh(name, hru_poly, shape_key, df, time_index=None, output_dir=None):
    '''
    Plot a parameter
    '''
    
    # Get extent information
    minx, miny, maxx, maxy = hru_poly.geometry.total_bounds

    # Get a Pandas dataframe from the netcdf dataset for a single timestep
    xdf_df = df[name][time_index].to_dataframe(name=name)

    # cmap, norm = set_colormap(name, param_var, min_val=cparam.minimum, max_val=cparam.maximum)
    cmap, norm = set_colormap(name, xdf_df[name])

    crs_proj = get_projection(hru_poly)

    # This takes care of multipolygons that are in the NHM geodatabase/shapefile
    geoms_exploded = hru_poly.explode().reset_index(level=1, drop=True)

    print('Writing first plot')
    df_mrg = geoms_exploded.merge(xdf_df[name], left_on=shape_key, right_index=True, how='left')

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(30, 20))

    ax = plt.axes(projection=crs_proj)
    ax.coastlines()
    ax.gridlines()
    ax.set_extent([minx, maxx, miny, maxy], crs=crs_proj)

    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapper.set_array(df_mrg[name])
    plt.colorbar(mapper, shrink=0.6)

    time_str = pd.to_datetime(df['time'][time_index].values).isoformat()
#     time_str = df['time'].iloc[time_index].isoformat()
    plt.title(f'Variable: {name},  Date: {time_str}')    

    col = plot_polygon_collection(ax, df_mrg.geometry, values=df_mrg[name], colormap=cmap,
                                  norm=norm, linewidth=0.0)

    if output_dir is not None:
#         if is_monthly:
#             plt.savefig(f'{output_dir}/{name}_{time_index+1:02}.png', dpi=150, bbox_inches='tight')

#             for tt in range(1, 12):
#                 print(f'    Index: {tt}')
#                 param_var = self.get_dataframe(name).iloc[:, tt].to_frame(name=name)
#                 df_mrg = geoms_exploded.merge(xdf_df[name], left_on=shape_key, right_index=True, how='left')

#                 if is_monthly:
#                     ax.set_title(f'Variable: {name},  Month: {tt+1}')

#                 col.set_array(df_mrg[name])
#                 # fig
#                 plt.savefig(f'{output_dir}/{name}_{tt+1:02}.png', dpi=150, bbox_inches='tight')
#         else:
        plt.savefig(f'{output_dir}/{name}_{time_index+1:02}.png', dpi=150, bbox_inches='tight')


# %%
# Open the netcdf file of NHM output variables
xdf = xarray.open_mfdataset(filename, chunks={'hruid': 1040}, combine='by_coords')
xdf = xdf.set_index(hruid='hruid')

# %%
hru_poly = shapefile_hrus(geofile, layer_name=layer_name)

# %%
plot_cbh('tmax', hru_poly, shape_key, xdf, time_index=365, output_dir='/Users/pnorton/tmp')

# %%
# # %%timeit -n1 -r1
# a = None
# 12.8 s - mostly taken up by the min/max of the variable data
# 1:05 - no dask
#  :56.7 - with dask 2w2t1G

# NOTE: dim changed from nhru to hru
# Get the min and max values for the variable
max_val = xdf[the_var].max(dim='time').max(dim='hru').values
min_val = xdf[the_var].min(dim='time').min(dim='hru').values

# norm = PowerNorm(gamma=0.05)
# norm = LogNorm(vmin=min_val, vmax=max_val)

if min_val == 0.:
    if the_var in ['net_ppt', 'net_rain', 'net_snow', 'prcp']:
        cmap.set_under(color='None')
        norm = LogNorm(vmin=0.000001, vmax=max_val)
    else:
        norm = Normalize(vmin=0.000001, vmax=max_val)
else:
    if the_var in ['tmax_hru', 'tmin_hru', 'tmax', 'tmin']:
        norm = Normalize(vmin=-max_val, vmax=max_val)
    else:
        norm = Normalize(vmin=min_val, vmax=max_val)

# %%
