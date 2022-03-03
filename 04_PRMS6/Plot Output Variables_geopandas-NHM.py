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
#     display_name: Python [conda env:bandit]
#     language: python
#     name: conda-env-bandit-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %% [markdown]
# # Plot by-HRU NHM output variables

# %%
import cartopy.crs as ccrs
from cartopy.io import shapereader
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm, PowerNorm
from matplotlib.collections import PatchCollection
import geopandas

from matplotlib.patches import Polygon
import shapely

import matplotlib as mpl

import pyproj as prj
from osgeo import ogr

import xarray

import netCDF4 as cf
import pandas as pd
from datetime import datetime
from calendar import monthrange
from collections import namedtuple
import numpy as np
import os

# %%
work_dir = '/Volumes/parker_rocks/NHM_output/netcdf'
filename = '{}/pkwater_equiv.nc'.format(work_dir)

# %%
# ### Get extent information from the national HRUs shapefile

# Need two shapefiles 1) in projected coordinates, 2) in geographic coordinates
# If gdal is installed can create geographic coordinates from projected with:
#   ogr2ogr -t_srs epsg:4326 output_wgs84.shp input.shp

# shpfile = '/Users/pnorton/Projects/National_Hydrology_Model/extraction_requests/20180307_red_river/GIS/HRU_subset_nad83.shp'
# shpfile_extent = '/Users/pnorton/Projects/National_Hydrology_Model/extraction_requests/20180307_red_river/GIS/HRU_subset_usaea.shp'

shpfile = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple/nhruNationalIdentifier.shp'
shpfile_extent = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp'
# shpfile_extent = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/NHM_GF_reduced.gdb'
# shpfile_extent = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_usaea/nhruNationalIdentifier.shp'

# Name of attribute to use. Change to match the name of the HRU id attribute in the shapefile
shape_key='hru_id_nat'

# Use gdal/ogr to get the extent information
# Shapefile can be in projected coordinates
# Driver can be: OpenFileGDB or ESRI Shapefile
inDriver = ogr.GetDriverByName("ESRI Shapefile")
inDataSource = inDriver.Open(shpfile_extent, 0)
inLayer = inDataSource.GetLayer()
# inLayer = inDataSource.GetLayerByName('nhruNationalIdentifier')
extent = inLayer.GetExtent()

# Get the spatial reference information from the shapefile
spatial_ref = inLayer.GetSpatialRef()

# Create transformation object using projection information from the shapefile 
xform = prj.Proj(spatial_ref.ExportToProj4())

west, east, south, north = extent
pad = 100000.    # amount to pad the extent values with (in meters)
#east += pad
#west -= pad
#south -= pad
#north += pad

LL_lon, LL_lat = xform(west, south, inverse=True)
UR_lon, UR_lat = xform(east, north, inverse=True) 
print('\tExtent: ({0:f}, {1:f}, {2:f}, {3:f})'.format(west, east, south, north))
print('\tExtent: (LL: [{}, {}], UR: [{}, {}])'.format(LL_lon, LL_lat, UR_lon, UR_lat))

extent_dms = [LL_lon, UR_lon, LL_lat, UR_lat]

# Matplotlib basemap requires the map center (lon_0, lat_0) be in decimal degrees
# and yet the corners of the extent can be in projected coordinates
cen_lon, cen_lat = xform((east+west)/2, (south+north)/2, inverse=True)

print('cen_lon: {}'.format(cen_lon))
print('cen_lat: {}'.format(cen_lat))

# %%
print(spatial_ref)

# %%
# Read the shapefile
hru_df = geopandas.read_file(shpfile_extent)


# %%
# Open the netcdf file of NHM output variables
xdf = xarray.open_mfdataset(filename, chunks={'hru': 1040})
# xdf = xdf.set_index(hru='hru')

# xdf2 = xdf.set_index(nhru='nhm_id')
# xdf_df = xdf2['pkwater_equiv'][2].to_dataframe()
# xdf_df.drop('time', axis=1, inplace=True)

# df_mrg = hru_df.merge(xdf_df, left_on='hru_id_nat', right_index=True, how='left')

# df_mrg.plot(column='pkwater_equiv', legend=True, cmap='jet')

# %%
xdf['hru']

# %%
# Create the colormap
# cmap = 'BrBG' #'GnBu_r' # for snow
# cmap = 'GnBu_r'
cmap = 'jet'

# create the colormap if a list of names is given, otherwise
# use the given colormap
lscm = mpl.colors.LinearSegmentedColormap
if isinstance(cmap,(list,tuple)):
    cmap = lscm.from_list('mycm', cmap)
else:
    cmap = plt.get_cmap(cmap)
    
missing_color = '#ff00cb'   # pink/magenta 

norm = PowerNorm(gamma=0.05)



# %%
def plot_polygon_collection(ax, geoms, values=None, colormap='Set1',  facecolor=None, edgecolor=None,
                            alpha=0.5, linewidth=1.0, **kwargs):
    """ Plot a collection of Polygon geometries """
    # from https://stackoverflow.com/questions/33714050/geopandas-plotting-any-way-to-speed-things-up
    patches = []

    for poly in geoms:

        a = np.asarray(poly.exterior)
        if poly.has_z:
            poly = shapely.geometry.Polygon(zip(*poly.exterior.xy))

        patches.append(Polygon(a))

    patches = PatchCollection(patches, facecolor=facecolor, linewidth=linewidth, edgecolor=edgecolor, alpha=alpha, **kwargs)

    if values is not None:
        patches.set_array(values)
        patches.set_cmap(colormap)

    ax.add_collection(patches, autolim=True)
    ax.autoscale_view()
    return patches


# %%
the_var = 'pkwater_equiv'
xdf_df = xdf[the_var][2].to_dataframe()
xdf_df.head()

# %%
the_var = 'pkwater_equiv'

print('Build dataset')
# This takes care of multipolygons that are in the NHM geodatabase/shapefile
geoms_exploded = hru_df.explode().reset_index(level=1, drop=True)

xdf_df = xdf[the_var][2].to_dataframe()
df_mrg = geoms_exploded.merge(xdf_df, left_on='hru_id_nat', right_index=True, how='left')

print('Setup plot area')
lcc_proj = ccrs.LambertConformal(central_latitude=40., standard_parallels=(20, 60))
# aea_proj = ccrs.AlbersEqualArea(central_longitude=-96., central_latitude=23., standard_parallels=(29.5, 45.5))

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(30,20))

ax = plt.axes(projection=lcc_proj)
ax.coastlines()
ax.gridlines()
ax.set_extent(extent_dms)

mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
mapper.set_array(df_mrg[the_var])
plt.colorbar(mapper, shrink=0.6)

plt.title('Variable: {},  Date: {}'.format(the_var, xdf_df['time'].iloc[0].isoformat()))

print('Call plot_polygon_collection()')
col = plot_polygon_collection(ax, df_mrg.geometry, values=df_mrg[the_var], colormap=cmap, norm=norm)

# %%
xdf_df = xdf[the_var][4100].to_dataframe()
df_mrg = geoms_exploded.merge(xdf_df, left_on='hru_id_nat', right_index=True, how='left')

ax.set_title('Variable: {},  Date: {}'.format(the_var, xdf_df['time'].iloc[0].isoformat()))
col.set_array(df_mrg[the_var])
fig

# %%

# %%
