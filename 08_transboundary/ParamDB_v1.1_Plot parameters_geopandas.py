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
#     display_name: Python [conda env:gis]
#     language: python
#     name: conda-env-gis-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import cartopy.crs as ccrs

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm, PowerNorm
from matplotlib.collections import PatchCollection

# NOTE: requires geopandas >= 0.7.0
import geopandas

from matplotlib.patches import Polygon
import shapely

import matplotlib as mpl

import pyproj as prj
from osgeo import ogr
import numpy as np

from pyPRMS.ParamDb import ParamDb
from xml.etree import cElementTree as ElementTree

# %%
work_dir = '/Users/pnorton/tmp/tmp_paramdb'
shpfile = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
layer_name = 'nhru_v11'

# Name of attribute to use. Change to match the name of the HRU id attribute in the shapefile
shape_key='nhru_v11'

pdb = ParamDb(paramdb_dir=work_dir, verbose=True, verify=True)


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
# Read the shapefile
hru_df = geopandas.read_file(shpfile, layer=layer_name)

# %%
if hru_df.crs.name == 'USA_Contiguous_Albers_Equal_Area_Conic_USGS_version':
    print('Overriding USGS crs')
    hru_df.crs = 'EPSG:5070'
    
print(hru_df.crs)

# Get extent information
minx, miny, maxx, maxy = hru_df.geometry.total_bounds

# %%
hru_df

# %%
# How to get the coordinate-related information

# hru_df.crs.prime_meridian
# print(hru_df.crs.coordinate_operation.to_wkt(pretty=True))
xx = hru_df.crs.coordinate_operation.params
print(xx)

# %%
# Select a variable to plot and starting time_index (if needed)
the_var = 'sro_to_dprst_imperv'
time_index = 0

# Use for nhru-dimensioned parameters
param_var = pdb.parameters.get_dataframe(the_var).iloc[:]

# Use for nhru-by-nmonths dimensioned parameters
# param_var = pdb.parameters.get_dataframe(the_var).iloc[:, time_index].to_frame(name=the_var)

param_var.head()

# %%
# Create the colormap
# cmap = 'BrBG' #'GnBu_r' # for snow
# cmap = 'GnBu_r'
# cmap = 'jet'

if the_var in ['tmax_allsnow', 'tmax_allrain_offset']:
    cmap = 'RdBu_r'
elif the_var in ['net_ppt', 'net_rain', 'net_snow']:
    cmap = 'YlGnBu'
elif the_var in ['tmax_cbh_adj', 'tmin_cbh_adj']:
    cmap = 'coolwarm'
elif the_var in ['dprst_frac']:
    cmap = 'tab20b'
else:
    cmap = 'jet'

# create the colormap if a list of names is given, otherwise
# use the given colormap
lscm = mpl.colors.LinearSegmentedColormap
if isinstance(cmap,(list,tuple)):
    cmap = lscm.from_list('mycm', cmap)
else:
    cmap = plt.get_cmap(cmap)
    
missing_color = '#ff00cb'   # pink/magenta 

# Get the min and max values for the variable
max_val = param_var.max().max()
min_val = param_var.min().min()

# Temporary override for tmax_allrain_offset
# min_val = 0.451507
# max_val = 10.334641

# norm = PowerNorm(gamma=0.05)
# norm = LogNorm(vmin=min_val, vmax=max_val)

if min_val == 0.:
    if the_var in ['net_ppt', 'net_rain', 'net_snow']:
        cmap.set_under(color='None')
        norm = LogNorm(vmin=0.000001, vmax=max_val)
    else:
        norm = Normalize(vmin=0.000001, vmax=max_val)
else:
    if the_var in ['tmax_allsnow', 'tmax_allrain_offset']:
        norm = Normalize(vmin=min_val, vmax=max_val)
    elif the_var in ['tmax_cbh_adj', 'tmin_cbh_adj']:
        norm = Normalize(vmin=-max_val, vmax=max_val)
    else:
        norm = Normalize(vmin=min_val, vmax=max_val)


# %%
print(min_val)
print(max_val)

# %%
# This takes care of multipolygons that are in the NHM geodatabase/shapefile
geoms_exploded = hru_df.explode().reset_index(level=1, drop=True)
# print(geoms_exploded)

# xdf_df = xdf[the_var][2].to_dataframe()
df_mrg = geoms_exploded.merge(param_var, left_on=shape_key, right_index=True, how='left')

aa = {}
for yy in hru_df.crs.coordinate_operation.params:
    aa[yy.name] = yy.value

if '9822' in hru_df.crs.coordinate_operation.method_code:
    # Albers Equal Area 
    crs_proj = ccrs.AlbersEqualArea(central_longitude=aa['Longitude of false origin'], 
                                    central_latitude=aa['Latitude of false origin'], 
                                    standard_parallels=(aa['Latitude of 1st standard parallel'], 
                                                        aa['Latitude of 2nd standard parallel']),
                                    false_easting=aa['Easting at false origin'],
                                    false_northing=aa['Northing at false origin'])
elif '9802' in hru_df.crs.coordinate_operation.method_code:
    # Lambert Conformal Conic
    crs_proj = ccrs.LambertConformal(central_latitude=aa['Latitude of false origin'], 
                                     central_longitude=aa['Longitude of false origin'],
                                     standard_parallels=(aa['Latitude of 1st standard parallel'], 
                                                         aa['Latitude of 2nd standard parallel']),
                                     false_easting=aa['Easting at false origin'],
                                     false_northing=aa['Northing at false origin'])
else:
    pass
    # We're gonna crash
    
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(30,20))

ax = plt.axes(projection=crs_proj)
ax.coastlines()
ax.gridlines()

ax.set_extent([minx, maxx, miny, maxy], crs=crs_proj)

mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
mapper.set_array(df_mrg[the_var])
plt.colorbar(mapper, shrink=0.6)

# Select one plot title (leave the others commented out)
# plt.title('Variable: {},  Date: {}'.format(the_var, xdf_df['time'].iloc[0].isoformat()))
plt.title(f'Variable: {the_var},  Month: {time_index+1}')
# plt.title('Variable: {}'.format(the_var))

col = plot_polygon_collection(ax, df_mrg.geometry, values=df_mrg[the_var], colormap=cmap, norm=norm, linewidth=0.0)

# %%
# After making the first plot you can plot other time steps much more
# quickly by adjusting the time_index and running this code.
# This only applies to parameters which are nmonths-by-nhru.
time_index = 7
param_var = pdb.parameters.get_dataframe(the_var).iloc[:,time_index].to_frame(name=the_var)
df_mrg = geoms_exploded.merge(param_var, left_on=shape_key, right_index=True, how='left')

# Select one title based on parameter type
ax.set_title(f'Variable: {the_var},  Month: {time_index+1}')
# ax.set_title('Variable: {},  Date: {}'.format(the_var, xdf_df['time'].iloc[0].isoformat()))

col.set_array(df_mrg[the_var])
fig

# %%
