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
from matplotlib.collections import LineCollection, PatchCollection

# NOTE: requires geopandas >= 0.7.0
import geopandas as gpd

from matplotlib.patches import Polygon
import shapely

import matplotlib as mpl

import pyproj as prj
from osgeo import ogr
import numpy as np

from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.plot_helpers import get_projection

from xml.etree import cElementTree as ElementTree

# %%
# work_dir = '/Users/pnorton/tmp/tmp_paramdb'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/2020-04-15_red_river_v2'
filename = f'{workdir}/myparam.param'

# shpfile = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
# layer_name = 'nsegment_v11'

# Name of attribute to use. Change to match the name of the HRU id attribute in the shapefile
# shape_key='nsegment_v11'

geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
layer_name = 'nsegment_v11'
seg_shape_key = 'nsegment_v11'


# %%
def plot_geopandas_as_collection(ax, geoms, colormap='Set1',  facecolor=None, edgecolor=None,
                                 alpha=0.5, linewidth=1.0, **kwargs):
    """ Plot a collection of Polygon geometries """
    # from https://stackoverflow.com/questions/33714050/geopandas-plotting-any-way-to-speed-things-up
    patches = []

    for poly in geoms:
        a = np.asarray(poly.exterior)
        if poly.has_z:
            a = shapely.geometry.Polygon(zip(*poly.exterior.xy))

        patches.append(Polygon(a))

    patches = PatchCollection(patches, facecolors=facecolor, linewidth=linewidth, edgecolors=edgecolor,
                              alpha=alpha, **kwargs)

    ax.add_collection(patches, autolim=True)
    ax.autoscale_view()
    return patches


# %%
def plot_line_collection(ax, geoms, values=None, vary_width=False, vary_color=True, colors=None, colormap='Set1', edgecolor=None,
                         alpha=0.5, linewidth=1.0, **kwargs):
    """ Plot a collection of line geometries """
    lines = []
    for geom in df_mrg.geometry:
        a = np.asarray(geom.coords)

        if geom.has_z:
            a = shapely.geometry.LineString(zip(*geom.xy))

        lines.append(shapely.geometry.LineString(a))

    if vary_width:
        lwidths = ((values / values.max()).to_numpy() + 0.01) * linewidth     
        lines = LineCollection(lines, linewidths=lwidths, colors=colors, alpha=alpha, **kwargs)
    elif vary_color:
        lines = LineCollection(lines, linewidth=linewidth, alpha=alpha, **kwargs)
        
    if vary_color and values is not None:
        lines.set_array(values)
        lines.set_cmap(colormap)

    ax.add_collection(lines, autolim=True)
    ax.autoscale_view()
    return lines


# %%
pdb = ParameterFile(filename, verbose=True, verify=True)
# pdb = ParamDb(paramdb_dir=work_dir, verbose=True, verify=True)

# %%
# Read the segment shapefile
# f'{workdir}/GIS/HRU_subset.shp', layer_name=None, shape_key='nhru_v11'
seg_shape_key = 'nsegment_v'
gdf = gpd.read_file(f'{workdir}/GIS/Segments_subset.shp')

minx, miny, maxx, maxy = gdf.geometry.total_bounds
# gdf = gpd.read_file(geodatabase, layer=layer_name)

# %%
# Read the HRU shapefile
hru_shape_key = 'nhru_v11'
hru_poly = gpd.read_file(f'{workdir}/GIS/HRU_subset.shp')

minx, miny, maxx, maxy = hru_poly.geometry.total_bounds

hru_geoms_exploded = hru_poly.explode().reset_index(level=1, drop=True)

# %%
# Load the parameter data
param_name = 'seg_length'
param_var = pdb.parameters.get_dataframe(param_name).iloc[:]
param_var.head()

# %%
# Merge the parameter data with the geospatial information
df_mrg = gdf.merge(param_var, left_on=seg_shape_key, right_index=True, how='left')

# %%
cmap = 'jet'

max_val = param_var.max().max()
min_val = param_var.min().min()
print(f'min: {min_val}; max: {max_val}')

norm = Normalize(vmin=min_val, vmax=max_val)

# %%
crs_proj = get_projection(gdf)
print(crs_proj)

# %%

# %%
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(30, 20))

ax = plt.axes(projection=crs_proj)
ax.coastlines()
ax.gridlines()
ax.set_extent([minx, maxx, miny, maxy], crs=crs_proj)

mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
mapper.set_array(df_mrg[param_name])

plt.colorbar(mapper, shrink=0.6)

plt.title('Variable: {}'.format(param_name))

hru_poly = plot_geopandas_as_collection(ax, hru_geoms_exploded.geometry, colormap=cmap, norm=norm, linewidth=0.5,
                                        facecolor='none', edgecolor='silver')

col = plot_line_collection(ax, df_mrg.geometry, values=df_mrg[param_name], vary_width=False,
                           vary_color=True, colors='blue', colormap=cmap, norm=norm, 
                           alpha=1.0, linewidth=3.0)


# %%

# %%
def func1(a=False, **kwargs):
    print(f'func1: a={a}')
    print(kwargs)
    func2(**dict(kwargs, world='goodbye', globe='big'))
    print('back in func1')
    print(kwargs)
    
def func2(b=True, world='Hello', **kwargs):
    print(f'func2: b={b}; world={world}')
    print(kwargs)


# %%
func1()

# %%
func1(b=False, john='dude', world='spasm')

# %%

# %%
