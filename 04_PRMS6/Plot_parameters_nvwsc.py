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
#     display_name: Python [conda env:prms_py3]
#     language: python
#     name: conda-env-prms_py3-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
# %matplotlib inline

# from __future__ import (absolute_import, division, print_function)
# from future.utils import iteritems

import pandas as pd
from datetime import datetime
from calendar import monthrange
from collections import namedtuple

import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize

# mpl.use('Agg')
# from matplotlib.backends.backend_pdf import PdfPages

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
import os
import pyproj as prj
from osgeo import ogr
import sys

from pyPRMS.ParameterFile import ParameterFile

# %% [markdown]
# ### Load the parameter file

# %%
# filename = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/03_pyPRMS_prototype/crap_new.param'
filename = '/Users/pnorton/Projects/National_Hydrology_Model/tmp/nvwsc/Desch.params6'

# Load parameter file
pfile = ParameterFile(filename)

# %% [markdown]
# ### Select parameter to plot

# %%
cparam = 'hru_strmseg_down_id'
var_df = pfile.parameters.get_dataframe(cparam)

# Set the min and max values allowed - right now we just take the min and max of the dataframe values
min_val = var_df.min().min()
max_val = var_df.max().max()
val_rng = max_val - min_val

cblabel = cparam

# %% [markdown]
# ### Get extent information and select the colormap

# %%
# Create the colormap
# See: https://matplotlib.org/tutorials/colors/colormaps.html
# cmap = 'BrBG' #'GnBu_r' # for snow
# cmap = 'GnBu_r'
# cmap = 'jet'
cmap = 'nipy_spectral'

# create the colormap if a list of names is given, otherwise
# use the given colormap
lscm = mpl.colors.LinearSegmentedColormap
if isinstance(cmap,(list,tuple)):
    cmap = lscm.from_list('mycm', cmap)
else:
    cmap = plt.get_cmap(cmap)
    
missing_color = '#ff00cb'   # pink/magenta 

# ### Get extent information from the national HRUs shapefile
# Need two shapefiles 1) in projected coordinates, 2) in geographic coordinates
# If gdal is installed can create geographic coordinates from projected with:
#   ogr2ogr -t_srs epsg:4326 output_wgs84.shp input.shp
shpfile = '/Users/pnorton/Projects/National_Hydrology_Model/tmp/nvwsc/GIS/hru_params_v3_wgs84.shp'
shpfile_extent = '/Users/pnorton/Projects/National_Hydrology_Model/tmp/nvwsc/GIS/hru_params_v3.shp'

# Name of attribute to use. Change to match the name of the HRU id attribute in the shapefile
shape_key = 'HRU_ID'

# Use gdal/ogr to get the extent information
# Shapefile can be in projected coordinates
# Driver can be: OpenFileGDB or ESRI Shapefile
inDriver = ogr.GetDriverByName("ESRI Shapefile")
inDataSource = inDriver.Open(shpfile_extent, 0)
inLayer = inDataSource.GetLayer()
extent = inLayer.GetExtent()

# Get the spatial reference information from the shapefile
spatial_ref = inLayer.GetSpatialRef()

# Create transformation object using projection information from the shapefile 
xform = prj.Proj(spatial_ref.ExportToProj4())

west, east, south, north = extent

# If padding to the extent boundaries as needed
pad = 10000.    # amount to pad the extent values (in meters)
east += pad
west -= pad
south -= pad
north += pad

LL_lon, LL_lat = xform(west, south, inverse=True)
UR_lon, UR_lat = xform(east, north, inverse=True) 
print('\tExtent: ({0:f}, {1:f}, {2:f}, {3:f})'.format(west, east, south, north))
print('\tExtent: (LL: [{}, {}], UR: [{}, {}])'.format(LL_lon, LL_lat, UR_lon, UR_lat))

# Matplotlib basemap requires the map center (lon_0, lat_0) be in decimal degrees
# and yet the corners of the extent can be in projected coordinates
cen_lon, cen_lat = xform((east+west)/2, (south+north)/2, inverse=True)

# %% [markdown]
# ### Plot the parameter on the map

# %%
if var_df.shape[1] > 1:
    print('Currently unable to handle 2D parameters')
else:
    # Create the map figure
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(30,20))
    # ax = axes.flatten()
    ax = axes

    # Load the basemap
    m = Basemap(width=east-west, height=north-south, resolution='c',
                projection='laea', lat_0=cen_lat, lon_0=cen_lon, ax=ax)

    # Add basemap features; see: https://matplotlib.org/basemap/users/geography.html
    m.drawstates()
    
    # draw parallels.
    m.drawparallels(np.arange(0.,90,2.), labels=[1,0,0,0], fontsize=20)

    # draw meridians
    m.drawmeridians(np.arange(180.,360.,2.), labels=[0,0,0,1], fontsize=20)
    m.drawmapboundary()

    # ------------------------------------------------------------------
    # Use basemap to read and draw the shapefile
    # Two variables are added to the basemap, m.nhruDd and m.nhruDd_info
    #     m.nhruDd contains the lines of the borders
    #     m.nhruDd_info contains the info on the hru, like the name
    print('Read shapefile...')
    m.readshapefile(os.path.splitext(shpfile)[0], 'nhruDd', drawbounds=False)

    print('Create dataframe')
    df_poly = pd.DataFrame({'shapes': [Polygon(np.array(ss), closed=True) for ss in m.nhruDd],
                            'id': [aa[shape_key] for aa in m.nhruDd_info]})
    df_poly = df_poly.merge(var_df, left_on='id', right_index=True, how='left')

    print('Patch Collection')
    pc = PatchCollection(df_poly.shapes, zorder=2)
    norm = Normalize(vmin=min_val, vmax=max_val)

    print('facecolor')
    pc.set_facecolor(cmap(norm(df_poly[cblabel].fillna(-99).values)))
    pc.set_edgecolor('face')
    ax.add_collection(pc)

    print('mapping...')
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapper.set_array(df_poly[cblabel])
    plt.colorbar(mapper, shrink=0.4)

    # Output a png
    # fig.savefig(outfile, dpi=250, bbox_inches='tight',debug=True)
    # fig.clf()
    # plt.close()

    plt.show()

# %%
var_df.tail()

# %%
a = 1
print('{:02d}'.format(a))

# %%
