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

# %%
from future.utils import iteritems


import matplotlib as mpl
mpl.use('Agg')
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
import pyproj as prj
import os
from osgeo import ogr

import datetime
import pandas as pd
from collections import OrderedDict

import pyPRMS.ParameterFile as nhm
from pyPRMS.Cbh import Cbh
from pyPRMS.prms_helpers import dparse

CBH_VARNAMES = ['prcp', 'tmin', 'tmax']
CBH_INDEX_COLS = [0, 1, 2, 3, 4, 5]

# %%
cblabel = 'tmin'

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20170912_CONUS'
cbh_file = '{}/{}_onemonth.cbh'.format(workdir, cblabel)
param_file = '{}/myparam.param'.format(workdir)

# %%
st_date = datetime.datetime(1980, 1, 1)
en_date = datetime.datetime(1980, 1, 1)

# %%
# Read the nhm_id from the parameter file
pdb = nhm.ParameterFile(param_file)
df_nhm = pdb.parameters.get_dataframe('nhm_id')

# %%
# Read the CBH file
df_cbh = Cbh(filename=cbh_file, st_date=st_date, en_date=en_date)
df_cbh.read_cbh_full()

# Rename columns with NHM HRU ids
ren_dict = {k + 6: v for k, v in enumerate(df_nhm['nhm_id'].tolist())}

# NOTE: The rename is an expensive operation
df_cbh.data.rename(columns=ren_dict, inplace=True)

df_cbh.data.reset_index(drop=True, inplace=True)
df_cbh2 = df_cbh.data.transpose()
df_cbh2.rename(columns={0: cblabel}, inplace=True)
df_cbh2.head()

# %%

# %%
# Setup output to a pdf file
# outpdf = PdfPages('map_%s.pdf' % (cblabel))
outfile = 'map_%s.png' % cblabel

df = df_cbh2

# %%
# Create the colormap
cmap = 'jet' #'GnBu_r' # for snow

# create the colormap if a list of names is given, otherwise
# use the given colormap
lscm = mpl.colors.LinearSegmentedColormap
if isinstance(cmap,(list,tuple)):
    cmap = lscm.from_list('mycm', cmap)
else:
    cmap = plt.get_cmap(cmap)
    
missing_color = '#ff00cb'   # pink/magenta 

# ### Read in parameter file and create dataframe of selected parameter from the NHM parameter database
# #### <I>Min and max values for range are currently hardcoded or selected from the range of parameter values. <br>It would be better to read accepted range information from an xml file for the parameters.</I>

# TODO: Lookup dimensions for given parameter use that to select segment IDs or HRU IDs 
#       for the index column. If parameter is 2D (e.g. nhru x nmonths) then use the
#       second dimension to loop and create one plot/map for each value of the second
#       dimension.


# Set the min and max values allowed - right now we just take the min and max of the dataframe values
min_val = df.min().min()
max_val = df.max().max()
val_rng = max_val - min_val

# ### Get extent information from the national HRUs shapefile
shpfile = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple/nhruNationalIdentifier.shp'
shpfile_extent = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_usaea/nhruNationalIdentifier.shp'

# Name of attribute to use. Change to match the name of the HRU id attribute in the shapefile
shape_key='hru_id_nat'

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
pad = 100000.    # amount to pad the extent values with (in meters)
#east += pad
#west -= pad
#south -= pad
#north += pad

LL_lon, LL_lat = xform(west, south, inverse=True)
UR_lon, UR_lat = xform(east, north, inverse=True) 
print('\tExtent: ({0:f}, {1:f}, {2:f}, {3:f})'.format(west, east, south, north))
print('\tExtent: (LL: [{}, {}], UR: [{}, {}])'.format(LL_lon, LL_lat, UR_lon, UR_lat))

# Matplotlib basemap requires the map center (lon_0, lat_0) be in decimal degrees
# and yet the corners of the extent can be in projected coordinates
cen_lon, cen_lat = xform((east+west)/2, (south+north)/2, inverse=True)

# ### Create the map figure
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(30,20))
# ax = axes.flatten()
ax = axes

# Load the basemap
m = Basemap(width=east-west, height=north-south, resolution='c',
            projection='laea', lat_0=cen_lat, lon_0=cen_lon, ax=ax)
#m = Basemap(llcrnrlon=west, llcrnrlat=south, urcrnrlon=east, urcrnrlat=north, resolution='c',
            #projection='laea', lat_0=(south+north)/2, lon_0=(east+west)/2, ax=ax)

# draw parallels.
m.drawparallels(np.arange(0.,90,10.), labels=[1,0,0,0], fontsize=20)

# draw meridians
m.drawmeridians(np.arange(180.,360.,10.), labels=[0,0,0,1], fontsize=20)
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
df_poly = df_poly.merge(df, left_on='id', right_index=True, how='left')

print('Patch Collection')
pc = PatchCollection(df_poly.shapes, zorder=2)
norm = Normalize()

print('facecolor')
# set_color sets both the edgecolor and the facecolor
pc.set_color(cmap(norm(df_poly[cblabel].fillna(-99).values)))
ax.add_collection(pc)

print('mapping...')
mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
mapper.set_array(df_poly[cblabel])
plt.colorbar(mapper, shrink=0.4)

# Output a png
fig.savefig(outfile, dpi=250, bbox_inches='tight',debug=True)
fig.clf()
plt.close()


# %%
df.min()

# %%
df_poly.min()

# %%
# Read the CBH file
df_crap = Cbh(filename=cbh_file, st_date=st_date, en_date=en_date)
df_crap.read_cbh_full()

# %%
dd = df_crap.data
dd.head()

# %%
gg = dd.columns.tolist()
gg.sort()
print(gg)

# %%
len(gg)

# %%
len(set(gg))

# %%
ren_dict

# %%
