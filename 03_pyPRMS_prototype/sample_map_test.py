#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize

# from matplotlib.backends.backend_pdf import PdfPages

# from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import pyproj as prj
import os
from osgeo import ogr
import sys

# import pyPRMS.NhmParamDb as nhm
import pyPRMS.ParameterFile as nhm

# Create the colormap
# cmap = 'GnBu_r' # for snow
cmap = 'jet'

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
cblabel = 'soil2gw_max'

# Setup output to a pdf file
# outpdf = PdfPages('map_%s.pdf' % (cblabel))
outfile = 'map_%s.png' % cblabel

paramdb_dir = '/Users/pnorton/Projects/National_Hydrology_Model/paramDb/nhmparamdb'
param_file = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20170912_CONUS/myparam.param'

pdb = nhm.ParameterFile(param_file)
# pdb = nhm.NhmParamDb(paramdb_dir)

df = pdb.parameters.get_dataframe(cblabel)
#param_data = pdb.get(cblabel).data
#param_hrus = pdb.get('nhm_id').data

# Create a DataFrame of the parameter
#df = pd.DataFrame(param_data.reshape([1, param_data.size]), columns=param_hrus)

#df = pd.Series(param_data, index=param_hrus)

# Convert series to dataframe for merging with shapefile attributes
#df_a = pd.DataFrame(df)
#df_a.reset_index(inplace=True)
#df_a.rename(columns={'index': 'hru', 0: 'paramval'}, inplace=True)

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
    
# Output for pdf files
# outpdf.savefig()
# outpdf.close()

# Interactive plots
# plt.show()

