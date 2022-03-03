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
#     display_name: Python [conda env:idp_bandit]
#     language: python
#     name: conda-env-idp_bandit-py
# ---

# %%
# # %matplotlib inline

import matplotlib as mpl
mpl.use('Agg')
# from matplotlib.backends.backend_pdf import PdfPages

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import os
from osgeo import ogr
import sys

import pyPRMS.NhmParamDb as nhm
reload(nhm)

# %%
# Create the colormap
# cmap = 'YlGnBu_r'
cmap = 'GnBu_r' # for snow
#cmap = 'OrRd'  # for liquid
#cmap = 'seismic'
#cmap = ['Green', 'W','Sienna']

# create the colormap if a list of names is given, otherwise
# use the given colormap
lscm = mpl.colors.LinearSegmentedColormap
if isinstance(cmap,(list,tuple)):
    cmap = lscm.from_list('mycm', cmap)
else:
    cmap = plt.get_cmap(cmap)
    
missing_color = '#ff00cb'   # pink/magenta 

# %% [markdown]
# ### Read in parameter file and create dataframe of selected parameter from the NHM parameter database

# %% [markdown]
# #### <I>Min and max values for range are currently hardcoded or selected from the range of parameter values. <br>It would be better to read accepted range information from an xml file for the parameters.</I>

# %%
# TODO: Lookup dimensions for given parameter use that to select segment IDs or HRU IDs 
#       for the index column. If parameter is 2D (e.g. nhru x nmonths) then use the
#       second dimension to loop and create one plot/map for each value of the second
#       dimension.
cblabel = 'gwflow_coef'

# Setup output to a pdf file
# outpdf = PdfPages('map_%s.pdf' % (cblabel))
outfile = 'map_%s.png' % cblabel

paramdb_dir = '/Users/pnorton/Projects/National_Hydrology_Model/paramDb/nhmparamdb'

pdb = nhm.NhmParamDb(paramdb_dir)

param_data = pdb.get(cblabel).data
param_hrus = pdb.get('nhm_id').data

# Create a DataFrame of the parameter
df = pd.DataFrame(param_data.reshape([1, param_data.size]), columns=param_hrus)

df = pd.Series(param_data, index=param_hrus)
df_dict = df.to_dict()

# Set the min and max values allowed - right now we just take the min and max of the dataframe values
min_val = df.min().min()
max_val = df.max().max()
val_rng = max_val - min_val
# min_val = 0.0
# max_val = 1.0

# %%
param_data = pdb.parameters.get(cblabel).data
param_hrus = pdb.parameters.get('nhm_id').data

# Create a DataFrame of the parameter
df = pd.DataFrame(param_data.reshape([1, param_data.size]), columns=param_hrus)

df = pd.Series(param_data, index=param_hrus)
df_dict = df.to_dict()

# Set the min and max values allowed - right now we just take the min and max of the dataframe values
min_val = df.min().min()
max_val = df.max().max()
val_rng = max_val - min_val

# %%
print(df[1000])

# %% [markdown]
# ### Get extent information from the national HRUs shapefile

# %%
shpfile = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple/nhruNationalIdentifier.shp'
# shpfile = ''
# Name of attribute to use. Change to match the name of the HRU id attribute in the shapefile
shape_key='hru_id_nat'

# Use gdal/ogr to get the extent information
# Shapefile must must unprojected and supply lat/lon values
# Driver can be: OpenFileGDB or ESRI Shapefile
inDriver = ogr.GetDriverByName("ESRI Shapefile")
inDataSource = inDriver.Open(shpfile, 0)
inLayer = inDataSource.GetLayer()
extent = inLayer.GetExtent()

west, east, south, north = extent
pad = 0.02    # amount to pad the extent values with
east += pad
west -= pad
south -= pad
north += pad
print('\tExtent: ({0:f}, {1:f}, {2:f}, {3:f})'.format(west, east, south, north))

# %% [markdown]
# ### Create the map figure

# %%
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(30,20))
# ax = axes.flatten()
ax = axes

# Load the basemap
m = Basemap(llcrnrlon=west, llcrnrlat=south, urcrnrlon=east, urcrnrlat=north, 
            resolution='c', projection='laea', 
            lat_0=(south+north)/2, lon_0=(east+west)/2, ax=ax)

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
m.readshapefile(os.path.splitext(shpfile)[0], 'nhruDd', drawbounds=False)

# tmp_idx = df.index.tolist()
# tmpdict = {el:0 for el in tmp_idx}

for nhruDd_borders, nhruDd_info in zip(m.nhruDd, m.nhruDd_info):
    index = nhruDd_info[shape_key]
#     index -= 1    # Only needed when the hru ids in the parameter file are zero-based
    
    sys.stdout.write('\rIndex: {}'.format(index))
    sys.stdout.flush()
    
    #skip those that aren't in the dataset without complaints
    if index in df_dict:
#     if index in df.index:
        #set the color for each region
#         val = df.loc[index].values[0]
#         val = df[index]
        val = df_dict[index]

        if pd.isnull(val):
            # Record exists but the value is NaN
            color = missing_color
        elif val > max_val:
            # set color for patch to highlight HRUs with parameter value
            # out of range.
            pass
#             color = '#FF1493'
#             print "> max_val:", index, val
        else:
            color = cmap((val - min_val) / val_rng)
    else:
        # The record is totally missing. This shouldn't happen but sometimes it does.
        # This highlights those HRUs with missing information.
        color = '#ff3300'    # Red
        print "Missing Value for nhru =", index

    #extract the x and y of the countours and plot them
    xx, yy = zip(*nhruDd_borders)
    ax.fill(xx, yy, facecolor=color, edgecolor=color)
    #     patches = ax.fill(xx, yy, facecolor=color, edgecolor=color)


#generate a synthetic colorbar starting from the maximum and minimum of the dataset
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", "2%", pad="3%")

#axc, kw = mpl.colorbar.make_axes(ax)
norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)

cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
cb1.set_label(cblabel, fontsize=18)
# cb1.set_ticks()
# cb1.ax.tick_params(labelsize=20)
ax.patch.set_facecolor('0.93')

# Output a png
fig.savefig(outfile, dpi=250, bbox_inches='tight',debug=True)
fig.clf()
plt.close()
    
# Output for pdf files
# outpdf.savefig()
# outpdf.close()

# Interactive plots
# plt.show()

# %%

# %%
m.nhruDd_info[0]

# %%
m.nhruDd_info[0:-1]

# %%
