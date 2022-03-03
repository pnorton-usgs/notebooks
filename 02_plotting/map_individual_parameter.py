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
#     display_name: Python [default]
#     language: python
#     name: python2
# ---

# %%

# %%
# #%matplotlib inline

import matplotlib as mpl
mpl.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

# import prms_lib
import pyPRMS.ParameterFile as pf

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

# %%
# %%time
# Load the data to plot
# params = prms_lib.parameters('/media/scratch/PRMS/regions/r10U/input/params/daymet.params')
params = prms_lib.parameters('/media/scratch/PRMS/regions/r10U/daymet.params.expanded2')

thevar = params.get_var('sro_to_dprst_imperv')['values']
# df = pd.DataFrame(thevar, columns=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])
df = pd.DataFrame(thevar, columns=['value'])
df.index.name = 'HRU'

cblabel = 'sro_to_dprst_imperv'
missing_color = '#00BFFF'


# %%
# df = pd.DataFrame(thevar, columns=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])
# df.index.name = 'HRU'

# df = df[df > 0.0]
print df.min().min()
print df.max().max()
print df.loc[0].value
print df.head()

# %%
# months = ['January', 'February', 'March', 'April', 'May', 'June',
#           'July', 'August', 'September', 'October', 'November', 'December']

# Name of shapefile
shpfile='/media/scratch/PRMS/notebooks/nhru_10U/nhru_10U_simpl'

# Name of attribute to use
shape_key='hru_id_reg'

# Setup output to a pdf file
outpdf = PdfPages('map_%s.pdf' % (cblabel))

#fig, ax = plt.subplots(1,figsize=(20,30))
#ax = plt.gca()
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(30,20))
# ax = axes.flatten()
ax = axes

extent = (50, 42, -95, -114)  # Extent for r10U
north, south, east, west = extent

# min_val = df.min().min()
# max_val = df.max().max()
min_val = 0.0
max_val = 1.0

# for mm in xrange(12):
# print months[mm]
# sys.stdout.flush()

# Subset the data
# Series_data = df.iloc[:,mm]
Series_data = df.iloc[:]

ax.set_title(cblabel);

#     print "Loading basemap..."
# Load the basemap
m = Basemap(llcrnrlon=west,llcrnrlat=south,urcrnrlon=east,urcrnrlat=north, resolution='c',
            projection='laea', lat_0=(south+north)/2, lon_0=(east+west)/2, ax=ax)

# draw parallels.
m.drawparallels(np.arange(0.,90,10.), labels=[1,0,0,0], fontsize=20)

# draw meridians
m.drawmeridians(np.arange(180.,360.,10.), labels=[0,0,0,1], fontsize=20)
m.drawmapboundary()

# ------------------------------------------------------------------
# use basemap to read and draw the shapefile
# Two variables are added to the basemap, m.nhruDd and m.nhruDd_info
#     print 'Loading shapefile...'
m.readshapefile(shpfile,'nhruDd',drawbounds=False);

#     print 'Color HRUs...'
# m.nhruDd contains the lines of the borders
# m.nhruDd_info contains the info on the hru, like the name
for nhruDd_borders, nhruDd_info in zip(m.nhruDd, m.nhruDd_info):
    index = nhruDd_info[shape_key]
    index -= 1
    
    #skip those that aren't in the dataset without complaints
    if index in Series_data.index:
        #set the color for each region
        val = Series_data.loc[index].value

        if pd.isnull(val):
            # Record exists but the value is NaN
            color = missing_color
        elif val > max_val:
            color = '#FF1493'
            print "> max_val:", index, val
        else:
            color = cmap((val - min_val) / (max_val - min_val))
    else:
        # The record is totally missing
        color = '#ff3300'
        print "Missing Value for nhru =", index, val

    #extract the x and y of the countours and plot them
    xx, yy = zip(*nhruDd_borders)
    patches = ax.fill(xx, yy, facecolor=color, edgecolor=color)

#generate a synthetic colorbar starting from the maximum and minimum of the dataset
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", "2%", pad="3%")

#axc, kw = mpl.colorbar.make_axes(ax)
norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)

cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
cb1.set_label(cblabel)
#cb1.set_ticks()
#cb1.ax.tick_params(labelsize=20)
ax.patch.set_facecolor('0.93')
    
outpdf.savefig()
outpdf.close()
#plt.show()

# %%
