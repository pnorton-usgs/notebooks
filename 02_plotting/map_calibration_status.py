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
#     display_name: Python 2
#     language: python
#     name: python2
# ---

# %%

# %%
# %matplotlib inline

import matplotlib as mpl
#mpl.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob
import numpy as np
import os
import pandas as pd
import sys

import prms_lib
import prms_cfg

# %%
configfile = '/media/scratch/PRMS/calib_runs/pipestem_1/basin.cfg'
basinConfigFile = 'basin.cfg'
runid = '2016-04-05_1017'

cblabel = 'Calibration_status'

cfg = prms_cfg.cfg(configfile)
base_calib_dir = cfg.base_calib_dir
print os.getcwd()

# Read the basins_file. Single basin per line
bfile = open(cfg.get_value('basins_file'), "r")
basins = bfile.read().splitlines()
bfile.close()
print 'Total basins:', len(basins)

# success list
s_list = [el for el in glob.glob('%s/*/runs/%s/.success' % (base_calib_dir, runid))]
s_list += [el for el in glob.glob('%s/*/runs/%s/.warning' % (base_calib_dir,runid))]
s_list += [el for el in glob.glob('%s/*/runs/%s/.error' % (base_calib_dir,runid))]
s_list += [el for el in glob.glob('%s/*/runs/%s/.retry' % (base_calib_dir,runid))]

b_list = []
stat = {'.success': 1, '.warning': 2, '.error': 3, '.retry': 4}
for ee in s_list:
    tmp = ee.split('/')
    ri = tmp.index('runs')
    hru = int(tmp[ri-1].split('_')[1]) + 1

    b_list.append([hru, stat[tmp[-1]]])
    
df = pd.DataFrame(b_list, columns=['HRU', 'status'])
df.sort(columns=['HRU'], inplace=True)
df.set_index(['HRU'], inplace=True)

# %%
# Name of shapefile
shpfile = '/media/scratch/PRMS/notebooks/shapefiles/upper_pipestem_ll'
#shpfile='/media/scratch/PRMS/notebooks/nhru_10U/nhru_10U_simpl'

# Name of attribute to use
shape_key='hru_id_loc'
# shape_key='hru_id_reg'

# Setup output to a pdf file
# outpdf = PdfPages('%s_map_v%d.pdf' % (cblabel, theversion))

fig, ax = plt.subplots(1,figsize=(20,30))

extent = (47.6, 47, -98.5, -99.9)  # Extent for Pipestem creek watershed
#extent = (50, 42, -95, -114)  # Extent for r10U
north, south, east, west = extent

#     sys.stdout.flush()
    
# Subset the data
Series_data = df.iloc[:]

ax.set_title('Calibration status');

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

max_val = 4.
min_val = 0.

#     print 'Color HRUs...'
# m.nhruDd contains the lines of the borders
# m.nhruDd_info contains the info on the hru, like the name
for nhruDd_borders, nhruDd_info in zip(m.nhruDd, m.nhruDd_info):
    index = nhruDd_info[shape_key]

    #skip those that aren't in the dataset without complaints
    if index in df.index:
        #set the color for each region
        val = df.loc[index].values[0]

        if pd.isnull(val):
            # Record exists but the value is NaN
            color = missing_color
        elif val == 1:
            color = '#008000'
        elif val == 2:
            color = '#FF8C00'
        elif val == 3:
            color = '#FF0000'
        elif val == 4:
            color = '#FF1493'
    else:
        # The record is totally missing
        color = '#C0C0C0'

    #extract the x and y of the countours and plot them
    xx, yy = zip(*nhruDd_borders)
    patches = ax.fill(xx, yy, facecolor=color, edgecolor='grey')

#generate a synthetic colorbar starting from the maximum and minimum of the dataset
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", "2%", pad="3%")

# Create a colormap of discrete colors
cmap = mpl.colors.ListedColormap(['#c0c0c0','#008000','#ff8c00','#ff0000','#FF1493'])
bounds = [0,1,2,3,4,5]

cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, boundaries=bounds, spacing='uniform', ticks=[0.5, 1.5, 2.5, 3.5, 4.5])

# Add colorbar tick label text, provides a lot of control over placement
cb1.ax.get_yaxis().set_ticks([])
thelabels = ['Pending/Running', 'Success', 'Warning', 'Error', 'Cancelled']
for j, lab in enumerate(thelabels):
    cb1.ax.text(.5, (2 * j + 1) / float(len(thelabels)*2), lab, ha='center', va='center', weight='bold', color='white', rotation=90)

# Adjust the placement of the colorbar text
# cb1.ax.get_yaxis().labelpad = -20

# Change the tick labels and rotate them
# Doesn't seem to provide control over centering vertically and horizontally
#cb1.ax.set_yticklabels(['Pending/Running', 'Success', 'Warning', 'Error'], va='center', ha='left', rotation=90)

cb1.set_label(cblabel)
ax.patch.set_facecolor('0.93')
    
# outpdf.savefig()
# outpdf.close()
plt.show()

# %%
