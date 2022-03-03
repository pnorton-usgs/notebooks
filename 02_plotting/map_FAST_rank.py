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

import prms_lib

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

stats = ['amp', 'ar_1', 'cv', 'kurt', 'mean', 'phase', 'skew']
procs = ['gwres_flow', 'hru_actet', 'hru_outflow', 'infil', 
         'snowmelt', 'soil_most', 'sroff', 'ssres_flow']

cstat = 4
cproc = 4

# Get the national model HRU id's from the input parameter file
params = prms_lib.parameters('/Users/pnorton/Projects/National_Hydrology_Model/FAST/rPipestem/input/daymet.params.expanded')

thevar = params.get_var('nhm_id')['values']
df_hrus = pd.DataFrame(thevar, columns=['HRU'])

# Load the data to plot
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/FAST/rPipestem/output/sens'
filename = 'sens_{}_{}.csv'.format(stats[cstat], procs[cproc])

df = pd.read_csv('{}/{}'.format(workdir, filename))
df.index.name = 'HRU'

# Restrict to partial variances that explain at least 1% of the total variance
df_tmp = df[df >= 0.01]

# Rank the sensitive parameters by HRU, for parameters with sensitivity >= 1%
df_ranked = df_tmp.rank(axis=1, ascending=False)

by_param = df_ranked.unstack().unstack()

# Remove rows where all entries are NaN
by_param = by_param.dropna(axis=0, how='all', thresh=None, inplace=False).copy()

xx = by_param.unstack()
tmp1 = pd.DataFrame(xx[xx == 1.0]).reset_index().level_1
# tmp1.index.name = 'HRU'


rnk1 = pd.concat([df_hrus, tmp1], axis=1).reset_index().drop(['index'], axis=1).set_index(['HRU'])
print rnk1

rank_colors = {'cecn_coef': '#339933', 'dprst_depth_avg': '#336633',  
               'dprst_et_coef': '#cc33cc', 'dprst_frac': '#cc0000', 
               'dprst_seep_rate_open': '#cccc33' , 'emis_noppt': '#ff6600' ,
               'fastcoef_lin': '#9933ff', 'freeh2o_cap': '#333366', 
               'gwflow_coef': '#3399ff', 'potet_sublim': '#33cccc',
               'pref_flow_den': '#00cc66', 'radmax': '#6666ff',
               'slowcoef_lin': '#ff3300', 'snowinfil_max': '#ccff00', 
               'soil_moist_max': '#cccc99', 'sro_to_dprst_perv': '#ff99cc',
               'tmax_index': '#cc6699'}

# df = pd.DataFrame(thevar, columns=['value'])
# df.index.name = 'HRU'

cblabel = '{}_{}'.format(stats[cstat], procs[cproc])
missing_color = '#00BFFF'


# %%

# %%
# months = ['January', 'February', 'March', 'April', 'May', 'June',
#           'July', 'August', 'September', 'October', 'November', 'December']

# Name of shapefile
shpfile='/Users/pnorton/Projects/National_Hydrology_Model/notebooks/gis_pipestem/upper_pipestem_ll'

# Name of attribute to use
shape_key='hru_id'

# Setup output to a pdf file
outpdf = PdfPages('map_ranked_{}_{}.pdf'.format(stats[cstat], procs[cproc]))

#fig, ax = plt.subplots(1,figsize=(20,30))
#ax = plt.gca()
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(30,20))
ax = axes.flatten()
# ax = axes

extent = (47.6, 47, -98.5, -99.9)  # Extent for Pipestem creek watershed
# extent = (50, 42, -95, -114)  # Extent for r10U
north, south, east, west = extent

# min_val = df.min().min()
# max_val = df.max().max()
min_val = 0.0
max_val = 1.0

# ax.set_title(cblabel);

for rr in range(3):
    print 'rank:', rr
    # Get the current ranked results
    xx = by_param.unstack()
    tmp1 = pd.DataFrame(xx[xx == rr+1]).reset_index().level_1
    rnk1 = pd.concat([df_hrus, tmp1], axis=1).reset_index().drop(['index'], axis=1).set_index(['HRU'])
    Series_data = rnk1.iloc[:]
    
    # Load the basemap
    m = Basemap(llcrnrlon=west, llcrnrlat=south, urcrnrlon=east, urcrnrlat=north, resolution='c',
                projection='laea', lat_0=(south+north)/2, lon_0=(east+west)/2, ax=ax[rr])

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

        # Skip those HRUs that aren't in the dataset without complaint
        if index in Series_data.index:
            # Set the color for each region
            if pd.isnull(val):
                # Record exists but the value is NaN
                color = missing_color
            else:
                val = Series_data.loc[index].values[0]
                color = rank_colors[val]
        else:
            # The record is totally missing
            color = '#ff3300'
            print "Missing Value for nhru =", index, val

        # Extract the x and y of the countours and plot them
        xx, yy = zip(*nhruDd_borders)
        patches = ax[rr].fill(xx, yy, facecolor=color, edgecolor='grey')

    # Generate a synthetic colorbar starting from the maximum and minimum of the dataset
#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes("right", "2%", pad="3%")

#     #axc, kw = mpl.colorbar.make_axes(ax)
#     norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)

#     cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
#     cb1.set_label(cblabel)
    #cb1.set_ticks()
    #cb1.ax.tick_params(labelsize=20)
    ax[rr].patch.set_facecolor('0.93')
    
outpdf.savefig()
outpdf.close()
#plt.show()

# %%
