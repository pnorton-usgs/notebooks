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

import numpy as np
import os
from osgeo import ogr
import pandas as pd
import sys


# %%
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


# %%
runid = '2016-04-06_1609'
output_type = 'png'  # one of pdf or png

output_filename = '/media/scratch/PRMS/calib_runs/pipestem_1/%s_map_pbias' % runid

# Name of shapefile
shpfile = '/media/scratch/PRMS/notebooks/shapefiles/upper_pipestem_ll.shp'
#shpfile='/media/scratch/PRMS/notebooks/nhru_10U/nhru_10U_simpl'

# Name of attribute to use
shape_key='hru_id_loc'
#shape_key='hru_id_reg'

# Extent information (set to None to automatically detect)
#   Order = east, west, south, north
extent = None

# %%
if not extent or len(extent) != 4:
    print('Using extent information from shapefile')

    # Use gdal/ogr to get the extent information
    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(shpfile, 0)
    inLayer = inDataSource.GetLayer()
    extent = inLayer.GetExtent()

west, east, south, north = extent
pad = 0.02
east += pad
west -= pad
south -= pad
north += pad
print('\tExtent: (%f, %f, %f, %f)' % (west, east, south, north))    

# best_ofs = ['OF_AET']
# ofs_to_plot = ['OF_AET', 'OF_SWE', 'OF_runoff']

best_ofs = ['OF_AET', 'OF_SWE', 'OF_runoff', 'OF_comp']
ofs_to_plot = ['OF_AET', 'OF_SWE', 'OF_runoff']

# Load the information about the best pareto sets
# df = pd.read_csv('/media/scratch/PRMS/notebooks/02_test_and_prototype/%s_best.csv' % runid, index_col=0)
df = pd.read_csv('/media/scratch/PRMS/calib_runs/pipestem_1/%s_best.csv' % runid, index_col=0)
vmin = df[['OF_AET', 'OF_runoff', 'OF_SWE']].min().min()
vmax = df[['OF_AET', 'OF_runoff', 'OF_SWE']].max().max()
vmin -= (vmin % 50.)
vmax -= (vmax % 50.)

# Load the data to plot
# varname = 'Composite'
# best_of = 'OF_comp'
# of_to_plot = 'OF_AET'
# theversion = 1



# %%
# Create the colormap
cmap = 'seismic'

# create the colormap if a list of names is given, otherwise
# use the given colormap
# lscm = mpl.colors.LinearSegmentedColormap
# if isinstance(cmap,(list,tuple)):
#     cmap = lscm.from_list('mycm', cmap)
# else:
#     cmap = plt.get_cmap(cmap)
    

midp = 1.0 - (vmax / (vmax + abs(vmin)))
print midp
shifted_cmap = shiftedColorMap(mpl.cm.seismic, midpoint=midp, name='shifted')

cblabel = 'Percent bias'

missing_color = '#00BFFF'  # for missing values

# %%

if output_type == 'pdf':
    # Setup output to a pdf file
    outpdf = PdfPages('%s.pdf' % output_filename)

    fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(30,20))
    ax = axes.flatten()

axii = 0
for jj, curr_best in enumerate(best_ofs):
    print 'curr_best:', curr_best
    sys.stdout.flush()
    
    df2 = df[df['best'] == curr_best].copy()
    df2.reset_index(inplace=True)
    df2.set_index('HRU', inplace=True)
    
    for ii, curr_of in enumerate(ofs_to_plot):
        print '\tcurr_of:', curr_of
        sys.stdout.flush()
        
        if output_type == 'png':
            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6,6))
            ax = [axes]
            
        ax[axii].set_title('%s Percent Bias for %s' % (curr_of, curr_best));

        df3 = df2[curr_of]
        #     print "Loading basemap..."
        # Load the basemap
        m = Basemap(llcrnrlon=west, llcrnrlat=south, urcrnrlon=east, urcrnrlat=north, resolution='c',
                    projection='laea', lat_0=(south+north)/2, lon_0=(east+west)/2, ax=ax[axii])
        
        m.drawlsmask(land_color='#666666', ocean_color='aqua', lakes=True)
        
        # draw parallels.
        m.drawparallels(np.arange(0.,90,10.), labels=[1,0,0,0], fontsize=20)

        # draw meridians
        m.drawmeridians(np.arange(180.,360.,10.), labels=[0,0,0,1], fontsize=20)
        m.drawmapboundary()

        # ------------------------------------------------------------------
        # use basemap to read and draw the shapefile
        # Two variables are added to the basemap, m.nhruDd and m.nhruDd_info
        #     print 'Loading shapefile...'
        m.readshapefile(os.path.splitext(shpfile)[0], 'nhruDd', drawbounds=False);

        # find minimum and maximum of the dataset to normalize the colors
#         max_val = 200.
#         min_val = -200.

        #     print 'Color HRUs...'
        # m.nhruDd contains the lines of the borders
        # m.nhruDd_info contains the info on the hru, like the name
        for nhruDd_borders, nhruDd_info in zip(m.nhruDd, m.nhruDd_info):
            index = nhruDd_info[shape_key]

            #skip those that aren't in the dataset without complaints
            if index in df3.index:
                #set the color for each region
                val = df3.loc[index]

                if pd.isnull(val):
                    # Record exists but the value is NaN
                    color = missing_color
                else:
#                     color = shifted_cmap((val - min_val) / (max_val - min_val))
                    color = shifted_cmap((val - vmin) / (vmax - vmin))
            else:
                # The record is totally missing
                color = '#c0c0c0'

            #extract the x and y of the countours and plot them
            xx, yy = zip(*nhruDd_borders)
            patches = ax[axii].fill(xx, yy, facecolor=color, edgecolor='#666666')

        #generate a synthetic colorbar starting from the maximum and minimum of the dataset
        divider = make_axes_locatable(ax[axii])
        cax = divider.append_axes("right", "2%", pad="3%")

        #axc, kw = mpl.colorbar.make_axes(ax)
#         norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        cb1 = mpl.colorbar.ColorbarBase(cax, cmap=shifted_cmap, norm=norm)
        cb1.set_label(cblabel)
        #cb1.set_ticks()
        #cb1.ax.tick_params(labelsize=20)
        ax[axii].patch.set_facecolor('0.93')
        
        if output_type == 'png':
            plt.savefig('%s_%s_%s.png' % (output_filename, curr_best, curr_of), dpi=300)
#             plt.show()
            plt.close()
        else:    
            axii += 1

if output_type == 'pdf':
    outpdf.savefig()
    outpdf.close()

plt.show()

# %%
aa = df[['OF_AET', 'OF_runoff', 'OF_SWE']].min().min()
bb = df[['OF_AET', 'OF_runoff', 'OF_SWE']].max().max()

# %%
print ((aa - 25) % 50.) - (aa-25)

# %%
print ((bb + 25) % 50.) - (bb+25)

# %%
bb % 50

# %%
bb - (bb% 50)

# %%
aa - (aa % 50)

# %%

# %%
