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
# %matplotlib inline

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
import pandas as pd
import os
from osgeo import ogr
import sys
# import pickle

import pyPRMS.ParameterFile as pf
reload(pf)

# %%
# Create the colormap
cmap = 'GnBu_r' # for snow
# cmap = 'YlGnBu_r'
#cmap = 'OrRd'  # for liquid
#cmap = 'seismic'
#cmap = ['Green', 'W','Sienna']

# Create the colormap if a list of names is given, otherwise use the given colormap
lscm = mpl.colors.LinearSegmentedColormap
if isinstance(cmap,(list,tuple)):
    cmap = lscm.from_list('mycm', cmap)
else:
    cmap = plt.get_cmap(cmap)
    
missing_color = '#ff00cb'   # pink/magenta 

# %% [markdown]
# ### Read in parameter file and create dataframe of selected parameter from a PRMS parameter file

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
outfile = os.path.normpath('map_pfile_%s.png' % cblabel)
print(outfile)

filename = os.path.normpath('/Users/pnorton/Projects/National_Hydrology_Model/datasets/regions/subset_testing/red_river.param')
paramset = pf.ParameterFile(filename)

df = paramset.parameterset.parameters.get_DataFrame(cblabel)
# df = paramset.parameterset.get_DataFrame(cblabel)
print(df.head())

# Set the min and max values allowed - right now we just take the min and max of the dataframe values
min_val = df.min().min()
max_val = df.max().max()

# %%
df.shape

# %%

# %% [markdown]
# ### Get extent information from the national HRUs shapefile

# %%
# shpfile = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/regions/subset_testing/HRU_subset_NAD/HRU_subset.shp'
shpfile = os.path.normpath('/Users/pnorton/Projects/National_Hydrology_Model/datasets/regions/subset_testing/red_river/HRU_subset_NAD/HRU_subset.shp')

pfile = 'HRU_subset.pickle'

# Name of attribute to use. Change to match the name of the HRU id attribute in the shapefile
shape_key='hru_id_nat'

# Use gdal/ogr to get the extent information
# Shapefile must must unprojected and supply lat/lon values
inDriver = ogr.GetDriverByName("ESRI Shapefile")
inDataSource = inDriver.Open(shpfile, 0)
inLayer = inDataSource.GetLayer()
extent = inLayer.GetExtent()

west, east, south, north = extent
pad = 0.3    # amount to pad the extent with
east += pad
west -= pad
south -= pad
north += pad
print('\tExtent: ({0:f}, {1:f}, {2:f}, {3:f})'.format(west, east, south, north))

# %% [markdown]
# ### Create the map figure

# %%
if df.shape[1] > 1:
    print('Currently unable to handle 2D parameters')
else:
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(30,20))
    # ax = axes.flatten()
    ax = axes

    # if os.path.isfile(pfile):
    #     print('Using pickled map')
    #     m = pickle.load(open(pfile, 'rb'))
    # else:
    # Load the basemap
    m = Basemap(llcrnrlon=west,llcrnrlat=south,urcrnrlon=east,urcrnrlat=north, resolution='c',
                projection='laea', lat_0=(south+north)/2, lon_0=(east+west)/2, ax=ax)
    #     pickle.dump(m, open(pfile,'wb'), -1)  # pickle it

    # Draw parallels
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
    pc.set_facecolor(cmap(norm(df_poly[cblabel].fillna(-99).values)))
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
print(df_poly.head())

# %%
# %timeit value = df.iloc[5]
# df

# %%
# %timeit value = df_dict[56254]
# df_dict[56254]

# %%
import os

thepath = '/Users/pnorton/'

print(os.path.normpath(thepath))

# %%
df.head()

# %%
df.rename(columns=['hru', 'paramval'], inplace=True)
df.head()

# %%
df_a = pd.DataFrame(df)
print(df_a.head())
df_a.reset_index(inplace=True)
print(df_a.head())
df_a.rename(columns={'index': 'hru', 0: 'paramval'}, inplace=True)
print(df_a.head())

# %%

# %%
