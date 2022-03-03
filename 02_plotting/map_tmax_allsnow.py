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
# # %matplotlib inline

# %%
from IPython import parallel
c = parallel.Client()
view = c.load_balanced_view()

# %%
c.ids

# %%
# #%%px --local

# %matplotlib inline

import matplotlib as mpl
#mpl.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap

#from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# import pylab


def chart_info(Series_data):
    """create a basemap from a shapefile and the information contained in a series coded by hru id
    
    Arguments
    ===========
    Series_data: a pandas Series
                This series should contains the data for each region with the same code 
                that can be found under the shapefile
                
    ## THESE ARE HARD-CODED BELOW
    shapefile: string
                The location of the shapefile that contains the information about the geography
    extension: 4-tuple
                the extension of the projection, in the format north, south, east, west
    shape_key: string
                The name of the field in the shapefile that indicate the region. 
                The values of this field should match those in Series_data
    title: string, optional
                Title of the plot
    cmap: pylab Colormap or list of colour names, optional
                This gives the colormap that will be used to colorize the plot
    
    Returns
    ===========
    ax: pylab.Axes
                The axes on which the plot has been drawn
    """
    #NAME OUTPUT FILES AND TITLE NAME
#     yr = Series_data.name.replace("Unnamed: ","")
    title = 'tmax_allrain'
#     outfile = '/cxfs/projects/usgs/water/mows/NHM/khakala/python_scripts/nrel_updated_' + str(yr) + '.png'
    
    # Name of shapefile
    shpfile='/media/scratch/PRMS/notebooks/nhru_10U/nhru_10U_simpl'
    #shpfile='/media/scratch/PRMS/GIS/all_nhru_nad/nhruNationalIdentifier'
    
    #CHANGE COLORMAP
    #cmap = ['Green', 'W','Sienna']
    #cmap = 'YlGnBu_r'
    cmap = 'seismic'
    ##################
    
    #extent=(50, 22, -64, -119) # Extent for USA
    extent = (50, 42, -95, -114)  # Extent for r10U
    shape_key='hru_id_reg'
    
    #Start function
    fig, ax = plt.subplots(1,figsize=(20,30))
    ax = plt.gca()
    
    # create the colormap if a list of names is given, otherwise
    # use the given colormap
    lscm = mpl.colors.LinearSegmentedColormap
    if isinstance(cmap,(list,tuple)):
        cmap = lscm.from_list('mycm', cmap)
    else:
        cmap = plt.get_cmap(cmap)
    
    # create a new basemap with the given extent
    # TODO: allow to create more general projections
    print "Loading basemap..."
    north, south, east, west = extent
    m = Basemap(llcrnrlon=west,llcrnrlat=south,urcrnrlon=east,urcrnrlat=north, resolution='c',
            projection='laea', lat_0=(south+north)/2, lon_0=(east+west)/2)
    
    # use basemap the read and draw the shapefile
    # it will add two variables to the basemap, m.nhruDd and m.nhruDd_info
    print 'Loading shapefile...'
    m.readshapefile(shpfile,'nhruDd',drawbounds=False);
    
    # draw parallels.
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=20)
    
    # draw meridians
    meridians = np.arange(180.,360.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=20)
    m.drawmapboundary()
    
    # find minimum and maximum of the dataset to normalize the colors
    max_val = 40.
    min_val = 24.
#    max_val = Series_data.max()
#    min_val = Series_data.min()
    
    # cycle through states, color each one.
    # m.nhruDd contains the lines of the borders
    # m.nhruDd_info contains the info on the hru, like the name
    print 'Color HRUs...'
    for nhruDd_borders, nhruDd_info in zip(m.nhruDd, m.nhruDd_info):
        index = nhruDd_info[shape_key]
    
        #skip those that aren't in the dataset without complaints
        if index in Series_data:
            #set the color for each region
            val = Series_data.loc[index]
            
            if pd.isnull(val):
                # Record exists but the value is NaN
                color = '#ff0099'
            elif val > 50.:
                color = '#DB7093'
            elif val < -10.:
                color = '#8A2BE2'
            else:
                color = cmap((val - min_val) / (max_val - min_val))
        else:
            # The record is totally missing
            color = '#ff3300'
            
        #extract the x and y of the countours and plot them
        xx, yy = zip(*nhruDd_borders)
        patches = ax.fill(xx, yy, facecolor=color, edgecolor=color)
        
    ax.set_title(title);
    
    #generate a synthetic colorbar starting from the maximum and minimum of the dataset
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "2%", pad="3%")
    #axc, kw = mpl.colorbar.make_axes(ax)
    norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
    #norm = mpl.colors.Normalize(vmin=-2.5, vmax=2.5)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cb1.set_label('Temperature for all-rain event')
    cb1.ax.tick_params(labelsize=20)
    ax.patch.set_facecolor('0.93')
    #cb1.set_clim(-2.5,2.5)
#     fig.savefig('/media/scratch/PRMS/notebooks/crap.png', dpi=250, bbox_inches='tight')
#     plt.close(fig)
    #return ax


# %%
#Load .csv file to color shapefile
#data_AnnualTot = pd.read_csv('/cxfs/projects/usgs/water/mows/NHM/khakala/python_scripts/AllConus_nrel_updated.csv')
df = pd.read_csv('/media/scratch/PRMS/notebooks/02_test_and_prototype/hru_ranges_rain_v2.csv', index_col=0)
df.head()

# %%
# Select the month to plot on the basemap
data = df.iloc[:,2]

#print data.head(15)

#if 14 in data:
#    print data.loc[14]
    
min_val = data.min()
max_val = data.max()
print "min:", min_val
print "max:", max_val
#print 'max:', max_val
#print 'min:', min_val
#val = max_val
#print (val - min_val) / (max_val - min_val)
#
#pd.isnull(val)


# %%
chart_info(data)

# %%
# %%time
extent = (50, 42, -95, -114)  # Extent for r10U
north, south, east, west = extent
m = Basemap(llcrnrlon=west,llcrnrlat=south,urcrnrlon=east,urcrnrlat=north, resolution='c',
        projection='laea', lat_0=(south+north)/2, lon_0=(east+west)/2)

# %%
# %%time
shpfile='/media/scratch/PRMS/notebooks/nhru_10U/nhru_10U'
m.readshapefile(shpfile,'nhruDd',drawbounds=False);

# %%
# %%time
import shapefile

sf = shapefile.Reader('%s.shp' % shpfile)
recs    = sf.records()
shapes  = sf.shapes()
Nshp    = len(shapes)

# %%
# From: http://stackoverflow.com/questions/15968762/shapefile-and-matplotlib-plot-polygon-collection-of-shapefile-coordinates

from osgeo import ogr
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

# Extract first layer of features from shapefile using OGR
ds = ogr.Open('/media/scratch/PRMS/notebooks/nhru_10U/nhru_10U.shp')
nlay = ds.GetLayerCount()
lyr = ds.GetLayer(0)

# Get extent and calculate buffer size
ext = lyr.GetExtent()
xoff = (ext[1] - ext[0]) / 50
yoff = (ext[3] - ext[2]) / 50

# Prepare figure
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(ext[0]-xoff,ext[1]+xoff)
ax.set_ylim(ext[2]-yoff,ext[3]+yoff)

paths = []
lyr.ResetReading()

# Read all features in layer and store as paths
for feat in lyr:
    geom = feat.geometry()
    codes = []
    all_x = []
    all_y = []
    for i in range(geom.GetGeometryCount()):
        # Read ring geometry and create path
        r = geom.GetGeometryRef(i)
        x = [r.GetX(j) for j in range(r.GetPointCount())]
        y = [r.GetY(j) for j in range(r.GetPointCount())]
        # skip boundary between individual rings
        codes += [mpath.Path.MOVETO] + \
                     (len(x)-1)*[mpath.Path.LINETO]
        all_x += x
        all_y += y
    path = mpath.Path(np.column_stack((all_x,all_y)), codes)
    paths.append(path)

# %%

# %%
view.map(chart_info, data)

# %%
# From: http://stackoverflow.com/questions/15968762/shapefile-and-matplotlib-plot-polygon-collection-of-shapefile-coordinates

from osgeo import ogr
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

# Extract first layer of features from shapefile using OGR
ds = ogr.Open('/media/scratch/PRMS/notebooks/nhru_10U/nhru_10U.shp')
nlay = ds.GetLayerCount()
lyr = ds.GetLayer(0)

# Get extent and calculate buffer size
ext = lyr.GetExtent()
xoff = (ext[1] - ext[0]) / 50
yoff = (ext[3] - ext[2]) / 50

# Prepare figure
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(ext[0]-xoff,ext[1]+xoff)
ax.set_ylim(ext[2]-yoff,ext[3]+yoff)

paths = []
lyr.ResetReading()

# Read all features in layer and store as paths
for feat in lyr:
    geom = feat.geometry()
    codes = []
    all_x = []
    all_y = []
    for i in range(geom.GetGeometryCount()):
        # Read ring geometry and create path
        r = geom.GetGeometryRef(i)
        x = [r.GetX(j) for j in range(r.GetPointCount())]
        y = [r.GetY(j) for j in range(r.GetPointCount())]
        # skip boundary between individual rings
        codes += [mpath.Path.MOVETO] + \
                     (len(x)-1)*[mpath.Path.LINETO]
        all_x += x
        all_y += y
    path = mpath.Path(np.column_stack((all_x,all_y)), codes)
    paths.append(path)

# Add paths as patches to axes
for path in paths:
    patch = mpatches.PathPatch(path, \
            facecolor='blue', edgecolor='black')
    ax.add_patch(patch)

ax.set_aspect(1.0)
plt.show()

# %%
