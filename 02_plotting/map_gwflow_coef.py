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
# from IPython import parallel
# c = parallel.Client()
# view = c.load_balanced_view()

# %%
# c.ids

# %%
# # %%px --local

# %matplotlib inline
import matplotlib as mpl
#from matplotlib.backends.backend_pdf import PdfPages
#mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import pylab
from mpl_toolkits.basemap import Basemap

def chart_info(Series_data):
    """create a basemap from a shapefile and the information contained in a series coded by hru id
    
    Arguments
    ===========
    Series_data: a pandas Series
                This series should contains the data for each region with the same code that can be found under the shapefile
    ## THIS ARE HARD-CODED BELOW
    shapefile: string
                The location of the shapefile that contains the information about the geography
    extension: 4-tuple
                the extension of the projection, in the format north, south, east, west
    shape_key: string
                tha name of the field in the shape file that indicate the region. The values of this field
                should match those on the Series_data
    title: string, optional
                title of the plot
    cmap: pylab colormap or list of colour names, optional
                this gives the colormap that will be used to colorize the plot
    
    Returns
    ===========
    ax: pylab.Axes
                the axes on which the plot has been drawn
    """
    #NAME OUTPUT FILES AND TITLE NAME
#     yr = Series_data.name.replace("Unnamed: ","")
    title = 'gwflow_coef'
#     outfile = '/cxfs/projects/usgs/water/mows/NHM/khakala/python_scripts/nrel_updated_' + str(yr) + '.png'
    
    #Name all inputs:
    #shapefile='/cxfs/projects/usgs/water/mows/NHM/nhmHruShapefiles/hrusAllConus/hrusAllConusDd'
    shapefile='/Users/pnorton/Projects/National_Hydrology_Model/GIS/tmp/all_hru/nhruNationalIdentifier'
    
    #CHANGE COLORMAP
    #cmap = ['w','b']
    #cmap = ['Green', 'W','Sienna']
    cmap = ['W','Sienna']
    ##################
    
    extension=(50, 22, -64, -119) #Extent for USA
    shape_key='hru_id_nat'
    
    #Start function
    fig, ax = pylab.subplots(1,figsize=(20,30))
    ax = pylab.gca()
    
    # create the colormap if a list of names is given, otherwise
    # use the given colormap
    lscm = mpl.colors.LinearSegmentedColormap
    if isinstance(cmap,(list,tuple)):
        cmap = lscm.from_list('mycm',cmap)
    
    #create a new basemap with the given extension
    #TODO: allow to create more general projections
    north, south, east, west = extension
    m = Basemap(llcrnrlon=west,llcrnrlat=south,urcrnrlon=east,urcrnrlat=north,
            projection='laea', lat_0=(south+north)/2, lon_0=(east+west)/2)
    
    #use basemap the read and draw the shapefile
    #it will add two variables to the basemap, m.nhruDd and m.nhruDd_info
    m.readshapefile(shapefile,'nhruDd',drawbounds=False);
    
    # draw parallels.
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=20)
    
    # draw meridians
    meridians = np.arange(180.,360.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=20)
    m.drawmapboundary()
    
    #find minimum and maximum of the dataset to normalize the colors
    #max_pop = Series_data.max()*1.0
    #min_pop = Series_data.min()*1.0
    #max_pop = 0.6 for updated PET
    max_pop = 1000
    min_pop = 0
    
    # cycle through states, color each one.
    # m.states contains the lines of the borders
    # m.states_info contains the info on the region, like the name
#     for nhruDd_borders, nhruDd_info in zip(m.nhruDd, m.nhruDd_info):
#         index = nhruDd_info[shape_key]
    
#         #skip those that aren't in the dataset without complaints
#         if index not in Series_data:
#             continue
            
#         #set the color for each region
#         pop = Series_data[index-1]# Need index-1 because panda starts rows at 0 not 1
#         color =  cmap((pop-min_pop) / (max_pop-min_pop))
        
#         #extract the x and y of the countours and plot them
#         xx,yy = zip(*nhruDd_borders)
#         patches = ax.fill(xx,yy,facecolor=color,edgecolor=color)
        
    ax.set_title(title);
    
    #generate a synthetic colorbar starting from the maximum and minimum of the dataset
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "2%", pad="3%")
    #axc, kw = mpl.colorbar.make_axes(ax)
    norm = mpl.colors.Normalize(vmin=min_pop, vmax=max_pop)
    #norm = mpl.colors.Normalize(vmin=-2.5, vmax=2.5)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cb1.set_label('percent')
    cb1.ax.tick_params(labelsize=20)
    ax.patch.set_facecolor('0.93')
    #cb1.set_clim(-2.5,2.5)
#     fig.savefig(outfile, dpi=250, bbox_inches='tight')
    #pyplot.close(fig)
    #return ax


# %%
#Load .csv file to color shapefile
data_AnnualTot = pd.read_csv('/cxfs/projects/usgs/water/mows/NHM/khakala/python_scripts/AllConus_nrel_updated.csv')

# %%
#data_AnnualTot

# %%
time = ['Unnamed: 1', 'Unnamed: 2', 'Unnamed: 3', 'Unnamed: 4', 'Unnamed: 5', 'Unnamed: 6', 'Unnamed: 7', 'Unnamed: 8', 'Unnamed: 9', 'Unnamed: 10', 'Unnamed: 11', 'Unnamed: 12']
# data = [data_AnnualTot[time[0]], data_AnnualTot[time[1]]]
data = [data_AnnualTot[time[0]], data_AnnualTot[time[1]], data_AnnualTot[time[2]], data_AnnualTot[time[3]], data_AnnualTot[time[4]], data_AnnualTot[time[5]], data_AnnualTot[time[6]], data_AnnualTot[time[7]], data_AnnualTot[time[8]], data_AnnualTot[time[9]], data_AnnualTot[time[10]], data_AnnualTot[time[11]]]

# %%
chart_info(None)

# %%
# data

# %%
view.map(chart_info, data)

# %%
# type(data)

# %%
# chart_info(data)

# %%
#data_AnnualTot

# %%
#max_test = data_AnnualTot.max()*1.0

# %%
#min_test = data_AnnualTot.min()*1.0

# %%
#max_test

# %%
#min_test

# %%
