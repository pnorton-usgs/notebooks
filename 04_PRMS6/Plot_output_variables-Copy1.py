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
# %matplotlib inline

from __future__ import (absolute_import, division, print_function)
from future.utils import iteritems

import netCDF4 as cf
import pandas as pd
from datetime import datetime
from calendar import monthrange
from collections import namedtuple

import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize, LogNorm, PowerNorm

# mpl.use('Agg')
# from matplotlib.backends.backend_pdf import PdfPages

from mpl_toolkits.axes_grid1 import make_axes_locatable

# see https://github.com/matplotlib/basemap/issues/354
from mpl_toolkits.basemap import Basemap


import matplotlib.pyplot as plt

import numpy as np
import os
import pyproj as prj
from osgeo import ogr
import sys

# %%
# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/src/fortran/prms6/tests/red_river_of_the_south/output'
# filename = '{}/nhru_out_meanmonthly.nc'.format(work_dir)

# work_dir = '/Users/pnorton/tmp'
# filename = '{}/BCCA_0-125deg_tasmin_day_ACCESS1-0_historical_r1i1p1.nc'.format(work_dir)

# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190425_red_river'
# filename = '{}/daymet_v3_cbh.nc'.format(work_dir)

# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v3_1980-01-01_2016-12-31/tmp'
# filename = '{}/daymet_v3_cbh_tmin_20150101-20151231.nc'.format(work_dir)

# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/src/fortran/converters'
# filename = '{}/tmaxf.nc'.format(work_dir)

# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/NHM_output/netcdf'
# filename = '{}/dprst_stor_hru.nc'.format(work_dir)

work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/conus/output'
filename = '{}/summary_daily.nc'.format(work_dir)

# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks'
# filename = '{}/test1.nc'.format(work_dir)

first_time = True

# %%
fhdl = cf.Dataset(filename, 'r')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the size of the unlimited record dimension
dimensions = fhdl.dimensions
recdims = namedtuple('recdim', 'name, size')

print("======= Dimensions =======")
for xx, yy in dimensions.items():
    if yy.isunlimited():
        print('%s: %d (unlimited)' % (xx, len(yy)))
        recdim = recdims(name=xx, size=len(yy))
    else:
        if xx in ['time']:
            # Found a limited record dimension
            print('Fixed dimension, {}'.format(xx))
            recdim = recdims(name=xx, size=len(yy))
        print('%s: %d' % (xx, len(yy)))

# %%
# print("======= Variables =======")
# for xx in fhdl.variables:
#     print(xx)

# %%
d0units = fhdl.variables[recdim.name].units
print(d0units)

try:
    d0calendar = fhdl.variables[recdim.name].calendar
    print(d0calendar)
except KeyError:
    print('calendar attribute missing, assuming standard calendar')
    d0calendar = 'standard'

# %%
cblabel = 'pkwater_equiv'
time_idx = 100

# Read the time values into a list
timelist = cf.num2date(fhdl.variables[recdim.name][:], units=d0units, calendar=d0calendar)
# print(timelist)

# Create dataframe of the output variable
data_var = fhdl.variables[cblabel]

# Set the min and max values allowed - right now we just take the min and max of the dataframe values
min_val = np.min(data_var)
max_val = np.max(data_var)
val_rng = max_val - min_val


# var_df = pd.DataFrame(data_var[time_idx, :], columns=[cblabel])

# # Create dataframe of nhm_id
# # nhm_id_var = fhdl.variables['hru']
# nhm_id_var = fhdl.variables['nhm_id']
# # nhm_id_var = fhdl.variables['hru_feature_id']
# nhm_df = pd.DataFrame(nhm_id_var[:], columns=['nhru'])

# # Create a DataFrame of the output variable
# # df = var_df.merge(nhm_df, left_index=True, right_index=True)
# # df.set_index('nhru', inplace=True)
# df = var_df.join(nhm_df).set_index('nhru').sort_index()
# print(df.head())
# print(df.info())

# %%
timelist[time_idx].isoformat()

# %%
# Create the colormap
# cmap = 'BrBG' #'GnBu_r' # for snow
# cmap = 'GnBu_r'
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


# # Set the min and max values allowed - right now we just take the min and max of the dataframe values
# min_val = df.min().min()
# # max_val = 10.
# max_val = df.max().max()
# val_rng = max_val - min_val

# ### Get extent information from the national HRUs shapefile

# Need two shapefiles 1) in projected coordinates, 2) in geographic coordinates
# If gdal is installed can create geographic coordinates from projected with:
#   ogr2ogr -t_srs epsg:4326 output_wgs84.shp input.shp

# shpfile = '/Users/pnorton/Projects/National_Hydrology_Model/extraction_requests/20180307_red_river/GIS/HRU_subset_nad83.shp'
# shpfile_extent = '/Users/pnorton/Projects/National_Hydrology_Model/extraction_requests/20180307_red_river/GIS/HRU_subset_usaea.shp'

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



# %%
time_idx = 75

var_df = pd.DataFrame(data_var[time_idx, :], columns=[cblabel])

# Create dataframe of nhm_id
# nhm_id_var = fhdl.variables['hru']
nhm_id_var = fhdl.variables['nhm_id']
# nhm_id_var = fhdl.variables['hru_feature_id']
nhm_df = pd.DataFrame(nhm_id_var[:], columns=['nhru'])

# Create a DataFrame of the output variable
# df = var_df.merge(nhm_df, left_index=True, right_index=True)
# df.set_index('nhru', inplace=True)
df = var_df.join(nhm_df).set_index('nhru').sort_index()
print(df.head())
print(df.info())

# %%

# %%
# first_time = True

if df.shape[1] > 1:
    print('Currently unable to handle 2D parameters')
else:
    if first_time:
        # Create the map figure
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(30,20))
        # ax = axes.flatten()
        ax = axes
#         ax = plt.gca()

        # Load the basemap
        m = Basemap(width=east-west, height=north-south, resolution='c',
                    projection='laea', lat_0=cen_lat, lon_0=cen_lon, ax=ax)

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

        print('Create dataframe from shapefile')
        df_poly = pd.DataFrame({'shapes': [Polygon(np.array(ss), closed=True) for ss in m.nhruDd],
                                'id': [aa[shape_key] for aa in m.nhruDd_info]})


    print('Create plotting dataframe')
    df_plot = df_poly.merge(df, left_on='id', right_index=True, how='left')
    
    print('Patch Collection')
    pc = PatchCollection(df_plot.shapes, cmap=cmap, match_original=True, zorder=2)
    norm = PowerNorm(gamma=0.05)

#     if first_time:
#         print('mapping...')
#         mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        
#     mapper.set_array(df_plot[cblabel])
#     plt.colorbar(mapper, shrink=0.4)

    if first_time:
        ax.add_collection(pc)
        
        
    if first_time:
        print('facecolor')
    #     pc.set_array(df_plot[cblabel].fillna(-99).values)
        pc.set_facecolor(cmap(norm(df_plot[cblabel].fillna(-99).values)))        
    else:
        ax.collections[0].set_facecolor(cmap(norm(df_plot[cblabel].fillna(-99).values)))
        
    plt.title('Variable: {},  Date: {}'.format(cblabel, timelist[time_idx].isoformat()))
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    
    first_time = False
#     plt.show()

# %%
if df.shape[1] > 1:
    print('Currently unable to handle 2D parameters')
else:
    if first_time:
        # Create the map figure
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(30,20))
        # ax = axes.flatten()
        ax = axes
#         ax = plt.gca()

        # Load the basemap
        m = Basemap(width=east-west, height=north-south, resolution='c',
                    projection='laea', lat_0=cen_lat, lon_0=cen_lon, ax=ax)

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

        print('Create dataframe from shapefile')
        df_poly = pd.DataFrame({'shapes': [Polygon(np.array(ss), closed=True) for ss in m.nhruDd],
                                'id': [aa[shape_key] for aa in m.nhruDd_info]})
    
        print('Patch Collection')
        pc = PatchCollection(df_poly.shapes, cmap=cmap, match_original=True, zorder=2)
        norm = PowerNorm(gamma=0.05)

    print('Create plotting dataframe')
    df_plot = df_poly.merge(df, left_on='id', right_index=True, how='left')

    if first_time:
        print('mapping...')
        mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        
    mapper.set_array(df_plot[cblabel])
    plt.colorbar(mapper, shrink=0.4)

    if first_time:
        ax.add_collection(pc)
        
        
    if first_time:
        print('facecolor')
    #     pc.set_array(df_plot[cblabel].fillna(-99).values)
        pc.set_facecolor(cmap(norm(df_plot[cblabel].fillna(-99).values)))        
    else:
        ax.collections[0].set_facecolor(cmap(norm(df_plot[cblabel].fillna(-99).values)))
        
    if first_time:
        fig.canvas.draw()
        
    plt.title('Variable: {},  Date: {}'.format(cblabel, timelist[time_idx].isoformat()))
    # Output a png
    # fig.savefig(outfile, dpi=250, bbox_inches='tight',debug=True)
    # fig.clf()
    # plt.close()

    
    
    first_time = False

# %%
time_idx = 0

var_df = pd.DataFrame(data_var[time_idx, :], columns=[cblabel])

# Create dataframe of nhm_id
# nhm_id_var = fhdl.variables['hru']
nhm_id_var = fhdl.variables['nhm_id']
# nhm_id_var = fhdl.variables['hru_feature_id']
nhm_df = pd.DataFrame(nhm_id_var[:], columns=['nhru'])

# Create a DataFrame of the output variable
# df = var_df.merge(nhm_df, left_index=True, right_index=True)
# df.set_index('nhru', inplace=True)
df = var_df.join(nhm_df).set_index('nhru').sort_index()
print(df.head())
print(df.info())

df_plot = df_poly.merge(df, left_on='id', right_index=True, how='left')
ax.collections[0].set_facecolor(cmap(norm(df_plot[cblabel].fillna(-99).values)))

fig.canvas.draw()
fig.canvas.flush_events()

# %%
if df.shape[1] > 1:
    print('Currently unable to handle 2D parameters')
else:
    # Create the map figure
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
    norm = PowerNorm(gamma=0.05)
#     norm = LogNorm(vmin=min_val, vmax=max_val)
#     norm = Normalize(vmin=min_val, vmax=max_val)

    print('facecolor')
    pc.set_facecolor(cmap(norm(df_poly[cblabel].fillna(-99).values)))
    ax.add_collection(pc)

    print('mapping...')
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapper.set_array(df_poly[cblabel])
    plt.colorbar(mapper, shrink=0.4)

    plt.title('Variable: {},  Date: {}'.format(cblabel, timelist[time_idx].isoformat()))
    # Output a png
    # fig.savefig(outfile, dpi=250, bbox_inches='tight',debug=True)
    # fig.clf()
    # plt.close()

    plt.show()

# %%
cmap(norm(df_plot[cblabel].fillna(-99).values))[0]

# %%

# %%

# %%

# %%

# %%

# %%
min_val

# %%
max_val

# %%
df.max()

# %%
df.idxmax()

# %%
df.mean()

# %%
df.median()

# %%
df.std()

# %%
df.info()

# %%
np.min(fhdl.variables[cblabel])

# %%
cmap(norm(df_plot[cblabel].fillna(-99).values))[2]

# %%
