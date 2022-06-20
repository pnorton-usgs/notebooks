# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %% [markdown]
# # Plot by-HRU NHM output variables

# %%
import cartopy.crs as ccrs
# from cartopy.io import shapereader
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import Normalize   # , LogNorm, PowerNorm
from matplotlib.collections import PatchCollection
import geopandas

from matplotlib.patches import Polygon
import shapely

import matplotlib as mpl

import pyproj as prj
from osgeo import ogr

# import xarray

# import netCDF4 as cf
import pandas as pd
# from datetime import datetime
# from calendar import monthrange
# from collections import namedtuple
import numpy as np
# import os


from pyPRMS.plot_helpers import get_projection

# %%
base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/src/code_lauren/NHM_CONUS_CALIBRATION/TandM_FIGURES/LAPTOP/T&Mfigures'
# filename = '{}/pkwater_equiv.nc'.format(work_dir)

gis_dir = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/NHM_v1.0'
shpfile = f'{gis_dir}/all_nhru_simple_aea_clipped/nhruNationalIdentifier.shp'
# shpfile_extent = f'{gis_dir}/all_nhru_simple_aea/nhruNationalIdentifier.shp'

states_file = f'{gis_dir}/states_outline_aea/states_outline.shp'

# %%
df = pd.read_csv(f'{base_dir}/DATA/RUNmean.csv', sep=' ', skipinitialspace=True)
df.index = df.index + 1
df.rename(columns={'x': 'RUN'}, inplace=True)
df.head()

# %%

# %%

# %%

# %%

# %%

# %%
# Read the shapefile
hru_df = geopandas.read_file(shpfile)

minx, miny, maxx, maxy = hru_df.geometry.total_bounds
# hru_df = hru_df.to_crs('epsg:4326')


crs_proj = get_projection(hru_df)

# if hru_df.crs.name == 'USA_Contiguous_Albers_Equal_Area_Conic_USGS_version':
#     print('Overriding USGS aea crs with EPSG:5070')
#     hru_df.crs = 'EPSG:5070'


# %%
hru_df.crs

# %%
states_df = geopandas.read_file(states_file)

# %%

# %%

# %%
min_val = 0.0
max_val = 200.0

# %%
# Create the colormap
# cmap = 'BrBG' #'GnBu_r' # for snow
# cmap = 'GnBu_r'
cmap = 'viridis'

# create the colormap if a list of names is given, otherwise
# use the given colormap
lscm = mpl.colors.LinearSegmentedColormap
if isinstance(cmap,(list,tuple)):
    cmap = lscm.from_list('mycm', cmap)
else:
    cmap = plt.get_cmap(cmap).copy()
    
missing_color = '#ff00cb'   # pink/magenta 

# norm = PowerNorm(gamma=0.05)
norm = Normalize(vmin=min_val, vmax=max_val)


# bnds = np.arange(min_val, max_val+2) - 0.5
# num_col = abs(max_val - min_val) + 1
# norm = colors.BoundaryNorm(boundaries=bnds, ncolors=num_col)


# %%
def plot_polygon_collection(ax, geoms, values=None, cmap='Set1', norm=None, facecolor=None, edgecolor='none',
                            alpha=1.0, linewidth=1.0, **kwargs):
    """ Plot a collection of Polygon geometries """
    # from https://stackoverflow.com/questions/33714050/geopandas-plotting-any-way-to-speed-things-up
    patches = []

    for poly in geoms:
        a = np.asarray(poly.exterior)
        if poly.has_z:
            poly = shapely.geometry.Polygon(zip(*poly.exterior.xy))

        patches.append(Polygon(a))

    patches = PatchCollection(patches, match_original=False,  
                              edgecolor='face', linewidth=linewidth, alpha=alpha, cmap=cmap, norm=norm, **kwargs)

    if values is not None:
        patches.set_array(values)
        # patches.set_cmap(cmap)

    ax.add_collection(patches, autolim=True)
    ax.autoscale_view()
    return patches


# %%
plt_vars = {'RUNmean': {'varname': 'RUN',
                        'file': 'DATA/RUNmean.csv', 
                        'title': 'Mean RUN',
                        'ticks': [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200],
                        'over': [200.0, 'yellow']},
            'AETmean': {'varname': 'AET',
                        'file': 'DATA/AETmean.csv',
                        'title': 'Mean AET',
                        'ticks': [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 , 0.11, 0.12, 0.13],
                        'under': [0.0, 'magenta']},
            'SCAmean': {'varname': 'SCA',
                        'file': 'DATA/SCAmean.csv',
                        'title': 'Mean SCA',
                        'ticks': [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
                        'under': [0.0, 'gray']},
            'RCHmean': {'varname': 'RCH',
                        'file': 'DATA/RCHmean.csv',
                        'title': 'Mean RCH',
                        'ticks': [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
                        'under': [0.0, 'magenta']},
            'SOMmean': {'varname': 'SOM',
                        'file': 'DATA/SOMmean.csv',
                        'title': 'Mean SOM',
                        'ticks': [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65],
                        'under': [0.0, 'magenta']}}


# %%
# %%time 
# platecarree = ccrs.PlateCarree(globe=ccrs.Globe(datum='NAD83'))

print('Build dataset')
# This takes care of multipolygons that are in the NHM geodatabase/shapefile
geoms_exploded = hru_df.explode(index_parts=True).reset_index(level=1, drop=True)

# fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(8, 24),
#                          subplot_kw=dict(projection=platecarree, frameon=False))
fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(8, 24),
                         subplot_kw=dict(projection=crs_proj, frameon=False))
ax = axes.flatten()
cc = 0

fig.suptitle('(A) Mean', fontsize=20)

for cvar, cdetail in plt_vars.items():
    print(f'{cc}: {cvar}')
    # Read the data
    df = pd.read_csv(f'{base_dir}/{cdetail["file"]}', sep=' ', skipinitialspace=True)
    df.index = df.index + 1
    df.rename(columns={'x': cdetail['varname']}, inplace=True)
    
    # Merge data with shapefile
    df_mrg = geoms_exploded.merge(df, left_on='hru_id_nat', right_index=True, how='left')

    # create the colormap if a list of names is given, otherwise
    # use the given colormap
    cmap = 'viridis'
    lscm = mpl.colors.LinearSegmentedColormap
    if isinstance(cmap,(list,tuple)):
        cmap = lscm.from_list('mycm', cmap)
    else:
        cmap = plt.get_cmap(cmap).copy()

    missing_color = '#ff00cb'   # pink/magenta 
    norm = Normalize(vmin=min(cdetail['ticks']), vmax=max(cdetail['ticks']))    
    
    print('Setup plot area')
    # ax[cc].coastlines()
    # ax[cc].gridlines()
    ax[cc].set_extent([minx, maxx, miny, maxy], crs=crs_proj)
    # ax[cc].set_extent([minx, maxx, miny, maxy], crs=platecarree)

    # ax[cc].outline_patch.set_linewidth(0)
    
    if 'over' in cdetail and 'under' in cdetail:
        cextend = 'both'
        cmap.set_over(cdetail['over'][1])
        cmap.set_under(cdetail['under'][1])
    elif 'over' in cdetail:
        cextend = 'max'
        cmap.set_over(cdetail['over'][1])
    elif 'under' in cdetail:
        cextend = 'min'
        cmap.set_under(cdetail['under'][1])
    else:
        cextend = 'neither'
        
    norm = mpl.colors.BoundaryNorm(cdetail['ticks'], cmap.N, extend=cextend)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapper.set_array(df_mrg[cdetail['varname']])

    cb = plt.colorbar(mappable=mapper, ax=ax[cc], shrink=0.8, ticks=cdetail['ticks'])
    cb.ax.tick_params(labelsize=10)
    
    ax[cc].set_title(cdetail['title'], fontsize=14)

    # print('Call plot_polygon_collection()')
    col = plot_polygon_collection(ax=ax[cc], geoms=df_mrg.geometry, 
                                  values=df_mrg[cdetail['varname']], **dict(cmap=cmap, norm=norm))

    states_df.plot(ax=ax[cc], facecolor='none', edgecolor='black')
    cc += 1

# plt.show()
plt.savefig(f'/Users/pnorton/tmp/the_mean.pdf', dpi=150, bbox_inches='tight')

# Close the figure so we don't chew up memory
fig.clf()
plt.close()

# %% [markdown]
# ### Single map

# %%
# %%time
cdetail = plt_vars['RUNmean']

df = pd.read_csv(f'{base_dir}/{cdetail["file"]}', sep=' ', skipinitialspace=True)
df.index = df.index + 1
df.rename(columns={'x': cdetail['varname']}, inplace=True)

# Merge data with shapefile
df_mrg = geoms_exploded.merge(df, left_on='hru_id_nat', right_index=True, how='left')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5),
                         subplot_kw=dict(projection=crs_proj, frameon=False))

fig.suptitle('(A) Mean', fontsize=20)

# create the colormap if a list of names is given, otherwise
# use the given colormap
cmap = 'viridis'
lscm = mpl.colors.LinearSegmentedColormap
if isinstance(cmap,(list,tuple)):
    cmap = lscm.from_list('mycm', cmap)
else:
    cmap = plt.get_cmap(cmap).copy()

missing_color = '#ff00cb'   # pink/magenta 
norm = Normalize(vmin=min(cdetail['ticks']), vmax=max(cdetail['ticks']))    

print('Setup plot area')
ax.set_extent([minx, maxx, miny, maxy], crs=crs_proj)
# ax[cc].set_extent([minx, maxx, miny, maxy], crs=platecarree)

if 'over' in cdetail and 'under' in cdetail:
    cextend = 'both'
    cmap.set_over(cdetail['over'][1])
    cmap.set_under(cdetail['under'][1])
elif 'over' in cdetail:
    cextend = 'max'
    cmap.set_over(cdetail['over'][1])
elif 'under' in cdetail:
    cextend = 'min'
    cmap.set_under(cdetail['under'][1])
else:
    cextend = 'neither'

norm = mpl.colors.BoundaryNorm(cdetail['ticks'], cmap.N, extend=cextend)
mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
mapper.set_array(df_mrg[cdetail['varname']])

cb = plt.colorbar(mappable=mapper, ax=ax, shrink=0.8, ticks=cdetail['ticks'])
cb.ax.tick_params(labelsize=10)

ax.set_title(cdetail['title'], fontsize=14)
# fig.title(ax[cc], cdetail['title'])

# print('Call plot_polygon_collection()')
col = plot_polygon_collection(ax=ax, geoms=df_mrg.geometry, 
                              values=df_mrg[cdetail['varname']], **dict(cmap=cmap, norm=norm))

states_df.plot(ax=ax, facecolor='none', edgecolor='black')

# %%

# %%

# %%

# %%

# %%

# %%

# %%
