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
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%

from collections import OrderedDict
import geopandas as gpd
import numpy as np
import pandas as pd
# import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy
import cartopy.crs as ccrs
# from pyproj import CRS

# %%
def get_exceedance_curve(ts, fclass):
    stat_name = f'{ts.name}_{fclass}'
    
    ranked_stat = sorted(ts[ts.notnull()])
    prob = np.arange(len(ranked_stat), dtype=float) + 1.0
    prob = (prob / (len(ranked_stat) + 1.0))
    
    # Return dataframe of exceedence curve values
    return pd.DataFrame({'exceedance': prob, stat_name: ranked_stat}, columns=['exceedance', stat_name])


# %%
calibrations = {'PRECAL': {'workdir': '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHRU/20210617_gm_PRECAL'},
                'byHRU': {'workdir': '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHRU/20210624_calib'},
                'byHW': {'workdir': '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHW/20211027_gm_byHW'},
                'byHWobs': {'workdir': '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHWobs/20211129_gm_byHWobs'}}

# NHMv1.0
# calib_name = 'byHRU'
# calib_dir = '20210624_calib'
# workdir = f'/Volumes/USGS_NHM1/calibrations/NHMv10/DAYMET_releases/{calib_name}'

# NHMv1.1
# calib_name = 'byHW'
# calib_dir = '20211027_gm_byHW'
# workdir = f'/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHW/{calib_dir}'

regions = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/HUC2_fixed.gpkg'

# %%
cal_df = {}
cal_ref_df = {}
cal_nonref_df = {}
cal_nan_df = {}
plot_info = []
plot_vars = ['ns', 'nslog', 'mon_ns']
plot_titles = {'ns': 'Nash-Sutcliffe (daily)',
               'nslog': 'Nash-Sutcliffe (log-daily)',
               'mon_ns': 'Nash-Sutcliffe (monthly)'}
fig_titles = {'ref': 'Reference streamgage Nash-Sutcliffe values',
              'nonref': 'Non-reference streamgage Nash-Sutcliffe values',
              'nonfalcone': 'Non-Falcone streamgage Nash-Sutcliffe values'}
gage_class = 'nonref'  # One of: 'ref', 'nonref', 'nonfalcone'

fig_title = fig_titles[gage_class]

for cal, ci in calibrations.items():
    print(cal)
    cal_df[cal] = pd.read_csv(f'{ci["workdir"]}/gage_stats_{cal}.csv', sep=',')
    
    cal_ref_df[cal] = cal_df[cal][cal_df[cal]['falcone_class'] == 'Ref']
    cal_nonref_df[cal] = cal_df[cal][cal_df[cal]['falcone_class'] == 'Non-ref']
    cal_nan_df[cal] = cal_df[cal][cal_df[cal]['falcone_class'].isna()]
    
    for vv in plot_vars:
        if gage_class == 'ref':
            plot_info.append({'calname': cal, 'datasource': cal_ref_df[cal], 'plot_var': vv, 'plot_title': plot_titles[vv]})
        elif gage_class == 'nonref':
            plot_info.append({'calname': cal, 'datasource': cal_nonref_df[cal], 'plot_var': vv, 'plot_title': plot_titles[vv]})
        else:
            plot_info.append({'calname': cal, 'datasource': cal_nan_df[cal], 'plot_var': vv, 'plot_title': plot_titles[vv]})


# %%

# %%
# Get HUC2 region outlines
region_df = gpd.read_file(regions)

selected_regions = ['01', '02', '03', '04', '05', '06', '07', '08',
                    '09', '10', '11', '12', '13', '14', '15', '16', 
                    '17', '18']
map_regions = region_df[region_df['HUC2'].isin(selected_regions)]

# %% [markdown]
# ## Plot maps of streamgage NS

# %%
nrows = len(calibrations)
ncols = len(plot_vars)
marker_size = 10.0

p_obj = {}

crs_proj = ccrs.PlateCarree()

fig = plt.figure(figsize=(30, 13))
fig.suptitle(fig_title, fontsize=20)

minx, miny, maxx, maxy = map_regions.geometry.total_bounds

# Adjust the northernmost extent
maxy -= 3.0

cmap = mpl.colors.ListedColormap(['blue', 'cyan', 'lime', 'orange', 'red'])
cmap.set_under('black')

bounds = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

for ii, dd in enumerate(plot_info):
    dsrc = dd['datasource']
    pvar = dd['plot_var']
    
    p_obj[ii] = plt.subplot(nrows, ncols, ii+1, projection=crs_proj, frame_on=False)
    p_obj[ii].set_extent([minx, maxx, miny, maxy], ccrs.Geodetic())
    p_obj[ii].patch.set_visible(False)
    p_obj[ii].set_title(f'{dd["calname"]}, {dd["plot_title"]}')
    p_obj[ii].coastlines()
    
    map_regions.plot(ax=p_obj[ii], alpha=0.9, facecolor="none", edgecolor='dimgrey', transform=ccrs.Geodetic())
    
    p_obj[ii].scatter(x=dsrc['longitude'], y=dsrc['latitude'], s=marker_size, c=dsrc[pvar], 
                      alpha=0.6, edgecolor="none", cmap=cmap, norm=norm, transform=ccrs.Geodetic())

# plt.subplots_adjust(bottom=0.001)
plt.subplots_adjust(hspace=0.1)

#            (left, bottom, width, height)
cax = plt.axes([0.3, 0.1, 0.4, 0.008])

cb = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
                  cax=cax,
                  extend='min',
                  shrink=0.35,
                  ticks=bounds,
                  spacing='proportional',
                  orientation='horizontal',
                  label='Nash-Sutcliffe Efficiency')

cb.set_label(label='Nash-Sutcliffe Efficiency', size='large', weight='bold')
cb.ax.tick_params(labelsize='large')

plt.savefig(f'/Users/pnorton/tmp/new_ns_{gage_class}_map.png', dpi=150, bbox_inches='tight')

# Close the figure so we don't chew up memory
fig.clf()

# %%

# %%

# %% [markdown]
# ## Exceedence curves

# %%
# Layout info at: https://matplotlib.org/stable/tutorials/intermediate/constrainedlayout_guide.html
fig, ax = plt.subplots(nrows=len(plot_vars), ncols=1, figsize=(10, 10), constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0.1, wspace=0.2)

plot_labels = ['ref', 'nonref']

color_set = {'PRECAL': 'darkgrey',
             'byHRU': 'royalblue',
             'byHW': 'coral',
             'byHWobs': 'darkgreen'}

fal_cls = {'ref': 'Ref', 'nonref': 'Non-ref'}
line_style = {'ref': 'solid', 'nonref': 'dashed'}

for idx, vv in enumerate(plot_vars):
    ax[idx].set_title(plot_titles[vv])
    
for cal, ci in calibrations.items():
    for idx, vv in enumerate(plot_vars):
        for ridx, rr in enumerate(plot_labels):
            objfun = get_exceedance_curve(cal_df[cal][cal_df[cal]['falcone_class'] == fal_cls[rr]].loc[:, vv], rr)
            
            if idx==0:
                objfun.plot(ax=ax[idx], x='exceedance', y=objfun.columns[1], ylim=(0.0, 1.0), 
                            color=color_set[cal], linestyle=line_style[rr], label=f'{cal}-{rr}')
                ax[idx].xaxis.label.set_visible(False)
                ax[idx].yaxis.set_label_text('Nash-Sutcliffe Goodness of Fit')                
            else:
                objfun.plot(ax=ax[idx], x='exceedance', y=objfun.columns[1], ylim=(0.0, 1.0), 
                            color=color_set[cal], linestyle=line_style[rr], legend=False)
                
        if idx == 1:
            ax[idx].xaxis.label.set_visible(False)

        ax[idx].grid(color='lightblue', linestyle='dotted')

plt.savefig(f'/Users/pnorton/tmp/new_ns_curves.png', dpi=150, bbox_inches='tight')

# Close the figure so we don't chew up memory
fig.clf()

# %%

# %%

# %%

# %%

# %%

# %%
# df = pd.read_csv(f'{workdir}/gage_stats_{calib_name}.csv', sep=',')
# df.head()

# %%
# ref_df = df[df['falcone_class'] == 'Ref']
# nonref_df = df[df['falcone_class'] == 'Non-ref']
# nan_df = df[df['falcone_class'].isna()]
# ref_df.plot(x='longitude', y='latitude', kind='scatter', c='ns', figsize=(20,10), colormap='turbo', vmin=-1.0, vmax=1.0)

# %% [markdown]
# ## Read the PRECAL stats

# %%
# calib_name = 'PRECAL'
# calib_dir = '20210617_gm_PRECAL'
# workdir = f'/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHRU/{calib_dir}'

# df_precal = pd.read_csv(f'{workdir}/gage_stats_{calib_name}.csv', sep=',')

# precal_ref_df = df_precal[df_precal['falcone_class'] == 'Ref']
# precal_nonref_df = df_precal[df_precal['falcone_class'] == 'Non-ref']
# precal_nan_df = df_precal[df_precal['falcone_class'].isna()]

# %%

# %%
region_df = gpd.read_file(regions)

selected_regions = ['01', '02', '03', '04', '05', '06', '07', '08',
                    '09', '10', '11', '12', '13', '14', '15', '16', 
                    '17', '18']
map_regions = region_df[region_df['HUC2'].isin(selected_regions)]

# %%


# map_regions.crs

# %%

# %%
# From: https://stackoverflow.com/questions/25428512/draw-a-map-of-a-specific-country-with-cartopy

# from cartopy.io import shapereader
# import numpy as np
# import geopandas
# import matplotlib.pyplot as plt

# import cartopy.crs as ccrs

# # get natural earth data (http://www.naturalearthdata.com/)

# # get country borders
# resolution = '10m'
# category = 'cultural'
# name = 'admin_0_countries'

# shpfilename = shapereader.natural_earth(resolution, category, name)

# # read the shapefile using geopandas
# df = geopandas.read_file(shpfilename)

# # read the german borders
# poly = df.loc[df['ADMIN'] == 'Germany']['geometry'].values[0]

# ax = plt.axes(projection=ccrs.PlateCarree())

# ax.add_geometries(poly, crs=ccrs.PlateCarree(), facecolor='none', 
#                   edgecolor='0.5')

# ax.set_extent([5, 16, 46.5, 56], crs=ccrs.PlateCarree())

# %% [markdown]
# ## Get Natural Earth countries shapefile

# %%
# from cartopy.io import shapereader

# resolution = '10m'
# category = 'cultural'
# name = 'admin_0_countries'

# shpfilename = shapereader.natural_earth(resolution, category, name)

# # read the shapefile using geopandas
# df = geopandas.read_file(shpfilename)

# # read the North America borders
# poly = df.loc[df['ADMIN'] == 'Germany']['geometry'].values[0]

# %%
minx, miny, maxx, maxy

# %%

# %%

# %%
# crs_proj = ccrs.LambertConformal()
# crs_proj = ccrs.AlbersEqualArea()
crs_proj = ccrs.PlateCarree()

# fig = plt.figure(figsize=(30, 20))
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(50, 20)) # , constrained_layout=True)

ax[0] = fig.add_axes([0, 0, 1, 1], projection=crs_proj, frameon=False)
ax[0].patch.set_visible(False)
ax[0].set_title('Nash-Sutcliffe (daily)')

ax[1] = fig.add_axes([0, 0, 1, 1], projection=crs_proj, frameon=False)
ax[1].patch.set_visible(False)
ax[1].set_title('Nash-Sutcliffe (log-daily)')

ax[2] = fig.add_axes([0, 0, 1, 1], projection=crs_proj, frameon=False)
ax[2].patch.set_visible(False)
ax[2].set_title('Nash-Sutcliffe (monthly)')

# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(30, 20))
# ax = fig.axes(projection=crs_proj)

minx, miny, maxx, maxy = map_regions.geometry.total_bounds
# ax.set_extent([-120, -66.5, 20, 50], ccrs.Geodetic())
ax[0].set_extent([minx, maxx, miny, maxy], ccrs.Geodetic())
ax[1].set_extent([minx, maxx, miny, maxy], ccrs.Geodetic())
ax[2].set_extent([minx, maxx, miny, maxy], ccrs.Geodetic())

# ax = plt.axes(projection=crs_proj)
ax[0].coastlines()
ax[1].coastlines()
ax[2].coastlines()
# ax.gridlines()

map_regions.plot(ax=ax[0], alpha=0.9, facecolor="none", edgecolor='k', transform=ccrs.Geodetic())
map_regions.plot(ax=ax[1], alpha=0.9, facecolor="none", edgecolor='k', transform=ccrs.Geodetic())
map_regions.plot(ax=ax[2], alpha=0.9, facecolor="none", edgecolor='k', transform=ccrs.Geodetic())
# norm = Normalize(vmin=-1.0, vmax=1.0)

# ref_df.plot(ax=ax, x='longitude', y='latitude', kind='scatter', c='ns', colormap='turbo', norm=norm, transform=ccrs.Geodetic())

cmap = mpl.colors.ListedColormap(['blue', 'cyan', 'lime', 'orange', 'red'])
cmap.set_under('black')

bounds = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)


# scatter data based on coords of centroid
# sc1 = ax[0].scatter(x=ref_df['longitude'], y=ref_df['latitude'], s=50, c=ref_df['ns'], alpha=0.6, cmap=cmap, norm=norm, transform=ccrs.Geodetic())
# sc2 = ax[1].scatter(x=ref_df['longitude'], y=ref_df['latitude'], s=50, c=ref_df['nslog'], alpha=0.6, cmap=cmap, norm=norm, transform=ccrs.Geodetic())
# sc3 = ax[2].scatter(x=ref_df['longitude'], y=ref_df['latitude'], s=50, c=ref_df['mon_ns'], alpha=0.6, cmap=cmap, norm=norm, transform=ccrs.Geodetic())

cb = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
#              cax=ax,
#              boundaries=[0] + bounds + [1],
                 extend='min',
                 shrink=0.35,
                 ticks=bounds,
                 spacing='proportional',
                 orientation='horizontal',
                 label='Nash-Sutcliffe Efficiency')


# cb = plt.colorbar(im, orientation="horizontal", pad=0.15)
cb.set_label(label='Nash-Sutcliffe Efficiency', size='large', weight='bold')
cb.ax.tick_params(labelsize='large')


# %% [markdown]
# ## Plot multiple maps in a single figure

# %%
nrows = 2
ncols = 3
marker_size = 10.0

gage_class = 'ref'  # One of: 'ref', 'nonref', 'nonfalcone'

if gage_class == 'ref':
    # Falcone reference gages
    fig_title = 'Reference streamgage Nash-Sutcliffe values'
    p_info = [{'calname': 'PRECAL', 'datasource': precal_ref_df, 'plot_var': 'ns', 'plot_title': 'Nash-Sutcliffe (daily)'},
              {'calname': 'PRECAL', 'datasource': precal_ref_df, 'plot_var': 'nslog', 'plot_title': 'Nash-Sutcliffe (log-daily)'},
              {'calname': 'PRECAL', 'datasource': precal_ref_df, 'plot_var': 'mon_ns', 'plot_title': 'Nash-Sutcliffe (monthly)'},
              {'calname': 'byHRU', 'datasource': ref_df, 'plot_var': 'ns', 'plot_title': 'Nash-Sutcliffe (daily)'},
              {'calname': 'byHRU', 'datasource': ref_df, 'plot_var': 'nslog', 'plot_title': 'Nash-Sutcliffe (log-daily)'},
              {'calname': 'byHRU', 'datasource': ref_df, 'plot_var': 'mon_ns', 'plot_title': 'Nash-Sutcliffe (monthly)'}]
elif gage_class == 'nonref':
    # Falcone non-reference gages
    fig_title = 'Non-reference streamgage Nash-Sutcliffe values'
    p_info = [{'calname': 'PRECAL', 'datasource': precal_nonref_df, 'plot_var': 'ns', 'plot_title': 'Nash-Sutcliffe (daily)'},
              {'calname': 'PRECAL', 'datasource': precal_nonref_df, 'plot_var': 'nslog', 'plot_title': 'Nash-Sutcliffe (log-daily)'},
              {'calname': 'PRECAL', 'datasource': precal_nonref_df, 'plot_var': 'mon_ns', 'plot_title': 'Nash-Sutcliffe (monthly)'},
              {'calname': 'byHRU', 'datasource': nonref_df, 'plot_var': 'ns', 'plot_title': 'Nash-Sutcliffe (daily)'},
              {'calname': 'byHRU', 'datasource': nonref_df, 'plot_var': 'nslog', 'plot_title': 'Nash-Sutcliffe (log-daily)'},
              {'calname': 'byHRU', 'datasource': nonref_df, 'plot_var': 'mon_ns', 'plot_title': 'Nash-Sutcliffe (monthly)'}]
else: 
    # Falcone non-reference gages
    fig_title = 'Non-Falcone streamgage Nash-Sutcliffe values'
    p_info = [{'calname': 'PRECAL', 'datasource': precal_nan_df, 'plot_var': 'ns', 'plot_title': 'Nash-Sutcliffe (daily)'},
              {'calname': 'PRECAL', 'datasource': precal_nan_df, 'plot_var': 'nslog', 'plot_title': 'Nash-Sutcliffe (log-daily)'},
              {'calname': 'PRECAL', 'datasource': precal_nan_df, 'plot_var': 'mon_ns', 'plot_title': 'Nash-Sutcliffe (monthly)'},
              {'calname': 'byHRU', 'datasource': nan_df, 'plot_var': 'ns', 'plot_title': 'Nash-Sutcliffe (daily)'},
              {'calname': 'byHRU', 'datasource': nan_df, 'plot_var': 'nslog', 'plot_title': 'Nash-Sutcliffe (log-daily)'},
              {'calname': 'byHRU', 'datasource': nan_df, 'plot_var': 'mon_ns', 'plot_title': 'Nash-Sutcliffe (monthly)'}]
    
p_obj = {}

crs_proj = ccrs.PlateCarree()

fig = plt.figure(figsize=(30, 13))
fig.suptitle(fig_title, fontsize=20)

minx, miny, maxx, maxy = map_regions.geometry.total_bounds

# Adjust the northernmost extent
maxy -= 3.0

cmap = mpl.colors.ListedColormap(['blue', 'cyan', 'lime', 'orange', 'red'])
cmap.set_under('black')

bounds = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

for ii, dd in enumerate(p_info):
    dsrc = dd['datasource']
    pvar = dd['plot_var']
    
    p_obj[ii] = plt.subplot(nrows, ncols, ii+1, projection=crs_proj, frame_on=False)
    p_obj[ii].set_extent([minx, maxx, miny, maxy], ccrs.Geodetic())
    p_obj[ii].patch.set_visible(False)
    p_obj[ii].set_title(f'{dd["calname"]}, {dd["plot_title"]}')
    p_obj[ii].coastlines()
    
    map_regions.plot(ax=p_obj[ii], alpha=0.9, facecolor="none", edgecolor='dimgrey', transform=ccrs.Geodetic())
    
    p_obj[ii].scatter(x=dsrc['longitude'], y=dsrc['latitude'], s=marker_size, c=dsrc[pvar], 
                      alpha=0.6, edgecolor="none", cmap=cmap, norm=norm, transform=ccrs.Geodetic())

# plt.subplots_adjust(bottom=0.001)
plt.subplots_adjust(hspace=0.1)

#            (left, bottom, width, height)
cax = plt.axes([0.3, 0.1, 0.4, 0.008])

cb = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
                  cax=cax,
                  extend='min',
                  shrink=0.35,
                  ticks=bounds,
                  spacing='proportional',
                  orientation='horizontal',
                  label='Nash-Sutcliffe Efficiency')

cb.set_label(label='Nash-Sutcliffe Efficiency', size='large', weight='bold')
cb.ax.tick_params(labelsize='large')

plt.savefig(f'/Users/pnorton/tmp/ns_{gage_class}_map.png', dpi=150, bbox_inches='tight')

# Close the figure so we don't chew up memory
fig.clf()

# %%

# %%

# %%

# %%

# %%
# byHRU stats
ns_ref_df = get_exceedance_curve(ref_df.loc[:, 'ns'], 'ref')
nslog_ref_df = get_exceedance_curve(ref_df.loc[:, 'nslog'], 'ref')
mon_ns_ref_df = get_exceedance_curve(ref_df.loc[:, 'mon_ns'], 'ref')

ns_nonref_df = get_exceedance_curve(nonref_df.loc[:, 'ns'], 'nonref')
nslog_nonref_df = get_exceedance_curve(nonref_df.loc[:, 'nslog'], 'nonref')
mon_ns_nonref_df = get_exceedance_curve(nonref_df.loc[:, 'mon_ns'], 'nonref')


# PRECAL stats
precal_ns_ref_df = get_exceedance_curve(df_precal[df_precal['falcone_class'] == 'Ref'].loc[:, 'ns'], 'ref')
precal_nslog_ref_df = get_exceedance_curve(df_precal[df_precal['falcone_class'] == 'Ref'].loc[:, 'nslog'], 'ref')
precal_mon_ns_ref_df = get_exceedance_curve(df_precal[df_precal['falcone_class'] == 'Ref'].loc[:, 'mon_ns'], 'ref')

precal_ns_nonref_df = get_exceedance_curve(df_precal[df_precal['falcone_class'] == 'Non-ref'].loc[:, 'ns'], 'nonref')
precal_nslog_nonref_df = get_exceedance_curve(df_precal[df_precal['falcone_class'] == 'Non-ref'].loc[:, 'nslog'], 'nonref')
precal_mon_ns_nonref_df = get_exceedance_curve(df_precal[df_precal['falcone_class'] == 'Non-ref'].loc[:, 'mon_ns'], 'nonref')

# %%

# %% [markdown]
# ## Plot exceedance curves of daily, log-daily, and monthly streamflow Nash-Sutcliff values

# %%
# Layout info at: https://matplotlib.org/stable/tutorials/intermediate/constrainedlayout_guide.html
fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10, 10), constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0.1, wspace=0.2)


calib_color = 'darkgrey'
# pre-calibration
ax[0].set_title('Nash-Sutcliffe (daily)')
precal_ns_ref_df.plot(ax=ax[0], x='exceedance', y=precal_ns_ref_df.columns[1], ylim=(0.0, 1.0), color=calib_color,
                      label='PRECAL-ref')

ax[1].set_title('Nash-Sutcliffe (log, daily)')
precal_nslog_ref_df.plot(ax=ax[1], x='exceedance', y=precal_nslog_ref_df.columns[1], ylim=(0.0, 1.0), color=calib_color,
                         legend=False)

ax[2].set_title('Nash-Sutcliffe (monthly)')
precal_mon_ns_ref_df.plot(ax=ax[2], x='exceedance', y=precal_mon_ns_ref_df.columns[1], ylim=(0.0, 1.0), color=calib_color,
                          legend=False)

ax[0].xaxis.label.set_visible(False)
ax[1].xaxis.label.set_visible(False)
ax[2].xaxis.label.set_visible(False)

ax[0].yaxis.set_label_text('Nash-Sutcliffe Goodness of Fit')
ax[1].yaxis.set_label_text('Nash-Sutcliffe Goodness of Fit')
ax[2].yaxis.set_label_text('Nash-Sutcliffe Goodness of Fit')


calib_color = 'royalblue'
# byHRU
ns_ref_df.plot(ax=ax[0], x='exceedance', y=ns_ref_df.columns[1], ylim=(0.0, 1.0), color=calib_color,
               label='byHRU-ref')
nslog_ref_df.plot(ax=ax[1], x='exceedance', y=nslog_ref_df.columns[1], ylim=(0.0, 1.0), color=calib_color,
                  legend=False)
mon_ns_ref_df.plot(ax=ax[2], x='exceedance', y=mon_ns_ref_df.columns[1], ylim=(0.0, 1.0), color=calib_color,
                   legend=False)

# Add in the non-ref gages

# pre-calibration
calib_color = 'darkgrey'
precal_ns_nonref_df.plot(ax=ax[0], x='exceedance', y=precal_ns_nonref_df.columns[1], ylim=(0.0, 1.0), color=calib_color, linestyle='dashed',
                         label='PRECAL-nonref')
precal_nslog_nonref_df.plot(ax=ax[1], x='exceedance', y=precal_nslog_nonref_df.columns[1], ylim=(0.0, 1.0), color=calib_color, linestyle='dashed',
                            legend=False)
precal_mon_ns_nonref_df.plot(ax=ax[2], x='exceedance', y=precal_mon_ns_nonref_df.columns[1], ylim=(0.0, 1.0), color=calib_color, linestyle='dashed',
                             legend=False)

# byHRU
calib_color = 'royalblue'
ns_nonref_df.plot(ax=ax[0], x='exceedance', y=ns_nonref_df.columns[1], ylim=(0.0, 1.0), color=calib_color, linestyle='dashed',
                  label='byHRU-nonref')
nslog_nonref_df.plot(ax=ax[1], x='exceedance', y=nslog_nonref_df.columns[1], ylim=(0.0, 1.0), color=calib_color, linestyle='dashed',
                     legend=False)
mon_ns_nonref_df.plot(ax=ax[2], x='exceedance', y=mon_ns_nonref_df.columns[1], ylim=(0.0, 1.0), color=calib_color, linestyle='dashed',
                      legend=False)


plt.savefig(f'/Users/pnorton/tmp/ns_curves.png', dpi=150, bbox_inches='tight')

# Close the figure so we don't chew up memory
fig.clf()

# %%

# %%

# %%
