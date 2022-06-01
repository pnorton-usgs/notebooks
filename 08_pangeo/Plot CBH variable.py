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
#     display_name: Python [conda env:holoview]
#     language: python
#     name: conda-env-holoview-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
from dask.distributed import Client
import geopandas
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
import xarray as xr

from pyPRMS.plot_helpers import set_colormap, get_projection, plot_polygon_collection

warnings.filterwarnings('ignore')

# %%
# Location of config information used for matplotlib
mpl.matplotlib_fname()

# %%
client = Client()
client

# %%
# work_dir = '/Volumes/USGS_NHM2/datasets/gridmet/gm_gf_v11_US_UNITS_filled'
work_dir = '/Volumes/USGS_NHM2/datasets/gridmet/gm_new_ncf_filled'

# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v3/tmp'
# filename = '{}/crap.nc'.format(work_dir)

geofile = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/NHM_v1.1/GIS/GFv1.1.gdb'
layer_name = 'nhru_v1_1_simp'
shape_key = 'nhru_v1_1'
# geofile = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
# layer_name = 'nhruv11_sim30'
# shape_key = 'nhru_v11'
# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v3_1980-01-01_2016-12-31/netcdf'
# filename = '{}/*.nc'.format(work_dir)

# %%

# %%
def shapefile_hrus(filename, layer_name=None):
    '''
    Read a shapefile or geodatabase that corresponds to HRUs
    '''

    hru_poly = geopandas.read_file(filename, layer=layer_name)

    if hru_poly.crs.name == 'USA_Contiguous_Albers_Equal_Area_Conic_USGS_version':
        print(f'Overriding USGS aea crs with EPSG:5070')
        hru_poly.crs = 'EPSG:5070'
#     self.__hru_shape_key = shape_key

    return hru_poly


def plot_cbh(name, hru_poly, shape_key, df, time_index=None, output_dir=None, **kwargs):
    '''
    Plot a parameter
    '''
    
    # Get extent information
    minx, miny, maxx, maxy = hru_poly.geometry.total_bounds

    time_str = ''

    if time_index is not None:
        # Get a Pandas dataframe from the netcdf dataset for a single timestep
        xdf_df = df[name][time_index].to_dataframe(name=name)
        time_str = pd.to_datetime(df['time'][time_index].values).isoformat()
    else:
        xdf_df = df[name].to_dataframe(name=name)

    # Mask out any zero values (used for checking NEGMASK)
    xdf_df[name] = xdf_df[name].mask(xdf_df[name] == 0)
    
    # cmap, norm = set_colormap(name, xdf_df[name], min_val=0.0, max_val=xdf_df.max().max(), cmap='tab20b', **kwargs)
    cmap, norm = set_colormap(name, xdf_df[name], min_val=0.0, max_val=xdf_df[name].max(), cmap='tab20b', **kwargs)
    # cmap, norm = set_colormap(name, xdf_df[name])

    # NOTE: 2022-03-09 - using set_bad will override ANY edgecolor setting.
    # cmap.set_bad('darkgrey', 1.0)
    
    crs_proj = get_projection(hru_poly)

    # This takes care of multipolygons that are in the NHM geodatabase/shapefile
    geoms_exploded = hru_poly.explode(index_parts=True).reset_index(level=1, drop=True)

    print('Writing first plot')
    df_mrg = geoms_exploded.merge(xdf_df[name], left_on=shape_key, right_index=True, how='left')

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(30, 20))

    ax = plt.axes(projection=crs_proj)
    ax.coastlines()
    ax.gridlines()
    ax.set_extent([minx, maxx, miny, maxy], crs=crs_proj)

    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapper.set_array(df_mrg[name])
    plt.colorbar(mapper, shrink=0.6)

    if time_index is not None:
        plt.title(f'Variable: {name},  Date: {time_str}')
    else:
        plt.title(f'Variable: {name}')
        
    col = plot_polygon_collection(ax, df_mrg.geometry, values=df_mrg[name],
                                  **dict(kwargs, cmap=cmap, norm=norm))

    if output_dir is not None:
        if time_index is not None:
            plt.savefig(f'{output_dir}/{name}_{time_index+1:02}.png', dpi=150, bbox_inches='tight')
        else:
            plt.savefig(f'{output_dir}/{name}.png', dpi=150, bbox_inches='tight')


# %%

# %%
# Open the netcdf file of NHM output variables
xdf = xr.open_mfdataset(f'{work_dir}/*_climate_*.nc', chunks={}, combine='by_coords', parallel=True)
xdf = xdf.assign_coords(hruid=xdf.hruid)
# xdf = xdf.set_index(hruid='hruid')

# xdf = xdf.sum(dim='time')

# %%
aa = xdf['tmax'][0].to_dataframe(name='tmax')
aa

# %%
aa['tmax'].max()

# %%
xdf

# %%
xdf['tmax'].loc[dict(time=slice('2000-01-01', '2000-01-02'), hruid=96751)].values

# %%

# %%
hru_poly = shapefile_hrus(geofile, layer_name=layer_name)

# %%
plot_cbh('tmax', hru_poly, shape_key, xdf, time_index=0, edgecolor='none', linewidth=0.0)

# %%

# %%
# %%time
xdf = xr.open_mfdataset(f'{work_dir}/*_climate_*.nc', chunks={}, combine='by_coords', parallel=True)
xdf = xdf.assign_coords(hruid=xdf.hruid)

# Create a mask; 1 = (tmax < tmin), 0 - otherwise
xdf['NEGMASK'] = xdf.tmax < xdf.tmin

xdf = xdf.drop_vars(['hru_lat', 'hru_lon', 'prcp', 'rhavg', 'rhmax', 'rhmin', 'ws', 'tmax', 'tmin'])

# Sum the mask values
xdf = xdf.sum(dim='time')
xdf

# %%

# %%
# %%time
xdf.load()

# %%

# %%

# %%
hru_poly = shapefile_hrus(geofile, layer_name=layer_name)

# %%
plot_cbh('NEGMASK', hru_poly, shape_key, xdf, edgecolor='none', linewidth=0.0)
# plot_cbh('NEGMASK', hru_poly, shape_key, xdf, time_index=360) # , output_dir='/Users/pnorton/tmp')

# %%
xdf.NEGMASK.plot.hist(bins=[1,2,5,10])

# %%

# %%

# %%

# %% [markdown]
# ## Load time-series of occurrences of tmax < tmin

# %%
# %%time
xdf_ts = xr.open_mfdataset(f'{work_dir}/gm_climate_*.nc', chunks={}, combine='by_coords', parallel=True)
xdf_ts = xdf_ts.assign_coords(hruid=xdf_ts.hruid)

# Create a mask; 1 = (tmax < tmin), 0 - otherwise
xdf_ts['NEGMASK'] = xdf_ts.tmax < xdf_ts.tmin

xdf_ts = xdf_ts.drop_vars(['hru_lat', 'hru_lon', 'prcp', 'rhavg', 'rhmax', 'rhmin', 'ws', 'tmax', 'tmin'])

# Sum the mask values
# xdf = xdf.sum(dim='time')
xdf_ts.load()

# %% [markdown]
# ## Get list of hruid where NEGMASK is greater than a threshold

# %%
# %%time
# Get the hruids where NEGMASK is greater than a threshold
thold = 10

neg_hrus = np.argwhere(xdf.NEGMASK.values > thold)
df_hrus = xdf.isel(hruid=neg_hrus[:,0])

gt_thold_hruid = df_hrus.hruid.values

print(f'Number of HRUs with at least {thold} occurrences of tmax > tmin: {gt_thold_hruid.size}')

# %% [markdown]
# ## Plot NEGMASK for a particular hruid

# %%
# %%time
chru = 102314
xdf_ts.NEGMASK.sel(hruid=chru).plot()

print(f'Number of occurrences (tmax > tmin): {xdf_ts.NEGMASK.sel(hruid=chru).sum().values}')

# %%
aa = xdf_ts.where(xdf_ts.NEGMASK > 0.0, drop=True)

# %%
aa

# %% tags=[]
# Print out the dates where tmax < tmin
for xx in aa.time.dt.strftime('%Y-%m-%d'):
    print(xx.values)

# %%



# %%

# %%

# %%

# %%
