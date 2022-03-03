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
#     display_name: Python [conda env:bandit_38]
#     language: python
#     name: conda-env-bandit_38-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import geopandas
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import xarray as xr

from pyPRMS.plot_helpers import set_colormap, get_projection, plot_polygon_collection

# %%
# workdir = '/Volumes/USGS_NHM1/gm_gf_v11_US_with_derived'
workdir = '/Volumes/USGS_NHM1/gm_gf_v11_US_filled'

filename = f'{workdir}/gm_climate_1979.nc'

src_fill_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/cbh_work/onhm-runners/cbh_filler'

map_filename = f'{src_fill_dir}/miss_to_pres_mapping.csv'


hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
hru_layer_name = 'nhruv11_sim30'
hru_shape_key='nhru_v11'

# %%
xdf = xr.open_mfdataset(filename, combine='by_coords', decode_cf=False, engine='netcdf4')

# %%
xdf

# %%

# %% [markdown]
# ### Read in the mapping file

# %%
map_df = pd.read_csv(map_filename, usecols=[1, 3], index_col=0)

# %%
map_df.info()

# %%
map_df.head()

# %%
# Get dict mapping missing IDs to fill IDs
miss_to_fill = map_df.to_dict()['nhru_v11_pres']

# %%
len(miss_to_fill)

# %%
fill_vars = ['tmin', 'tmax', 'prcp', 'rhavg']

for vv in fill_vars:
    print(vv)
    tmp_a = xdf[vv].values
    
    for mm, ff in miss_to_fill.items():
        tmp_a[:, mm-1] = tmp_a[:, ff-1]
        
    xdf[vv].values = tmp_a        

# %%
encoding = {'time': {'_FillValue': None},
            'hru_lon': {'_FillValue': None},
            'hru_lat': {'_FillValue': None}}

xdf.to_netcdf('/Users/pnorton/tmp/crap.nc', encoding=encoding)

# %%

# %% [markdown]
# ### Read in shapefile

# %%
hru_poly = geopandas.read_file(hru_geodatabase, layer=hru_layer_name)

if hru_poly.crs.name == 'USA_Contiguous_Albers_Equal_Area_Conic_USGS_version':
    print('Overriding USGS aea crs with EPSG:5070')
    hru_poly.crs = 'EPSG:5070'


# %%
def plot_cbh(name, hru_poly, shape_key, df, time_index=None, output_dir=None, **kwargs):
    '''
    Plot a parameter
    '''
    
    # Get extent information
    minx, miny, maxx, maxy = hru_poly.geometry.total_bounds

    # Get a Pandas dataframe from the netcdf dataset for a single timestep
    xdf_df = df[name][time_index].to_dataframe(name=name)

    # cmap, norm = set_colormap(name, param_var, min_val=cparam.minimum, max_val=cparam.maximum)
    cmap, norm = set_colormap(name, xdf_df[name])

    # Fill value = 9.96920996838687e+36
    xdf_df.mask(xdf_df == 9.96920996838687e+36)
    cmap.set_bad('darkgrey', 0.7)
    
    crs_proj = get_projection(hru_poly)

    # This takes care of multipolygons that are in the NHM geodatabase/shapefile
    geoms_exploded = hru_poly.explode().reset_index(level=1, drop=True)

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

    time_str = pd.to_datetime(df['time'][time_index].values).isoformat()
#     time_str = df['time'].iloc[time_index].isoformat()
    plt.title(f'Variable: {name},  Date: {time_str}')    

    col = plot_polygon_collection(ax, df_mrg.geometry, values=df_mrg[name], **dict(kwargs, cmap=cmap, norm=norm))    


# %%
plot_cbh('tmin', hru_poly, hru_shape_key, xdf, time_index=360)

# %%

# %%

# %%
