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

# %% [markdown]
# # NHGF Tools and applications demo: Creating a baseline dataset

# %%
import geopandas as gpd
import pickle
import time
import xagg as xa
import xarray as xr

# %% [markdown]
# Using soil moisture from the NLDAS Mosaic Land Surface Model L4 Monthly 0.125 x 0.125 degree V002 (NLDAS_MOS0125_M) product.
#
# ![title](NLDAS_MOS0125_M_002.png)

# %% [markdown]
# website link: https://disc.gsfc.nasa.gov/datasets/NLDAS_MOS0125_M_002/summary?keywords=NLDAS

# %%

# %%

# %%

# %%

# %%
nldas_dir = '/Volumes/USGS_NHM2/datasets/NLDAS/NLDAS_MOSAIC_L4_monthly'
# filename example: NLDAS_MOS0125_M.A200007.002.grb.SUB.nc4

geofile = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
layer_name = 'nhruv11_sim30'
shape_key = 'nhru_v11'

# %%
# Read the HRU geodatabase
gdf = gpd.read_file(geofile, layer=layer_name)


# NOTE: What's the projection of this product?
# print(f'Re-projecting {layer_name} to WGS84 (lat/lon)')
# gdf = gdf.to_crs(epsg=4326)

# %%
# Open a single month for generating the weights
ds_single = xr.open_dataset(f'{nldas_dir}/NLDAS_MOS0125_M.A200007.002.grb.SUB.nc4')

# %%
# uncomment below if first time through notebook to generate weights

start = time.perf_counter()
weightmap = xa.pixel_overlaps(ds_single,gdf)
end = time.perf_counter()
print(f'finished agg in {round(end-start, 2)} second(s)')

# save weights for future use
with open('nldas_weights.pickle', 'wb') as file:
    pickle.dump(weightmap, file)

# %%

# %%

# %%
ds_nldas = xr.open_mfdataset(f'{nldas_dir}/NLDAS_MOS0125_M.*', decode_cf=True, engine='netcdf4')

# xr.open_mfdataset(self.__src_path,
#                                                        chunks={'hruid': 1040}, combine='by_coords',
#                                                        decode_cf=True, engine='netcdf4')

# %%
ds_nldas

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
