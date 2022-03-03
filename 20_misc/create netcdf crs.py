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
import netCDF4 as nc
from pyproj import CRS

# %%
workdir = '/Volumes/USGS_NHM2/datasets/USGS_monthly_water_balance_model'
out_file = f'{workdir}/crap1.nc'

# %%
crs = CRS('epsg:4326')
cf_grid_mapping = crs.to_cf()
cf_grid_mapping

# %%
# Open netcdf file in append mode
nco = nc.Dataset(out_file, 'a')

# %%
# Add the crs information
crs_var = nco.createVariable('crs', 'i4')

for kk, vv in cf_grid_mapping.items():
    crs_var.setncattr(kk, vv)
    
# Add grid_mapping attribute to the variables
var_list = ['aet', 'pet', 'prcp', 'rain', 'runoff', 'snow', 'soilstorage', 'swe', 'tmean']

for vv in var_list:
    nco[vv].grid_mapping = 'crs'
    
# Remove _FillValue attribute for lat and lon coordinate variables
nco['lat'].delncattr('_FillValue')
nco['lon'].delncattr('_FillValue')

# %%
nco.close()

# %%
