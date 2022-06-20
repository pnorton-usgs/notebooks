# ---
# jupyter:
#   jupytext:
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

# %%
import pandas as pd
import numpy as np
import datetime
import xarray as xr

# %%
base_dir = '/Volumes/USGS_NHM2/NHM/NHM_v11/releases/gm_byHWobs_5.2.1'

varname = 'seginc_gwflow'
filename = f'{base_dir}/{varname}.nc'

output_file = '/Users/pnorton/tmp/test.csv'

# %% [markdown]
# ## Extract subset of segments from PRMS model output netcdf file

# %%
df = xr.open_dataset(filename, chunks={})

# This assigns the national nhm_id as the coordinate variable
df = df.assign_coords(nsegment=df.nhm_seg)
df

# %% [markdown]
# ### Define list of segments to extract

# %%
seg_list = [57478, 57467, 51960, 51957, 53016, 53025, 53563, 51140, 54561, 53630, 53798, 53792, 54970, 
            55007, 53558, 30916, 47790, 47932, 31181, 48615, 52940, 51366, 30626, 30619, 31588, 52201, 
            28798, 28631, 52210, 28801, 52181, 31559, 53114, 51988, 48429, 31718, 53434, 48500, 49106, 
            31774, 31883, 31887, 48528, 28379, 26709, 44426, 26715, 44577, 46397, 46976, 46930]

# %% [markdown]
# ### Extract selected segments for either the full available time or a date range

# %%
# %%time
# Pull subset of the segments for all timesteps
# df_ss = df.sel(dict(nsegment=seg_list))
# df_ss

# Pull subset of segment for a date range
st = datetime.datetime(2000, 1, 1)
en = datetime.datetime(2001, 12, 31)

df_ss = df.sel(dict(time=slice(st, en), nsegment=seg_list))
df_ss

# %% [markdown]
# ### Convert xarray to a pandas dateframe

# %%
# %%time
# Create a pandas dataframe of the variable you are interested in
df_pd = df_ss[varname].to_pandas()
df_pd

# %% [markdown]
# ### Write the results to a CSV file

# %%
# Output to a CSV file, use full floating point values
df_pd.to_csv(output_file, columns=seg_list, sep=',', index=True, header=True, chunksize=50)

# Output to a CSV file, reduce the precision of the floating point values
# df_pd.to_csv(output_file, columns=seg_list, sep=',', index=True, header=True, chunksize=50, float_format='%0.4f')

# %% [markdown]
# ### Plotting examples

# %%
# Plot all time steps for a single segment
df_pd.loc[:, 30916].plot()

# %%
# Plot a date range for a single segment
st = datetime.datetime(2000, 1, 1)
en = datetime.datetime(2000, 12, 31)

df_pd.loc[st:en, 30916].plot()

# %%

# %%

# %%
