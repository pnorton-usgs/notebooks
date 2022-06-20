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
import xarray as xr
import hvplot.xarray
import numpy as np
import warnings
import pandas as pd
from pygeohydro import NWIS, plot

warnings.filterwarnings('ignore')

# %%
work_dir = '/Volumes/USGS_NHM2/datasets/NWM'
filename = f'{work_dir}/nwmv21_nwis.nc'

filename_lookup = f'{work_dir}/nwm_linked_sites.csv'

# %%
nwis = NWIS()

# %% [markdown]
# ## Read the reach ID to streamgage lookup file

# %%
col_names = ['site_id', 'COMID', 'REACHCODE', 'REACH_meas']
col_types = [str, int, str, float]

cols = dict(zip(col_names, col_types))

df = pd.read_csv(filename_lookup, sep=',', dtype=cols)
df.head()

# %% [markdown]
# ## Read the NWM file

# %% tags=[]
df_nwm = xr.open_dataset(filename)
df_nwm

# %% [markdown]
# ## Check which reach IDs are either missing or have multiple streamgages associated with them

# %% tags=[]
# %%time
lim_distance = 1
bad_cnt = 0
multi_cnt = 0
no_gage = 0

bad_id = []
multi_id = {}
for ff in df_nwm.feature_id.values.tolist():
    try:
        cdf = df.loc[df.COMID == ff]
        
        if len(cdf) > 1:
            multi_cnt += 1
            multi_id[ff] = cdf['site_id'].tolist()
            # print(f'{ff}: {",".join(nldi_st_all.identifier.tolist())}')
            
        if len(cdf) == 0:
            no_gage += 1
    except:
        bad_cnt += 1
        bad_id.append(ff)
        
print(f'Number of missing reach IDs: {bad_cnt}')
print(f'Number of reach IDs without gages: {no_gage}')
print(f'fNumber of reach IDs with multiple streamgages: {multi_cnt}')

# %%
curr_reach = 3923

if curr_reach in multi_id:
    print(f'Reach: {curr_reach} has multiple ({len(multi_id)}) streamgage IDs. Using the first one')
    
curr_site = df['site_id'].loc[df.COMID == curr_reach].values[0]

# Plot the NWM values
df_nwm.streamflow.sel(time=slice("1980-01-01", "2020-12-31"), feature_id=curr_reach).plot()

# Plot the NWIS values
qobs = nwis.get_streamflow(curr_site, ("1980-01-01", "2020-12-31"), mmd=False)
_ = qobs.plot()

# %% [markdown]
# ## Plot the lat/lon points

# %%
# df_nwm['mask'] = df_nwm.streamflow.sel(time='1979-02-01T01:00:00')
# df_nwm['mask'][:] = 1

df_mask = df_nwm.drop_vars(['streamflow', 'crs', 'time']).to_dataframe()
df_mask.plot(x='longitude', y='latitude', kind='scatter', s=0.5)

# %% [markdown] tags=[]
# ## Ignore - experiment with using NLDI to lookup reach ID

# %%
from pygeohydro import NWIS, plot
from pynhd import NLDI
from pygeoutils import InvalidInputType

# %% tags=[]
# %%time
lim_distance = 1
bad_cnt = 0
multi_cnt = 0

bad_id = []
multi_id = {}
for ff in df.feature_id.values.tolist():
    try:
        nldi_st_all = NLDI().navigate_byid(fsource="comid", fid=str(ff), 
                                           navigation="upstreamMain", source="nwissite", distance=lim_distance)
        if len(nldi_st_all) > 1:
            multi_cnt += 1
            multi_id[ff] = nldi_st_all.identifier.tolist()
            # print(f'{ff}: {",".join(nldi_st_all.identifier.tolist())}')
    except InvalidInputType:
        bad_cnt += 1
        bad_id.append(ff)
        # print(f'ERROR: {ff}')

# %%
reach_id = '77832'
lim_distance = 1

nldi_flw_main = NLDI().navigate_byid(fsource="comid", fid=reach_id, 
                                     navigation="upstreamMain", source="flowlines", distance=lim_distance)

nldi_st_all = NLDI().navigate_byid(fsource="comid", fid=reach_id, 
                                   navigation="upstreamMain", source="nwissite", distance=lim_distance)

# %%
nldi_st_all

# %%
nldi_st_all.identifier.tolist()
