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
import datetime
import xarray as xr
import pandas as pd
import numpy as np

from Bandit.bandit_multi_locations import read_file
from pyPRMS.ParamDb import ParamDb

# %%
base_dir = '/Users/pnorton/Projects/National_Hydrology_Model'
paramdb_dir = f'{base_dir}/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS'
falcone_dir = f'/Users/pnorton/GIS/gagesII_additiona_data/basinchar_and_report_sept_2011'
nwis_dir = f'{base_dir}/datasets/bandit/poi_data'

output_filename = f'{base_dir}/calibrations/NHMv11/byHW_obs_sample/v11_calib_info/poi_info.csv'

nwis_file = f'{nwis_dir}/*_pois.nc'

st_date = datetime.datetime(1983, 1, 1)
en_date = datetime.datetime(2010, 12, 31)

# %% [markdown]
# ### Load POI station information

# %%
nwis_xdf = xr.open_mfdataset(nwis_file, decode_cf=True, combine='nested', concat_dim='poi_id', engine='netcdf4')
nwis_xdf

# %%
nwis_use_df = nwis_xdf[['poi_name', 'latitude', 'longitude', 'drainage_area', 'drainage_area_contrib']].to_pandas()
nwis_use_df.rename(columns={'drainage_area': 'da_obs', 'drainage_area_contrib': 'da_contrib_obs'}, inplace=True)
nwis_use_df.head()

# %%

# %% [markdown]
# ### Load parameter database

# %%
pdb = ParamDb(paramdb_dir, verbose=True, verify=True)

poi_to_seg = pdb.parameters.poi_to_seg
seg_to_poi = {vv: kk for kk, vv in poi_to_seg.items()}

# %%
poi_gage_id = pdb.parameters.get('poi_gage_id').data.tolist()
seg_cum_area = pdb.parameters.get_dataframe('seg_cum_area')

# %%
# Create list of POIs with NHM drainage area
poi_areas = {'poi_id': [],
             'da_nhm': []}

for xx in poi_gage_id:
    try:
        # Convert NHM acres to sq. mi.
        poi_areas['poi_id'].append(xx)
        poi_areas['da_nhm'].append(seg_cum_area.loc[poi_to_seg[xx]].values[0] * 0.0015625)
    except KeyError:
        print(f'{xx} has no POI')

# %%
# Create a dataframe of NHM drainage by POI
poi_areas_df = pd.DataFrame.from_dict(poi_areas)
poi_areas_df.set_index('poi_id', inplace=True)
poi_areas_df.head()

# %% [markdown]
# ### Load Falcone information

# %%
# col_names = ['STAID', 'CLASS', 'HYDRO_DISTURB_INDX']
# col_types = [str, str, int]
col_names = ['STAID', 'CLASS']
col_types = [str, str]
cols = dict(zip(col_names, col_types))

falcone_df = pd.read_excel(open(f'{falcone_dir}/gagesII_sept30_2011_conterm.xlsx', 'rb'), sheet_name='Bas_Classif', 
                           usecols=[0, 1], dtype=cols)
falcone_df.rename(columns={'STAID': 'poi_id', 'CLASS': 'falcone_class'}, inplace=True)
falcone_df.set_index('poi_id', inplace=True)

falcone_df.head()

# %%

# %%
# Create a dataframe which includes obs, falcone, and sim information
poi_info_df = nwis_use_df
poi_info_df = pd.merge(poi_info_df, falcone_df, how='left', left_index=True, right_index=True)
poi_info_df = pd.merge(poi_info_df, poi_areas_df, how='left', left_index=True, right_index=True)

# %%

# %%

# %% [markdown]
# ### Add computed fields

# %%
# Compute the correction factor for obs values
poi_info_df['da_corrfact_obs'] = poi_info_df['da_contrib_obs'] / poi_info_df['da_obs']

# Any NaN values should default to a correction factor of 1.0
poi_info_df['da_corrfact_obs'].fillna(value=1.0, inplace=True)

# %%
# Compute the actual drainage area (smallest value of da_obs and da_contrib_obs)
poi_info_df['da_actual_obs'] = poi_info_df[['da_obs', 'da_contrib_obs']].min(axis=1)

# %%
# Compute the model correction factor
poi_info_df['da_corrfact_sim'] = poi_info_df['da_actual_obs'] / poi_info_df['da_nhm']

# %%

# %%
# Write the file
poi_info_df.to_csv(output_filename, sep='\t', index=True)

# %%
poi_info_df[poi_info_df[da_corrfact_sim]]
