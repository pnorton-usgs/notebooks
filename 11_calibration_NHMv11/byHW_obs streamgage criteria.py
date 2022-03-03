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
hw_data_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/headwaters'
paramdb_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS'
falcone_dir = '/Users/pnorton/GIS/gagesII_additiona_data/basinchar_and_report_sept_2011'

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/poi_data'

nwis_file = f'{workdir}/*_pois.nc'

st_date = datetime.datetime(1980, 1, 1)
en_date = datetime.datetime(2019, 12, 31)

# %%
# xdf = xr.open_dataset(nwis_file, decode_cf=True, engine='netcdf4')
poi_xdf = xr.open_mfdataset(nwis_file, decode_cf=True, combine='nested', concat_dim='poi_id', engine='netcdf4')
poi_xdf

# %%
poi_xdf = xr.open_mfdataset(nwis_file, decode_cf=True, combine='nested', concat_dim='poi_id', engine='netcdf4')

poi_df = poi_xdf[['poi_name', 'latitude', 'longitude', 'drainage_area', 'drainage_area_contrib']].to_pandas()
poi_df.rename(columns={'drainage_area': 'da_obs', 'drainage_area_contrib': 'da_contrib_obs'}, inplace=True)

# Compute the correction factor for obs values
poi_df['da_ratio_obs'] = poi_df['da_contrib_obs'] / poi_df['da_obs']

# Any NaN values should default to a correction factor of 1.0
poi_df['da_ratio_obs'].fillna(value=1.0, inplace=True)

# Sometimes the full da and contributing da are swapped
poi_df['da_ratio_obs'] = np.where(poi_df['da_ratio_obs'] > 1.0, 1.0 / poi_df['da_ratio_obs'], poi_df['da_ratio_obs'])

poi_df['da_actual_obs'] = poi_df[['da_obs', 'da_contrib_obs']].min(axis=1)
poi_df.head()

# %%
df1 = poi_xdf['discharge'].loc[:, st_date:en_date].to_pandas()

# Get POIs which have no missing values in POR
missing_obs_cnt_df = df1.isnull().sum(axis=1)
missing_obs_cnt_df.name = 'missing_obs'

poi_contiguous = missing_obs_cnt_df[missing_obs_cnt_df == 0].index.tolist()
print(f'POIs with contiguous observations: {len(poi_contiguous)}')

# %%

missing_obs_cnt_df.head()

# %%
# da_nwis = xdf['drainage_area'].to_pandas()
# print(da_nwis.head())

# Convert area to acres
# da_nwis *= 640.
# print(da_nwis.head())

# %%
# da_stuff = xdf[['drainage_area', 'drainage_area_contrib']].to_pandas()

# %%
# da_stuff.head(50)

# %%
missing_obs_cnt_df['02157490']

# %% [markdown]
# ### Load headwater segments

# %%
segs_by_hw = read_file(f'{hw_data_dir}/hw_segs.csv')

segs_list = []
seg_to_hw = {}

for kk, vv in segs_by_hw.items():
    # kk is the headwater number
    # segs_list.append(kk)
    
    for xx in vv:
        seg_to_hw[xx] = kk
        segs_list.append(xx)
        
# Get set of unique segment ID
segs_set = set(segs_list)

print(f'Total number of segments: {len(segs_list)}')
print(f'Unique number of segments: {len(segs_set)}')

# %%

# %% [markdown]
# ### Load parameter database

# %%
pdb = ParamDb(paramdb_dir, verbose=True, verify=True)

poi_to_seg = pdb.parameters.poi_to_seg
seg_to_poi = {vv: kk for kk, vv in poi_to_seg.items()}

poi_gage_id = pdb.parameters.get('poi_gage_id').data.tolist()
seg_cum_area = pdb.parameters.get_dataframe('seg_cum_area')

# Create list of POIs in headwater domains
hw_poi_list = []

for xx in segs_set:
    try:
        hw_poi_list.append(seg_to_poi[xx])
    except KeyError:
        pass
        # print(f'{xx} has no POI')
        
print(f'Number of POIs in headwater areas: {len(hw_poi_list)}')



# Create list of POIs with the NHM drainage area
poi_areas = {'poi_id': [],
             'hw_id': [],
             'da_seg_cum': []}

for xx in poi_gage_id:
    try:
        if xx in hw_poi_list:
            # Convert NHM acres to sq. mi.
            poi_areas['poi_id'].append(xx)
            poi_areas['hw_id'].append(seg_to_hw[poi_to_seg[xx]])
            
            if poi_to_seg[xx] > 0:
                poi_areas['da_seg_cum'].append(seg_cum_area.loc[poi_to_seg[xx]].values[0] * 0.0015625)
            else:
                # I think this only happened in NHM v1.0
                poi_areas['da_seg_cum'].append(0)
    except KeyError:
        print(f'{xx} has no POI')
        
        
# Create a dataframe of NHM drainage by POI
poi_areas_df = pd.DataFrame.from_dict(poi_areas)
poi_areas_df.set_index('poi_id', inplace=True)
poi_areas_df.head() 

# %%
# Merge POI information with seg_cum_area
poi_info_df = poi_df
poi_info_df = pd.merge(poi_info_df, poi_areas_df, how='inner', left_index=True, right_index=True)
poi_info_df = pd.merge(poi_info_df, missing_obs_cnt_df, how='inner', left_index=True, right_index=True)
poi_info_df.head()

# Add drainage area ratio between nwis/hydat and NHM
poi_info_df['da_ratio'] = poi_info_df['da_actual_obs'] / poi_info_df['da_seg_cum']
poi_info_df['da_ratio'] = np.where(poi_info_df['da_ratio'] > 1.0, 1.0 / poi_info_df['da_ratio'], poi_info_df['da_ratio'])

# %%
poi_info_df.info()

# %%

# %%

# %% [markdown]
# ### Get POIs in headwaters that have contiguous records

# %%
hw_contig = set(hw_poi_list) & set(poi_contiguous)
print(f'Number of contiguous POIs in headwaters: {len(hw_contig)}')

# %%

# %%

# %% [markdown]
# ### Load Falcone information

# %%
col_names = ['STAID', 'CLASS', 'HYDRO_DISTURB_INDX']
col_types = [str, str, int]
cols = dict(zip(col_names, col_types))

falcone_df = pd.read_excel(open(f'{falcone_dir}/gagesII_sept30_2011_conterm.xlsx', 'rb'), sheet_name='Bas_Classif', 
                           usecols=[0, 1], dtype=cols)

falcone_df.rename(columns={'STAID': 'poi_id', 'CLASS': 'falcone_class'}, inplace=True)
falcone_df.set_index('poi_id', inplace=True)
falcone_df.info()

falcone_ids = falcone_df.index.tolist()
falcone_ref = falcone_df[falcone_df['falcone_class'] == 'Ref']
falcone_ref_ids = falcone_ref.index.tolist()

# %%

# %%
# Headwater POIs that are in Falcone dataset
hw_contig_falcone = set(falcone_ids) & hw_contig
print(f'Number of headwater POIs also in Falcone: {len(hw_contig_falcone)}')

hw_contig_ref = set(falcone_ref_ids) & hw_contig
print(f'Number of headwater POIs that are Falcone reference gages: {len(hw_contig_ref)}')

# %%
len(falcone_ref_ids)

# %%
poi_info_df = pd.merge(poi_info_df, falcone_df, how='inner', left_index=True, right_index=True)
poi_info_df.head()

# %%
poi_info_df.info()

# %%

# %%
# Remove POIs which are missing more than a given number of obs
df_reduce_1 = poi_info_df[poi_info_df['missing_obs'] < 3650]

# %%
df_reduce_1.info()

# %%
# Remove POIs that lack a DA 
df_reduce_2 = df_reduce_1[df_reduce_1['da_actual_obs'].notna()]

# %%
df_reduce_2.info()

# %%
df_reduce_3 = df_reduce_2[df_reduce_2['da_ratio'] >= 0.9]

# %%
df_reduce_3.info()

# %%
df_reduce_3.head()

# %%
df_reduce_3.to_csv('/Users/pnorton/tmp/nhm_v11_hwobs_1.csv', sep='\t', index=True)

# %%
df_reduce_1.loc['02160105']

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
hw_id_set = set(df_reduce_3['hw_id'].tolist())

# %%
len(hw_id_set)

# %%
# poi_areas_df = pd.DataFrame.from_dict(poi_areas)
# poi_areas_df.set_index('poi_id', inplace=True)

hwid = {'hw_id': list(hw_id_set)}

hwid_df = pd.DataFrame.from_dict(hwid)
hwid_df['use_obs'] = 1
hwid_df.set_index('hw_id', inplace=True)
hwid_df.head()

# %%
hwid_df.to_csv('/Users/pnorton/tmp/nhm_v11_hwobsid_1.csv', sep='\t', index=True)

# %%
hw2_df = pd.read_csv(f'{hw_data_dir}/hw_hrus.csv')
hw2_df.set_index('nhm_id', inplace=True)
hw2_df.info()

# %%
hw2_df = pd.merge(hw2_df, hwid_df, how='inner', left_index=False, left_on='hw_id', right_index=True)
hw2_df.head()


# %%
hw2_df.info()

# %%
hw2_df.to_csv('/Users/pnorton/tmp/nhm_v11_byHWobs_hrus.csv', sep='\t', index=True)

# %% [markdown]
# ## Create a hw_obs_segs file

# %%
segs_by_hw = read_file(f'{hw_data_dir}/hw_segs.csv')
hw_obs_fhdl = open(f'{hw_data_dir}/hw_obs_segs.csv', 'w')

hw_obs_fhdl.write('hw_id,start_seg,child_segs\n')

for kk, vv in segs_by_hw.items():
    # kk is the headwater number
    # segs_list.append(kk)
    if kk in hw_id_set:
        segs_str = [str(x) for x in vv]
        hw_obs_fhdl.write(f'{kk},{",".join(segs_str)}\n')

hw_obs_fhdl.close()

# %%

# %%

# %%

# %%
