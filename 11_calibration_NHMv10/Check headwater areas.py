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

# %%
base_dir = '/Users/pnorton/Projects/National_Hydrology_Model'
paramdb_dir = f'{base_dir}/datasets/paramdb_v10/paramdb_v10_daymet_CONUS'
falcone_dir = f'/Users/pnorton/GIS/gagesII_additiona_data/basinchar_and_report_sept_2011'
nwis_dir = f'{base_dir}/datasets/bandit/poi_data'

hw_data_dir = f'{base_dir}/datasets/headwaters_conus'

gagearea_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHW_obs_sample/helpers'
gageuse_file = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHW_obs_sample/byHW_obs_GAGEuse.txt'

output_filename = f'{base_dir}/calibrations/NHMv10/crap_poi_info.csv'

hw_file = f'{hw_data_dir}/2018-12-17_hwSegsVll.csv'

nwis_file = f'{nwis_dir}/*_pois.nc'

st_date = datetime.datetime(1981, 1, 1)
en_date = datetime.datetime(2010, 12, 31)

# %% [markdown]
# ### Load headwater file used for extractions

# %%
segs_by_hw = read_file(hw_file)

segs_list = []
seg_to_hw = {}

for kk, vv in segs_by_hw.items():
    # kk is the headwater ID
    # segs_list.append(kk)
    
    for xx in vv:
        seg_to_hw[xx] = kk
        segs_list.append(xx)
        
segs_set = set(segs_list)

print(len(segs_list))
print(len(segs_set))

# %%

# %%
segs_to_hw_dict = {}

for kk, vv in segs_by_hw.items():
    for xx in vv:
        segs_to_hw_dict[xx] = kk

print(len(segs_to_hw_dict))

# %%
segs_to_hw_dict

# %%

# %%

# %%

# %% [markdown]
# ### Load POI station information

# %%
poi_xdf = xr.open_mfdataset(nwis_file, decode_cf=True, combine='nested', concat_dim='poi_id', engine='netcdf4')

poi_df = poi_xdf[['poi_name', 'latitude', 'longitude', 'drainage_area', 'drainage_area_contrib']].to_pandas()
poi_df.rename(columns={'drainage_area': 'da_obs', 'drainage_area_contrib': 'da_contrib_obs'}, inplace=True)

poi_df['da_actual_obs'] = poi_df[['da_obs', 'da_contrib_obs']].min(axis=1)
poi_df.head()

# %%
poi_df.info()

# %%
df1 = poi_xdf['discharge'].loc[:, st_date:en_date].to_pandas()
df1.head()

# %%
# Get POIs which have no missing values in POR
miss_obs_count_df = df1.isnull().sum(axis=1)
poi_contiguous = miss_obs_count_df[miss_obs_count_df < 5800].index.tolist()
print(f'POIs with contiguous observations: {len(poi_contiguous)}')

# %%
miss_obs_count_df.head()

# %%
miss_obs_count_df['01031450']

# %%
poi_contiguous

# %%

# %%

# %%

# %% [markdown]
# ### Load parameter database

# %%
pdb = ParamDb(paramdb_dir, verbose=True, verify=True)

poi_to_seg = pdb.parameters.poi_to_seg
seg_to_poi = {vv: kk for kk, vv in poi_to_seg.items()}

poi_gage_id = pdb.parameters.get('poi_gage_id').data.tolist()
seg_cum_area = pdb.parameters.get_dataframe('seg_cum_area')

# Create list of POIs with NHM drainage area
poi_areas = {'poi_id': [],
             'da_seg_cum': []}

for xx in poi_gage_id:
    try:
        # Convert NHM acres to sq. mi.
        poi_areas['poi_id'].append(xx)
        
        if poi_to_seg[xx] > 0:
            poi_areas['da_seg_cum'].append(seg_cum_area.loc[poi_to_seg[xx]].values[0] * 0.0015625)
        else:
            poi_areas['da_seg_cum'].append(0)
    except KeyError:
        print(f'{xx} has no POI')
        
        
# Create a dataframe of NHM drainage by POI
poi_areas_df = pd.DataFrame.from_dict(poi_areas)
poi_areas_df.set_index('poi_id', inplace=True)
poi_areas_df.head()        

# %%
poi_areas_df.head(50)

# %%
# Build list of POIs for headwater areas
hw_poi_list = []

for xx in segs_set:
    try:
        hw_poi_list.append(seg_to_poi[xx])
    except KeyError:
        pass
        # print(f'{xx} has no POI')

# %%
len(hw_poi_list)

# %%
poi_areas_df.loc['01031300']
# '01031300' in poi_gage_id

# %%

# %%

# %% [markdown]
# ### Load GAGEareas

# %%
col_names = ['headwater', 'poi_id', 'falcone_class', 'disturb_idx', 'da_nwis', 'da_nhm', 'da_gf']
col_types = [str, str, str, "Int64", float, float, float]
cols = dict(zip(col_names, col_types))

gage_area_df = pd.read_csv(f'{gagearea_dir}/GAGEareas', sep=' ', names=col_names, dtype=cols, usecols=[0, 1, 2, 4, 5])
gage_area_df.set_index('poi_id', inplace=True)

# %%
gage_area_df

# %%
# gage_area_df[gage_area_df['da_nhm'] != gage_area_df['da_gf']]
# Only 11 occurences; tiny differences

# %%
# Combine POI area information and GAGEareas
# poi_info_df = poi_areas_df
# poi_info_df = pd.merge(poi_info_df, gage_area_df, how='left', left_index=True, right_index=True)
# poi_info_df.head()


# OR - combine POI area information and poi_df
poi_info_df = poi_areas_df
poi_info_df = pd.merge(poi_info_df, poi_df, how='left', left_index=True, right_index=True)
poi_info_df.head()

# %%
poi_info_df.info()

# %%
# Select rows where computed cumulative drainage area matches GAGEareas da_nhm
# This removes POIs where da_nhm is NaN
# poi_da_ok_df = poi_info_df[poi_info_df['da_seg_cum'] - poi_info_df['da_nhm'] < 0.001]
# poi_da_ok_df['da_ratio'] = poi_da_ok_df['da_nwis'] / poi_da_ok_df['da_seg_cum']

# poi_da_ok_df['da_ratio'] = np.where(poi_da_ok_df['da_ratio'] > 1.0, 1.0 / poi_da_ok_df['da_ratio'], poi_da_ok_df['da_ratio'])

# poi_da_ok_df

# %%
poi_info_df['da_ratio'] = poi_info_df['da_actual_obs'] / poi_info_df['da_seg_cum']
poi_info_df['da_ratio'] = np.where(poi_info_df['da_ratio'] > 1.0, 1.0 / poi_info_df['da_ratio'], poi_info_df['da_ratio'])

# %%
poi_info_df.head()

# %%

# %%
# POIs where gage drainage area and nhm area are within a percentage
poi_da_ok_df[poi_da_ok_df['da_ratio'] >= 0.50]

# %%
poi_da_ok_df[poi_da_ok_df['falcone_class'] == 'Ref'].info()

# %%
poi_info_df.loc['01031300']

# %% [markdown]
# ### Load GAGEuse file

# %%
col_names = ['headwater', 'hw_gage_num', 'poi_id', 'val_1', 'val_2', 'da_ratio', 'class', 'use_flag']
col_types = [str, "Int64", str, float, float, float, str, "Int64"]
cols = dict(zip(col_names, col_types))

gage_use_df = pd.read_csv(f'{gageuse_file}', sep='\s+', names=col_names, dtype=cols, usecols=[1, 2, 3, 4, 7])
gage_use_df.set_index('poi_id', inplace=True)

# %%
gage_use_df

# %%
# gage_use_df[gage_use_df['use_flag'] != 0].index.tolist()

# %%
poi_da_ok_df = pd.merge(poi_da_ok_df, gage_use_df, how='left', left_index=True, right_index=True)

# %%
poi_da_ok_df.loc['01031300']

# %%
pd.set_option('display.max_rows', None)
poi_da_ok_df[poi_da_ok_df['use_flag'].isnull()].info()

# %%
poi_da_ok_df[poi_da_ok_df['falcone_class'] == 'Ref'].count()

# %%
aa = poi_da_ok_df[poi_da_ok_df['use_flag'] == 0]
aa[aa['da_ratio'] >= 0.90].info()

# %%

# %% [markdown]
# ### Load Falcone gages

# %%
col_names = ['STAID', 'CLASS', 'HYDRO_DISTURB_INDX']
col_types = [str, str, int]
cols = dict(zip(col_names, col_types))

falcone_df = pd.read_excel(open(f'{falcone_dir}/gagesII_sept30_2011_conterm.xlsx', 'rb'), sheet_name='Bas_Classif', 
                           usecols=[0, 1], dtype=cols)

falcone_df.rename(columns={'STAID': 'poi_id', 'CLASS': 'falcone_class'}, inplace=True)
falcone_df.set_index('poi_id', inplace=True)
falcone_df.info()

# %%
falcone_df[falcone_df['falcone_class'] == 'Ref'].info()

# %%
falcone_df.head()

# %%
poi_info_df = pd.merge(poi_info_df, falcone_df, how='left', left_index=True, right_index=True)
poi_info_df.info()

# %%
poi_info_df[poi_info_df['poi_name'].isnull()]

# %%

# %%
falcone_all_list = falcone_df.index.tolist()

hw_in_falcone_all = []

for pp in hw_poi_list:
    if pp in falcone_all_list:
        hw_in_falcone_all.append(pp)
        
print(len(hw_in_falcone_all))

# %%
falcone_ref_list = falcone_df[falcone_df['falcone_class'] == 'Ref'].index.tolist()
# len(falcone_ref_list)

hw_in_falcone_ref = []

for pp in hw_poi_list:
    if pp in falcone_ref_list:
        hw_in_falcone_ref.append(pp)
        
print(len(hw_in_falcone_ref))

# %%
falcone_nonref_list = falcone_df[falcone_df['falcone_class'] == 'Non-ref'].index.tolist()
# len(falcone_ref_list)

hw_in_falcone_nonref = []

for pp in hw_poi_list:
    if pp in falcone_nonref_list:
        hw_in_falcone_nonref.append(pp)
        
print(len(hw_in_falcone_nonref))

# %%
# pdo_list = poi_da_ok_df[poi_da_ok_df['da_ratio'] >= 0.90].index.tolist()
pdo_list = poi_info_df[poi_info_df['da_ratio'] >= 0.90].index.tolist()

# %%
# Headwaters that have a_ratio >= 0.9
hw_da_ninety_list = []

for pp in hw_poi_list:
    if pp in pdo_list:
        hw_da_ninety_list.append(pp)
        
print(len(hw_da_ninety_list))

# %%

# %%
# Size of set that includes both da_ratio >= 0.9 and contiguous observations
poi_for_calib1 = set(pdo_list) & set(poi_contiguous)
print(len(poi_for_calib1))


# Remove pois that have gage_use > 0
gu_remove = gage_use_df[gage_use_df['use_flag'] != 0].index.tolist()

poi_for_calib2 = poi_for_calib1 - set(gu_remove)
print(len(poi_for_calib2))

# %%
aa[aa['falcone_class'] == 'Ref'].info()

# %%
aa[aa['falcone_class'] == 'Non-ref'].info()

# %%

# %%
# Map POIs back to headwater IDs
final_hw = {}
missing_from_hw = []

# for kk, vv in segs_by_hw.items():
for pp in poi_for_calib2:
    # Lookup the segment for the POI
    pseg = poi_to_seg[pp]
    
    try:
        if segs_to_hw_dict[pseg] not in final_hw:
            final_hw[segs_to_hw_dict[pseg]] = []

        final_hw[segs_to_hw_dict[pseg]].append(pp)
    except KeyError:
        print(f'{pp}, seg: {pseg} not in headwaters')
        missing_from_hw.append(pp)

# %%
len(missing_from_hw)

# %%
len(final_hw)

# %%
final_hw

# %%
hw_a = list(final_hw.keys())
hw_a.sort()

# %%
len(hw_a)

# %%
chk_file = open(f'{base_dir}/calibrations/NHMv10/byHW_obs_list.txt', 'r')

orig_hw = []

for xx in chk_file:
    orig_hw.append(int(xx.strip()))
    
len(orig_hw)

# %%
# Headwaters in the original headwater list but not in computed list
len(set(orig_hw) - set(hw_a))

# %%
# Headwaters in the computed list but not the original list
len(set(hw_a) - set(orig_hw))

# %%
# Headwaters common to both computed and original headwater lists
len(set(hw_a) & set(orig_hw))

# %%
set(hw_a) & set(orig_hw)

# %%
final_hw[27]

# %%

# %%
set(orig_hw) - set(hw_a)

# %%

# %%
# Read the byHW_obs 
hwobs_df = pd.read_csv(f'{base_dir}/calibrations/NHMv10/byHW_obs_list.txt')

hwobs_df['use_obs'] = 1
hwobs_df.set_index('hw_id', inplace=True)
hwobs_df.head()

# %%

# %%
new_segs_by_hw = read_file(f'{base_dir}/datasets/paramdb_v10/headwaters/hw_segs.csv')

new_segs_list = []

for kk, vv in new_segs_by_hw.items():
    # kk is the headwater ID
    # segs_list.append(kk)
    
    for xx in vv:
        new_segs_list.append(xx)
        
new_segs_set = set(new_segs_list)

print(len(new_segs_list))
print(len(new_segs_set))

# %%
print(len(segs_list))
print(len(segs_set))

# %%
segs_set - new_segs_set

# %%
new_segs_set - segs_set

# %%

# %%
seg_to_hru = pdb.parameters.seg_to_hru

if 0 in seg_to_hru:
    del seg_to_hru[0]

# %%
len(seg_to_hru)

# %%
seg_to_hru[8]

# %%

# %%

# %%

# %%
# Create list of POIs with the NHM drainage area
hw_hru = {'nhm_id': [],
          'hw_id': []}

for xx in segs_list:
    try:
        # hw_hru['nhm_id'].append(seg_to_hru[xx])
        
        for yy in seg_to_hru[xx]:
            hw_hru['nhm_id'].append(yy)    
            hw_hru['hw_id'].append(seg_to_hw[xx])
    except KeyError:
        print(f'seg: {xx} has no HRUs')
#     try:
#         if xx in hw_poi_list:
#             # Convert NHM acres to sq. mi.
#             poi_areas['poi_id'].append(xx)
#             poi_areas['hw_id'].append(seg_to_hw[poi_to_seg[xx]])
            
#             if poi_to_seg[xx] > 0:
#                 poi_areas['da_seg_cum'].append(seg_cum_area.loc[poi_to_seg[xx]].values[0] * 0.0015625)
#             else:
#                 # I think this only happened in NHM v1.0
#                 poi_areas['da_seg_cum'].append(0)
#     except KeyError:
#         print(f'{xx} has no POI')
        
        
# Create a dataframe of NHM drainage by POI
hw_hru_df = pd.DataFrame.from_dict(hw_hru)
hw_hru_df.set_index('nhm_id', inplace=True)
hw_hru_df.head() 

# %%
hw_hru_df.info()

# %%
hw_hru_df.to_csv('/Users/pnorton/tmp/nhm_v10_byHWobs_hrus.csv', sep='\t', index=True)

# %%
hw_hru_df = pd.merge(hw_hru_df, hwobs_df, how='inner', left_index=False, left_on='hw_id', right_index=True)

# %%
hw_hru_df.info()

# %%
hw_hru_df.to_csv('/Users/pnorton/tmp/nhm_v10_hw_to_nhmid.csv', sep='\t', index=True)

# %%
