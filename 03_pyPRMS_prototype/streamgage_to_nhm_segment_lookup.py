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
#     display_name: Python [conda env:bandit]
#     language: python
#     name: conda-env-bandit-py
# ---

# %%
import datetime
import msgpack
import numpy as np
import pandas as pd
import bandit_cfg as bc


# %%
# From bandit.py
def get_parameter(filename):
    with open(filename, 'rb') as ff:
        return msgpack.load(ff, use_list=False)


# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20170914_GreatPlains'

config = bc.Cfg('{}/bandit.cfg'.format(workdir))

# Location of the merged parameter database
merged_paramdb_dir = config.merged_paramdb_dir

poi_gage_segment = get_parameter('{}/poi_gage_segment.msgpack'.format(merged_paramdb_dir))['data']
poi_gage_id = get_parameter('{}/poi_gage_id.msgpack'.format(merged_paramdb_dir))['data']

# Create dictionary to lookup nhm_segment for a given poi_gage_id
poi_id_to_seg = dict(zip(poi_gage_id, poi_gage_segment))

# %%
# Open a file containing streamgage ids (in this case gagesII gages)
obs_col_names = ['STAID']
obs_col_types = [np.str_]
obs_cols = dict(zip(obs_col_names, obs_col_types))

df = pd.read_csv('{}/gpr_gagesII_gages_2659.txt'.format(workdir), usecols=obs_col_names, dtype=obs_cols)

# %%
# Build lists for outlet_ids corresponding to gage ids and another for gage ids with no entry in the NHM
outlet_ids = []
missing_stn = []

for ss in df.STAID.values:
    try:
        outlet_ids.append(poi_id_to_seg[ss])
#         print('{}: {}'.format(ss, poi_id_to_seg[ss]))
    except KeyError:
        print('WARNING: {} does not exist in poi_gage_id'.format(ss))
        missing_stn.append(ss)

# %%
# Update the bandit config with the outlet ids
config.update_value('outlets', outlet_ids)

# %%
print(config)

# %%
config.write('{}/bandit.cfg'.format(workdir))

# %%

# %%

# %%
len(missing_stn)

# %%
missing_stn

# %%
len(outlet_ids)

# %%
df.head()
df.to_csv('{}/streamgages.csv'.format(workdir), header=False, index=False)

# %%
aa = datetime.datetime(2000,9,30)
print(aa)

# %%
type(aa)

# %%
isinstance(aa, datetime.datetime)

# %%
str(datetime.datetime.now())

# %%
