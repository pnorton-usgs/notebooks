# ---
# jupyter:
#   jupytext:
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

# %%
import pandas as pd

# %%
base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11'

# %%
# paramdb_v11_gridmet_CONUS

# %%
# Read in all HRUs
fhdl = open(f'{base_dir}/paramdb_v11_gridmet_CONUS/nhm_id.csv', 'r', encoding='ascii')
rawdata = fhdl.read().splitlines()
fhdl.close()

it = iter(rawdata)
next(it)   # Skip header

hru_list = {}
for row in it:
    flds = row.split(',')
    hru_list[int(flds[1])] = dict(byHRU=1)

# %%

# %%
# Read in the byHW HRUs
fhdl = open(f'{base_dir}/headwaters/hw_hrus.csv', 'r', encoding='ascii')
rawdata = fhdl.read().splitlines()
fhdl.close()

it = iter(rawdata)
next(it)   # Skip header

byHW_hru_list = {}
for row in it:
    flds = row.split(',')
    hru_id = int(flds[0])
    hw_id = int(flds[1])
    
    if hru_id in hru_list:
        hru_list[hru_id]['byHW'] = 1
        if hw_id not in byHW_hru_list:
            byHW_hru_list[hw_id] = []
        byHW_hru_list[hw_id].append(hru_id)

# %%
# Read in the byHW HRUs
fhdl = open(f'{base_dir}/headwaters/HW_obs/nhm_v11_byHWobs_hrus.csv', 'r', encoding='ascii')
rawdata = fhdl.read().splitlines()
fhdl.close()

it = iter(rawdata)
next(it)   # Skip header

byHWobs_hru_list = {}
for row in it:
    flds = row.split('\t')
    hru_id = int(flds[0])
    hw_id = int(flds[1])
    
    if hru_id in hru_list:
        hru_list[hru_id]['byHWobs'] = 1
        if hw_id not in byHWobs_hru_list:
            byHWobs_hru_list[hw_id] = []
        byHWobs_hru_list[hw_id].append(hru_id)

# %%
# %%time
df = pd.DataFrame(hru_list, dtype='Int64').transpose()
df.byHRU = df.byHRU.fillna(0)
df.byHW = df.byHW.fillna(0)
df.byHWobs = df.byHWobs.fillna(0)
df['level'] = df.byHRU + df.byHW + df.byHWobs

# %%
df.head()

# %%
# Save results to csv
df.to_csv(f'{base_dir}/headwaters/nhm_v11_calibration_levels.csv', sep=',', index_label='nhm_id')

# %%

# %%
