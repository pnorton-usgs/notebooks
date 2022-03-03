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
import numpy as np
import pandas as pd
from collections import Counter
from pyPRMS.ParamDb import ParamDb
from Bandit.bandit_multi_locations import read_file

# %%
hw_data_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/headwaters'

paramdb_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS'

falcone_dir = '/Users/pnorton/GIS/gagesII_additiona_data/basinchar_and_report_sept_2011'
# falcone_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHW_obs_sample/helpers'

# %% [markdown]
# ## Load seg_outlets for the headwaters

# %%
hw_df = pd.read_csv(f'{hw_data_dir}/hw_area.csv', sep=',', index_col=0)

# %%
hw_df

# %%
hw_seg_list = hw_df['seg_outlet'].tolist()

# %% [markdown]
# ### Read the full set of seg_ids for headwater extractions

# %%
segs_by_hw = read_file(f'{hw_data_dir}/hw_segs.csv')

# %%
segs_by_hw

# %%
segs_list = []

for kk, vv in segs_by_hw.items():
    # kk is the headwater ID
    # segs_list.append(kk)
    
    for xx in vv:
        segs_list.append(xx)
        
segs_set = set(segs_list)

# %%
print(len(segs_list))
print(len(segs_set))

# %%

# %%
# How many segments occurred more than once?
cnt = Counter(segs_list)
segs_multi_list = [k for k, v in cnt.items() if v == 2]

# %%
len(segs_multi_list)

# %%

# %%

# %% [markdown]
# ## Load parameter database

# %%
pdb = ParamDb(paramdb_dir, verbose=True, verify=True)

# %%
poi_to_seg = pdb.parameters.poi_to_seg
seg_to_poi = {vv: kk for kk, vv in poi_to_seg.items()}

# %%
len(seg_to_poi)

# %%
hw_poi_list = []

for xx in segs_set:
    try:
        hw_poi_list.append(seg_to_poi[xx])
    except KeyError:
        print(f'{xx} has no POI')

# %%
len(hw_poi_list)

# %%
fhdl = open('/Users/pnorton/tmp/nhm_v11_pois.csv', 'w')

for pp in poi_to_seg.keys():
    fhdl.write(f'{pp}\n')
    
fhdl.close()

# %%
col_names = ['STAID', 'CLASS', 'HYDRO_DISTURB_INDX']
col_types = [str, str, int]
cols = dict(zip(col_names, col_types))

falcone_df = pd.read_excel(open(f'{falcone_dir}/gagesII_sept30_2011_conterm.xlsx', 'rb'), sheet_name='Bas_Classif', 
                           usecols=[0, 1, 3], dtype=cols)

falcone_df.info()

# %%
falcone_df.head()

# %%

# %%
falcone_ids = falcone_df['STAID'].tolist()

# %%
falcone_ids

# %%
hw_poi_list

# %%
in_both = set(falcone_ids) & set(hw_poi_list)

# %%
len(in_both)

# %%
