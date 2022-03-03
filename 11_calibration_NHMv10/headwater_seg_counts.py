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

# %%

# %%
import numpy as np
import pandas as pd
from collections import Counter

from Bandit.bandit_multi_locations import read_file

# %%
hw_data_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/headwaters_conus'

hw_file = f'{hw_data_dir}/2018-12-17_hwSegsVll.csv'

# %%
segs_by_hw = read_file(hw_file)

# %%
segs_list = []

for kk, vv in segs_by_hw.items():
    # kk is the headwater ID
    # segs_list.append(kk)
    
    for xx in vv:
        segs_list.append(xx)
        
segs_set = set(segs_list)

print(len(segs_list))
print(len(segs_set))

# %%
# How many segments occurred more than once?
cnt = Counter(segs_list)
segs_multi_list = [k for k, v in cnt.items() if v > 1]

print(len(segs_multi_list))

# %%
segs_multi_list

# %%
segs_by_hw

# %%
