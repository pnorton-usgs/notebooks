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
from pyPRMS.ParamDb import ParamDb

# %% [markdown]
# ### Calculate the total area for all headwaters in NHMv1.1

# %%

# %%
base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/headwaters'
hw_area_file = f'{base_dir}/hw_area.csv'

# %%
hwarea_df = pd.read_csv(hw_area_file, sep=',', header=0)
hwarea_df

# %%
hwarea_df.info()

# %%
total_area = hwarea_df['seg_cum_area_sqkm'].sum()

# %%
total_area / 2.59

# %%
total_area * 0.38610

# %%
