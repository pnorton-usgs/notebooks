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
#     display_name: Python [default]
#     language: python
#     name: python2
# ---

# %%
from __future__ import print_function
import pandas as pd
import numpy as np

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/Study_areas/Alaska/Solar Radiation'

# %%
df = pd.read_csv('{}/crb_daymet_swrad.csv'.format(workdir), comment='#', skiprows=[2], 
                 index_col=[0])
#                  parse_dates=[0], index_col=[0])

# %%
df.head()

# %%
df_langley = df * 2.0636285
df_langley.head()

# %%
df_langley.to_csv('{}/crb_daymet_solrad_langley.csv'.format(workdir), index=True, header=True)

# %%
df_ann = df.resample('A-DEC').mean()
df_ann.head()

# %%
df_ann.max().max()

# %%
