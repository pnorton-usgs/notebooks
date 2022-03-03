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
#     display_name: Python 2
#     language: python
#     name: python2
# ---

# %%
import pandas as pd
import numpy as np

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/files/tmp'
filename = 'SDCmatrix_wDEFAULTS'

# %%

df = pd.read_csv('{}/{}'.format(workdir, filename), sep=r'\s*', engine='python', header=None,
                 names=['nhmid', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11'],
                 index_col=['nhmid'])

# %%
df.info()

# %%
df.tail()

# %%
bb = df.unstack().unstack()
bb.head()

# %%
bb[57864]

# %%
nhmid = [57864, 57865, 57868, 57869, 57870,
         57873, 57874, 57875, 57878, 57879,
         57880, 57881, 57882, 57883]

pd.set_option('precision',4)

for xx in nhmid:
    for yy in bb[xx].values:
        print yy

# %%
