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
#     display_name: Python [conda env:py27]
#     language: python
#     name: conda-env-py27-py
# ---

# %%
import pandas as pd
import numpy as np

# %%
df = pd.read_csv('/Users/pnorton/tmp/prms_params_raw.txt', sep='\t') #, quotechar="'")

# hdl = open('/Users/pnorton/tmp/prms_params_raw.txt', 'r')

# %%
df.head()

# %%
aa = df.iloc[:,1].tolist()
aa.sort()

# %%
len(aa)

# %%
len(set(aa))

# %%
