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

# %%
df = pd.read_csv('/Users/pnorton/tmp/ts14.parval', header=None, sep=' ', skipinitialspace=True, skiprows=1,
                 usecols=[0,1], index_col=0)
print df.info()
print df.head()

# %%
try:
    print df.loc['RCH1_TS1'].values[0]
except KeyError:
    print 'using default'

# %%
df_riv = pd.read_csv('/Users/pnorton/tmp/ts1aFIXED.riv', header=None, sep=' ', skipinitialspace=True, skiprows=2,
                 usecols=[0,1,2,3,4,5])
print df_riv.head()
print df_riv.info()

# %%
