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

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import datetime
import numpy as np
import os
import pandas as pd

from pyPRMS.Cbh import Cbh

# %%

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/prms5/test4'
var = 'tmax'
filename = '{}/{}.cbh'.format(workdir, var)
filename_out = '{}/{}_c.cbh'.format(workdir, var)

# %%
df1 = Cbh(filename=filename)
df1.read_cbh_full()
df1b = df1.data

# Create list of column order for output
out_order = range(1, len(df1b.columns)+1)
out_order = [k+5 for k in out_order]

# Add additional column names to output order
for cc in ['second', 'minute', 'hour', 'day', 'month', 'year']:
    out_order.insert(0, cc)
print(out_order)

# %%

# %%
# Convert from Fahrenheit to Celsius
df2 = (df1b - 32.0) / 1.8

# Add additional columns required for CBH file format
df2['year'] = df2.index.year
df2['month'] = df2.index.month
df2['day'] = df2.index.day
df2['hour'] = 0
df2['minute'] = 0
df2['second'] = 0

df2.head()

# %%
out_cbh = open(filename_out, 'w')
out_cbh.write('Written by Bandit\n')
out_cbh.write('{} {}\n'.format(var, len(df1b.columns)))
out_cbh.write('########################################\n')
df2.to_csv(out_cbh, columns=out_order, sep=' ', index=False, header=False, encoding=None, chunksize=50)
out_cbh.close()

# %%

# %%

# %%
