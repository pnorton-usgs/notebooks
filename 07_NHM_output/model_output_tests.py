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
import pandas as pd
import numpy as np

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/NHM_output'
filename = '{}/nhru_hru_snow.zip'.format(workdir)

# %%

# df = pd.read_csv(cbh_file, sep=' ', skipinitialspace=True,
#                  skiprows=3, engine='c', memory_map=True,
#                  date_parser=dparse, parse_dates={'time': CBH_INDEX_COLS},
#                  index_col='time', header=None, na_values=[-99.0, -999.0, 'NaN', 'inf'])


df = pd.read_csv(filename, sep=',', engine='c', memory_map=True, header=1)

# %%
fhdl = open('someparam.csv', 'w')
fhdl.write('id,<param>\n')

for vars in stats:
    fhdl.write
