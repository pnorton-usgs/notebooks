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
import dask.dataframe as dd
import numpy as np
import pandas as pd
from dask.distributed import Client
import dask

# %%
client = Client()

# %%
client

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms5/tests/conus/output'
filename = f'{workdir}/nhru_out_tmaxc.csv'

# %%
# %%time
# Read but don't parse NA values
df = pd.read_csv(filename, sep=',', skipinitialspace=True, nrows=2, engine='c', memory_map=True,
                 index_col='Date', parse_dates=True, keep_default_na=False)

# df.info()
col_dtypes = {xx: np.float32 for xx in df.columns}

df = dd.read_csv(filename, sep=',', skipinitialspace=True, engine='c', memory_map=True, dtype=col_dtypes, 
                 parse_dates=True, index_col='Date', keep_default_na=False)

# %%
