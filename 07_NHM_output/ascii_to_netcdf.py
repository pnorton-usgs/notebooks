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
import pandas as pd
import numpy as np

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/model_output'
# filename = f'{workdir}/nhru_tmaxf_head.csv'

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms5/tests/conus/output'
# filename = f'{workdir}/nhru_out_tmaxc.csv'

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/byHRU_musk/output'
filename = f'{workdir}/nhru_gwres_flow.csv'

# %%
# %%time
df = pd.read_csv(filename, sep=',', skipinitialspace=True, nrows=2, engine='c', memory_map=True,
                 index_col='Date', parse_dates=True, na_values=[-99.0, -999.0])

# df.info()
col_dtypes = {xx: np.float32 for xx in df.columns}

df = pd.read_csv(filename, sep=',', skipinitialspace=True, engine='c', memory_map=True, dtype=col_dtypes, parse_dates=True, index_col='Date', na_values=[-99.0, -999.0])

# %%
# %%time
# Read but don't parse NA values
df = pd.read_csv(filename, sep=',', skipinitialspace=True, nrows=2, engine='c', memory_map=True,
                 index_col='Date', parse_dates=True, keep_default_na=False)

# df.info()
col_dtypes = {xx: np.float32 for xx in df.columns}

df = pd.read_csv(filename, sep=',', skipinitialspace=True, engine='c', memory_map=True, dtype=col_dtypes, parse_dates=True, index_col='Date', keep_default_na=False)

# %%
df.info()

# %%
df.head()

# %%
import read_output_var

# %%
nhru = 109951
nrow = 13242
# nrow = 99

# See: https://xbuba.com/questions/41864984

thedates = np.empty((10, nrow), dtype='c')
var_data = np.zeros((nhru, nrow), dtype=np.float32)
col_names = np.zeros((nhru), dtype=np.int32)

# %%
# %%time

# read_ascii(src, var_dates, var_data, col_names, ncol, nrow)
thedates, var_data, col_names = read_output_var.readfort.read_ascii(filename, thedates, var_data, col_names, nhru, nrow)

# %%
thedates.strides

# %%
var_data.shape

# %%
col_names

# %%

# %%
# See: https://stackoverflow.com/questions/59551264/join-array-of-strings-chars-into-single-string-with-separator-without-using-com

# x = thedates.astype(object)
x = thedates.transpose().astype(object)
x = np.sum(x, axis=1)
# np.sum(x.astype(object), axis = 1)

# x[:,:-1] += ' '
# x.sum(axis=1).reshape(-1, 1)

# %%
x

# %%
var_data.size * var_data.itemsize

# %%
thedates.size * thedates.itemsize

# %%
col_names.size * col_names.itemsize

# %%
b = x[0].decode().split('-')

# %%
import datetime

# c = [int(xx) for xx in b]
# datetime.datetime(*c)

# %%
datetime.datetime(*[int(xx) for xx in x[-1].decode().split('-')])

# %%
zz = [datetime.datetime.strptime(yy.decode(), '%Y-%m-%d') for yy in x]
print(zz)
# for yy in x:
#     zz = datetime.datetime.strptime(yy.decode(), '%Y-%m-%d')
#     print(zz)

# %%
# def get_datetime(time_str):
#     return datetime.strptime(time_str.tostring(), '%Y-%m-%d')

# timelist = [get_datetime(xx) for xx in time_str]

# %%
x

# %%
x[-1]

# %%
list(range(1,14+1))

# %%
import math

math.frexp(math.pi)

# %%
