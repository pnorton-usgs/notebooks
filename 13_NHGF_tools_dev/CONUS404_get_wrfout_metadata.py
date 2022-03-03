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
#     display_name: Python [conda env:holoview]
#     language: python
#     name: conda-env-holoview-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
from dask.distributed import Client
# import geopandas
import numpy as np
import pandas as pd
# import os
import xarray as xr
from collections import Counter

# %%
# client = Client(n_workers=2, threads_per_worker=2, memory_limit='1GB')
# client = Client(processes=False)
client = Client()
client

# %%
client.close()

# %%
work_dir = '/Volumes/USGS_NHM2/datasets/wrf-conus404/originals'

# Variables that are integrated over 60 minutes per hourly timestep
vars_60min_accum = ['ACDEWC', 'ACDRIPR', 'ACDRIPS', 'ACECAN', 'ACEDIR', 'ACETLSM', 'ACETRAN',
                    'ACEVAC', 'ACEVB', 'ACEVC', 'ACEVG', 'ACFROC', 'ACFRZC', 'ACGHB', 'ACGHFLSM',
                    'ACGHV', 'ACINTR', 'ACINTS', 'ACIRB', 'ACIRC', 'ACIRG', 'ACLHFLSM', 'ACLWDNLSM',
                    'ACLWUPLSM', 'ACMELTC', 'ACPAHB', 'ACPAHG', 'ACPAHLSM', 'ACPAHV', 'ACPONDING',
                    'ACQLAT', 'ACQRF', 'ACRAINLSM', 'ACRAINSNOW', 'ACRUNSB', 'ACRUNSF', 'ACSAGB',
                    'ACSAGV', 'ACSAV', 'ACSHB', 'ACSHC', 'ACSHFLSM', 'ACSHG', 'ACSNBOT', 'ACSNFRO',
                    'ACSNOWLSM', 'ACSNSUB', 'ACSUBC', 'ACSWDNLSM', 'ACSWUPLSM', 'ACTHROR', 'ACTHROS',
                    'ACTR', 'GRAUPEL_ACC_NC', 'PREC_ACC_NC', 'SNOW_ACC_NC']

# Variables that are accumulated from model start
vars_model_accum = ['ACLWDNB', 'ACLWDNBC', 'ACLWDNT', 'ACLWDNTC', 'ACLWUPB', 'ACLWUPBC',
                    'ACLWUPT', 'ACLWUPTC', 'ACSNOM', 'ACSWDNB', 'ACSWDNBC', 'ACSWDNT',
                    'ACSWDNTC', 'ACSWUPB', 'ACSWUPBC', 'ACSWUPT', 'ACSWUPTC']


# %%
# Read word map file for processing the description strings
fhdl = open('wrfout_words.txt', 'r', encoding='ascii')
rawdata = fhdl.read().splitlines()
fhdl.close()

it = iter(rawdata)
next(it)   # Skip header

word_map = {}
for row in it:
    flds = row.split('\t')
    if len(flds[2]) != 0:
        word_map[flds[0].replace('"', '')] = flds[2].replace('"', '')
    print(flds)

# %%

# %%

# %%

# %%
df = xr.open_dataset(f'{work_dir}/wrfout_d01_2020-09-30_00:00:00', decode_coords=False, chunks={})

# %%
df

# %%
# dim_names = [dd for dd in df.dims.keys()]
# dim_names

# dim_cnt = Counter(dim_names)
# dim_cnt = Counter()
attr_cnt = Counter()
word_cnt = Counter()
# dim_cnt['Time'] += 1

wrfout_vars = {}

for vv in list(df.keys()):
    cvar = df[vv]
    wrfout_vars[vv] = {}
    
    for cattr, val in cvar.attrs.items():
        if cattr in ['description', 'units', 'coordinates']:
            attr_cnt[cattr] += 1
            wrfout_vars[vv][cattr] = val
            
            if cattr == 'description':
                new_val = []
                for ww in val.split(' '):
                    if ww in word_map:
                        new_val.append(word_map[ww])
                    else:
                        new_val.append(ww)
                    word_cnt[ww] += 1
                    
#                 result = string[0].upper() + string[1:]
                outstr = ' '.join(new_val)
    
                if len(outstr) > 0:
                    outstr = outstr[0].upper() + outstr[1:]
                wrfout_vars[vv]['description_new'] = outstr
    
    wrfout_vars[vv]['datatype'] = cvar.encoding['dtype'].name
    wrfout_vars[vv]['dimensions'] = ' '.join(cvar.dims)
    
    if vv in vars_60min_accum:
        # Add a cell_methods field
        wrfout_vars[vv]['cell_methods'] = 'XTIME: sum (interval: 1 minute)'
    
#         for kk in list(const_df.keys()):
#             src_cvar = const_df[kk]

#             # Create list of dimensions, modify dim names as needed
#             cvar_dims = []
#             cvar_cnk = []
#             for xx in src_cvar.dims:

# %%
wrfout_vars

# %%
attr_cnt

# %%
out_df = pd.DataFrame(wrfout_vars).transpose()
out_df.head()

# %%
out_df.sort_index().to_csv('wrfout_metadata.txt', sep='\t', columns=['description_new', 'description', 'units', 'cell_methods', 'dimensions', 'coordinates', 'datatype'])

# %%
word_df = pd.DataFrame(word_cnt, index=[0]).transpose()
word_df.head()

# %%
word_df.to_csv('wrfout_words.csv', sep='\t')

# %%

# %%
fhdl = open('wrfout_words.txt', 'r', encoding='ascii')
rawdata = fhdl.read().splitlines()
fhdl.close()

it = iter(rawdata)
next(it)   # Skip header

word_map = {}
for row in it:
    flds = row.split('\t')
    if len(flds[2]) != 0:
        word_map[flds[0].replace('"', '')] = flds[2].replace('"', '')
    print(flds)
    

# %%
word_map['LATITUDE,']

# %%
len(flds[1])

# %%

# %%

# %%

# %%

# %%

# %%
