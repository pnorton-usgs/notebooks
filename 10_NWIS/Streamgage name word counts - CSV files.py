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
import re
# import xarray as xr
import numpy as np
import pandas as pd

from collections import OrderedDict

# %%
src_path = '/Users/pnorton/Program Development/FY2021_ReFED/for_report/data/trends_wy1960-2019'
filename = f'{src_path}/conus_annual_wy1960-2019_kendall.csv'

# %%
stn_col_names = ['site_no', 'station_nm', 'dec_lat_va', 'dec_long_va',
                 'drain_area_va', 'contrib_drain_area_va',
                 'pval', 'tau', 'trend', 'first_ten_yr', 'last_ten_yr', 'pct_chg']
stn_col_types = [np.str_, np.str_, float, float, float, float,
                 float, float, int, float, float, float]
stn_cols = dict(zip(stn_col_names, stn_col_types))

stations = pd.read_csv(filename, sep='\t', usecols=stn_col_names,
                       dtype=stn_cols)

stations.set_index('site_no', inplace=True)

# %%
stations.head()

# %%
stations['station_nm'] = stations['station_nm'].str.upper().str.title()

# %%
stations.head()

# %%
# df = pd.DataFrame({'A': ['bat', 'foo', 'bait'],
#                    'B': ['abc', 'bar', 'xyz']})
# df.replace(to_replace=r'^ba.$', value='new', regex=True)
#         A    B
# 0   new  abc
# 1   foo  new
# 2  bait  xyz

# %%

# %%

# %%

# %%

# %%
station_names = stations.loc[:, 'station_nm']
station_names.head()

# %%

# %%

# %%
stn_nm_list = stations.loc[:, 'station_nm'].tolist()

# word_list = OrderedDict()
word_list = {}

for xx in stn_nm_list:
    zz = re.split(',| |\(|\)|\"|-|\'', xx)
    # zz = re.split(',| |\(|\)|\.|\"|-|\'', xx)
    for yy in zz:
        lc = yy
        # lc = yy.lower()
        
        if lc not in word_list:
            word_list[lc] = 0
        word_list[lc] += 1


# %%
# del word_list['']
word_list

# %%
hdl = open(f'{src_path}/word_map2.csv', 'w')

for xx in word_list:
    hdl.write(f'{xx},\n')
    
hdl.close()

# %%
words = list(word_list.keys())
words.sort(key=str.lower)
len(words)

# %%
# key=str.lower

# %%

# %%
hdl = open(f'{src_path}/word_map.csv', 'r')

mod_words = {}

for xx in hdl:
    flds = xx.strip().split(',')
    if flds[1] != '':
        mod_words[flds[0]] = flds[1]
    
hdl.close()

# %%
len(mod_words)

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
# word_list = OrderedDict()
word_dict = {}
word_to_id = {}

for ss, xx in df_dict.items():
    words = re.split(',| |\(|\)|\.|\"|-|\'', xx)
    
    for yy in words:
        lc = yy.lower()
        if lc not in word_dict:
            word_dict[lc] = 0
        word_dict[lc] += 1

        if lc not in word_to_id:
            word_to_id[lc] = []
        word_to_id[lc].append(ss)
len(word_dict)

# %%
word_dict

# %%
del word_to_id['']
word_to_id

# %%
df.loc['01AD001']

# %%
df.loc['02365470']

# %%
df.loc['08390500']

# %%
gage_src['poi_name'].loc['06469400'].indexes

# %%
stations.head()

# %%
df = stations
df.drop(columns=['dec_lat_va', 'dec_long_va', 'drain_area_va', 
                 'contrib_drain_area_va', 'pval', 'tau', 'trend', 
                 'first_ten_yr', 'last_ten_yr', 'pct_chg'], inplace=True)
df.head()

# %%
df.to_csv(f'{src_path}/crap.csv', sep='\t', index=True)

# %%
hdl = open(f'{src_path}/word_map.csv', 'r')

mod_words = {}

for xx in hdl:
    flds = xx.strip().split(',')
    if flds[1] != '':
        mod_words[flds[0]] = flds[1]
    
hdl.close()

# %%

# %%

# %%

# %%
import re

def replace(match):
    return mod_words[match.group(0)]



# %%
stn = df.station_nm.tolist()
hdl = open(f'{src_path}/crap3.csv', 'w')

for xx in stn:
    print(f'Before: {xx}')

    yy = re.sub("|".join(fr'\b{re.escape(s)}\b' for s in mod_words), replace, re.sub(r'\.', '', xx))
    hdl.write(f'{yy}\n')
    print(f' After: {yy}')
    print('-'*40)
hdl.close()

# %%

# %%

# %%
df_dict = df.to_dict('index')

# %%
df_dict

# %%
hdl = open(f'{src_path}/crap4.csv', 'w')
hdl.write('site_no\tstation_nm\n')

for site, site_nm in df_dict.items():
    print(f'Before: {site_nm["station_nm"]}')

    yy = re.sub("|".join(fr'\b{re.escape(s)}\b' for s in mod_words), replace, re.sub(r'\.', '', site_nm["station_nm"]))
    hdl.write(f'{site}\t{yy}\n')
    print(f' After: {yy}')
    print('-'*40)
hdl.close()

# %%

# %%
# Read the first-pass cleaned streamgage names

stn_col_names = ['site_no', 'station_nm']
stn_col_types = [np.str_, np.str_]
stn_cols = dict(zip(stn_col_names, stn_col_types))

df_clean1 = pd.read_csv(f'{src_path}/streamgages_clean1.csv', sep='\t', usecols=stn_col_names,
                        dtype=stn_cols)

df_clean1.set_index('site_no', inplace=True)
df_clean1.head()

# %%
stn_clean = pd.merge(stations, df_clean1, how='inner', left_index=True, right_index=True)

# %%
stn_clean.head()

# %%
stn_clean.to_csv(f'{src_path}/streamgages_clean2.csv', sep='\t', index=True, 
                 columns=['station_nm_x', 'station_nm_y', 'dec_lat_va', 'dec_long_va', 
                          'drain_area_va', 'contrib_drain_area_va',
                          'pval', 'tau', 'trend', 'first_ten_yr', 'last_ten_yr', 'pct_chg'])

# %%
