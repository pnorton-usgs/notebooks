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
import xarray as xr

from collections import OrderedDict

# %%
src_path = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/poi_data/*_pois.nc'

# %%
gage_src = xr.open_mfdataset(src_path, chunks={'poi_id': 1040}, combine='nested', concat_dim='poi_id', decode_cf=True, engine='netcdf4')

# %%
gage_src.info()

# %%
df = gage_src['poi_name'].to_pandas()

# %%
df.head()

# %%
df_list = df.tolist()

# %%
# word_list = OrderedDict()
word_list = {}

for xx in df_list:
    zz = re.split(',| |\(|\)|\.|\"|-|\'', xx)
    for yy in zz:
        lc = yy.lower()
        
        if lc not in word_list:
            word_list[lc] = 0
        word_list[lc] += 1


# %%
# del word_list['']
word_list

# %%

# %%
words = list(word_list.keys())
words.sort(key=str.lower)
len(words)

# %%
# key=str.lower

# %%

# %%
words

# %%
df_dict = df.to_dict()
df_dict

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
