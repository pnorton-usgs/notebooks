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

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import datetime
import msgpack
import pandas as pd

from pyPRMS.Cbh import Cbh
from pyPRMS.prms_helpers import dparse

CBH_VARNAMES = ['prcp', 'tmin', 'tmax']
CBH_INDEX_COLS = [0, 1, 2, 3, 4, 5]
REGIONS = ['r01', 'r02', 'r03', 'r04']

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet'

cbh_file = '{}/r01_prcp.cbh.gz'.format(workdir)

st_date = datetime.datetime(1980, 1, 1)
en_date = datetime.datetime(1980, 2, 1)

idx_retrieve = {1: 1, 2: 2}

cc1 = Cbh(filename=cbh_file, indices=idx_retrieve, st_date=st_date, en_date=en_date)
# cc1 = Cbh(filename=cbh_file, st_date=st_date, en_date=en_date)
cc1.read_cbh()

# %%
cc1.data.columns.values.tolist()

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet/tmp'

infile = open('{}/r10L_prcp.cbh'.format(workdir), 'r')

outfile = open('{}/crap.cbh'.format(workdir), 'wb')

for xx in infile:
    outfile.write(xx)
infile.close()
outfile.close()

# %%
cc1.data.head()


# %%
def get_parameter(filename):
    with open(filename, 'rb') as ff:
        return msgpack.load(ff, use_list=False)

st_date = datetime.datetime(1980, 1, 1)
en_date = datetime.datetime(1980, 2, 1)

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet'
merged_paramdb_dir = '/Users/pnorton/Projects/National_Hydrology_Model/paramDb/merged_params_3'
var = 'prcp'

CBH_VARNAMES = ['prcp', 'tmin', 'tmax']
CBH_INDEX_COLS = [0, 1, 2, 3, 4, 5]
REGIONS = ['r01', 'r02', 'r03', 'r04']

hru_nhm_to_local = get_parameter('{}/hru_nhm_to_local.msgpack'.format(merged_paramdb_dir))
hru_nhm_to_region = get_parameter('{}/hru_nhm_to_region.msgpack'.format(merged_paramdb_dir))

cbh_hdl = Cbh(indices=hru_nhm_to_local, mapping=hru_nhm_to_region, var='prcp', st_date=st_date, en_date=en_date, 
              regions=['r02'])

cbh_hdl.read_cbh_multifile(workdir)
cbh_hdl.data.axes[1].rename('hru', inplace=True)



# %%
print(CBH_INDEX_COLS)
idx_retrieve = {1: 500, 2: 20}
aa = list(CBH_INDEX_COLS)
aa.extend(idx_retrieve)
print(aa)
print(CBH_INDEX_COLS)

# %%
print(hru_nhm_to_region)

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet/tmp'
cbh_file = '{}/r17_prcp.cbh'.format(workdir)
             
idx_retrieve = [6886, 6887, 6101, 6124, 6119, 6077, 6118, 
                6029, 6033, 6078, 6079, 6056, 6066, 6024, 
                6057, 6063, 5906, 5907, 5824, 5831, 5897, 5899]

load_cols = list(CBH_INDEX_COLS)
load_cols.extend([xx+5 for xx in idx_retrieve])

st_date = datetime.datetime(1980, 1, 1)
en_date = datetime.datetime(1980, 2, 1)

# Columns 0-5 always represent date/time information
df = pd.read_csv(cbh_file, sep=' ', skipinitialspace=True,
                 usecols=load_cols,
                 skiprows=3, engine='c', memory_map=True,
                 date_parser=dparse, parse_dates={'time': CBH_INDEX_COLS},
                 index_col='time', header=None, na_values=[-99.0, -999.0])


# %%
st_date = datetime.datetime(1980, 10, 1)
en_date = datetime.datetime(1980, 11, 1)

df2 = df[st_date:en_date]
df2.head(15)



# %%
df2.loc[:,6891]

# %%
gdp_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet/gdp'
gdp_file = '{}/prcp_DAYMET_Daily_r17_byHRU_1980-01-01_2014-12-31.csv'.format(gdp_dir)

gdp_data = pd.read_csv(gdp_file, na_values=[255.], header=0, skiprows=[0, 2])

gdp_data.rename(columns={gdp_data.columns[0]: 'time'}, inplace=True)
gdp_data['time'] = pd.to_datetime(gdp_data['time'])
gdp_data.set_index('time', inplace=True)

# %%
gdp_data.head(10)

# %%
gdp_data.info()

# %%
df2 = gdp_data[st_date:en_date]
df2 *= 0.0393701
df2.loc[:,'6886']

# %%
aa = [int(xx) for xx in gdp_data.columns.tolist()]

# %%
for xx in range(len(aa)-1):
    if aa[xx+1] - aa[xx] != 1:
        print('jump in sequence at {}'.format(xx))

# %%
