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
import datetime
import xarray as xr
import pandas as pd
import numpy as np

from pyPRMS.ControlFile import ControlFile
from Bandit.bandit_multi_locations import read_file
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ParamDb import ParamDb
from pyPRMS.Streamflow import Streamflow

# %%
workdir = '/Volumes/USGS_NHM1/bandit_jobs/20210408_v11_gridmet_headwaters/hw_0833'

ctl_filename = f'{workdir}/../control.default.bandit'

# %%
ctl = ControlFile(ctl_filename)

# %%
obs_filename = f'{workdir}/{ctl.get("data_file").values}'
param_filename = f'{workdir}/{ctl.get("param_file").values}'


# %%
def read_streamflow(filename):
    fhdl = open(filename, 'r')
    rawdata = fhdl.read().splitlines()
    fhdl.close()

    it = iter(rawdata)

    hdr_count = 0

    # Look for start of ID section
    for line in it:
        hdr_count += 1
        if line == '// ID':
            # print('Found ID line')
            break

    # Read the gage IDs
    gage_ids = []
    for line in it:
        flds = line.split(' ')
        hdr_count += 1

        if flds[0][0:4] == '////':
            # print('end of gage ids')
            break

        gage_ids.append(flds[1])

    # Skip the read of the header
    for line in it:
        hdr_count += 1
        if line[0:4] == '####':
            # print('end of header section')
            break


    # Read the streamflow data
    col_names = ['year', 'month', 'day', 'hour', 'min', 'sec']

    for ii in gage_ids:
        col_names.append(ii)

    df = pd.read_csv(filename, sep='\s+', header=None, names=col_names, 
                     skiprows=hdr_count, na_values=[-999.0],
                     parse_dates={'time': ['year', 'month', 'day']},
                     index_col='time')

    df.drop(columns=['hour', 'min', 'sec'], inplace=True)
    
    return df

# %%
df = read_streamflow(obs_filename)
df.info()

# %%
df.head()

# %%
df.notna().sum(axis=0)

# %%
# Reduce dataframe POIs which have more than 1 year of valid data
df.loc[:, df.notna().sum(axis=0) > 365]

# %%

# %%
st_date = datetime.datetime(*(ctl.get('start_time').values))
en_date = datetime.datetime(*(ctl.get('end_time').values))

# Restrict observations to model date range
df = df[st_date:en_date]

# %%

# %%

# %%
# Get the POI IDs which have 1 year or less of observation values
bad_poi_list = df.loc[:, df.notna().sum(axis=0) <= 365].columns.tolist()
bad_poi_list

# %%

# %%
params = ParameterFile(param_filename, verbose=True)

# %%
poi_ids = params.parameters.get('poi_gage_id')

# %%
poi_ids.data

# %%
# poi_list = ['02157000', '02160105']
params.remove_poi(bad_poi_list)

# %%
poi_ids.data

# %%
print(params.dimensions)

# %%
poi_to_seg = params.parameters.poi_to_seg0
poi_to_seg

# %%
# seg_cum_area = params.parameters.get_dataframe('seg_cum_area')
seg_cum_area = params.parameters.get('seg_cum_area').data
seg_cum_area

# %%
seg_cum_area[list(poi_to_seg.values())] / 640.0

# %%

# %%
selected_pois = ['02157470', '02160105', '02157490', '02157510']

df.loc[:, ['02157470', '02160105', '02157490', '02157510']]

# %%
df.index.second

# %%
ctl.get('start_time').values

# %%
ctl.get('start_time').values=[1980, 10, 2]

# %%
ctl.get('start_time').values

# %%
a_date = datetime.datetime(1980, 10, 3)
a_date

# %%
ctl.get('start_time').values = [a_date.year, a_date.month, a_date.day, 0, 0, 0]

# %%
ctl.get('start_time').values

# %%
import os
somedir = '/Volumes/USGS/something/somefile.nc'
dirpart = os.path.dirname(somedir)
nameparts = os.path.splitext(os.path.basename(somedir))
print(dirpart)
print(nameparts)

# %%
os.path.basename(somedir)

# %%
params.dimensions['npoigages'].size

# %%
os.listdir('/Volumes/USGS_NHM1/bandit_jobs/20210408_v11_gridmet_headwaters')

# %%
thedir = '/Volumes/USGS_NHM1/bandit_jobs/20210408_v11_gridmet_headwaters'
[name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name))]

# %%
import glob

headwater_dirs = glob.glob(f'/Volumes/USGS_NHM1/bandit_jobs/20210408_v11_gridmet_headwaters/hw_*')

# %%
just_hw = [os.path.basename(xx) for xx in headwater_dirs]

# %%
just_hw

# %%
just_hw.sort()

# %%
just_hw

# %%
odd_years = list(range(1979, 2019+1, 2))
odd_years

# %%
df_odd = df.copy()
df_odd['o_year'] = df_odd.index.year % 2
df_odd = df_odd[df_odd['o_year'] == 1]

# %%
st = datetime.datetime(1979, 1, 1)
en = datetime.datetime(1984, 12, 31)
df_odd[st:en]

# %%
df_odd.loc[:, df_odd.notna().sum(axis=0) <= 365].columns.tolist()

# %%
