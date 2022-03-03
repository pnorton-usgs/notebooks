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
import numpy as np
import datetime

# %%

# %%
workdir = '/Users/pnorton/USGS/Projects/National_Hydrology_Model/regions/r10U/error_bounds/SWE'
st = datetime.datetime(2003,10,31)
en = datetime.datetime(2010,9,30)

# %%
file1 = 'SWE_SNODAS_Daily_r10U_byHRU_2003-09-30_2014-06-13.csv'
ds1 = pd.read_csv('%s/%s' % (workdir, file1), na_values=['MEAN(m)'], header=True)
ds1.drop(0, inplace=True)

ds1.rename(columns={ds1.columns[0]:'thedate'}, inplace=True)
ds1['thedate'] = pd.to_datetime(ds1['thedate'])
ds1.set_index('thedate', inplace=True)
#ds1.head()

# %%
ds1_mth = ds1.resample('M', how='mean')
ds1_mth = ds1_mth[st:en]

# Convert meters to inches
ds1_mth = ds1_mth * 39.3701
ds1_mth.head()

# %%
file2 = 'wb.swe'
ds2 = pd.read_csv('%s/%s' % (workdir, file2), sep=' ', skipinitialspace=True, 
                  parse_dates={'thedate': [0, 1]}, index_col=['thedate'])

# The wb.swe file has messed up column headers for r10u, so renumber them.
ds2.rename(columns=lambda x: ds2.columns.get_loc(x)+1, inplace=True)

# Pandas messes up the day of the month when parsing dates using just a year and month.
# Resampling to month again will fix the dates without altering the values.
ds2 = ds2.resample('M', how='mean')
ds2 = ds2[st:en]

# Convert from millimeters to inches
ds2 = ds2 * 0.0393701

# %%

# %%
modis = pd.DataFrame(ds1_mth.ix[:,0])
modis.rename(columns={modis.columns[0]: 'modis'}, inplace=True)
#print modis.head()

wb = pd.DataFrame(ds2.ix[:,1])
wb.rename(columns={wb.columns[0]: 'wb'}, inplace=True)
#print wb.head()

ds_join = modis.join(wb)

ds_join['min'] = ds_join.min(axis=1)
ds_join['max'] = ds_join.max(axis=1)
ds_join.drop(['modis', 'wb'], axis=1, inplace=True)
ds_join['year'] = ds_join.index.year
ds_join['month'] = ds_join.index.month
ds_join.reset_index(inplace=True)
ds_join.drop(['thedate'], axis=1, inplace=True)
ds_join.to_csv('crap.swe', sep=' ', float_format='%0.5f', columns=['year', 'month', 'min', 'max'], 
               header=False, index=False)
print ds_join

# %%
