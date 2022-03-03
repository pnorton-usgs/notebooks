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
#     display_name: Python [conda env:py27]
#     language: python
#     name: conda-env-py27-py
# ---

# %%
import pandas as pd
import numpy as np
import datetime

# %%
srcdir = "/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v3/2017-09-18_daymet_v3"
dstdir = "/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v3/merged"

# %%
REGIONS = ['r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09',
           'r10L', 'r10U', 'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18']

vars = ['tmin', 'tmax', 'prcp']
# region = 'r01'
# var = vars[0]

# tmin_Daymet_daily_r09_byHRU_2001-01-01_2016-12-31.csv
date_ranges = ['1980-01-01_2000-12-31', '2001-01-01_2016-12-31']

# %%
for region in REGIONS:
    print(region)
    for var in vars:
        print('\t{}'.format(var))
        
        df_a = pd.read_csv('{}/{}_Daymet_daily_{}_byHRU_{}.csv'.format(srcdir, var, region, 
                                                                       date_ranges[0]), header=0, skiprows=[0, 2])

        df_a.rename(columns={df_a.columns[0]: 'time'}, inplace=True)
        df_a['time'] = pd.to_datetime(df_a['time'])
        df_a.set_index('time', inplace=True)
        # df_a.head()

        df_b = pd.read_csv('{}/{}_Daymet_daily_{}_byHRU_{}.csv'.format(srcdir, var, region, 
                                                                       date_ranges[1]), header=0, skiprows=[0, 2])

        df_b.rename(columns={df_b.columns[0]: 'time'}, inplace=True)
        df_b['time'] = pd.to_datetime(df_b['time'])
        df_b.set_index('time', inplace=True)

        # Drop first row in second file to remove duplicate date
        df_b.drop(df_b.index[0], inplace=True)
        # df_b.head()

        # Concatenate the two dataframes
        df_c = pd.concat([df_a, df_b])
        # print(df_c.head())

        # Reindex to the full date range we want. This fills in with NaNs
        # when a missing date is added.
        df2 = df_c.resample('D').mean()

        dr = pd.date_range('01-01-1980', '12-31-2016')
        df2 = df2.reindex(dr)
        # print(df2.head())

        # Interpolate values for the missing values (e.g. 12-31 for leap years)
        df3 = df2.interpolate(method='time')

        # Output the merged dataframe
        df3.to_csv('{}/{}_Daymet_daily_{}_byHRU_1980-01-01_2016-12-31.csv'.format(dstdir, var, region), 
                   index_label=['time']) 

# %%
df_c = pd.concat([df_a, df_b])
df_c.head()

# %%
dr = pd.date_range('01-01-1980', '12-31-2016')
df2 = df2.reindex(dr)

# df2 = df_c.resample('D').mean()
df2.head()

# %%
df2.info()

# %%
# Interpolate values for the missing values (e.g. 12-31 for leap years)
df3 = df2.interpolate(method='time')
# df3.iloc[365]

# %%
df3.info()

# %%
df3.tail()

# %%
df3.to_csv('{}/{}_Daymet_daily_{}_byHRU_1980-01-01_2016-12-31.csv'.format(dstdir, var, region)) 
                                                            

# %%
aa = None
if not aa:
    print('aa is False')

# %%
