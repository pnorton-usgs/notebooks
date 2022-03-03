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
def dparse(yr, mo, dy, hr, mn, sc):
    yr, mo, dy, hr, mn, sc = [int(x) for x in [yr, mo, dy, hr, mn, sc]]
    dt = datetime.datetime(yr, mo, dy, hr, mn, sc)
#     dt = datetime.datetime(int(yr), int(mo), int(dy), int(hr), hr(mn), int(sc))
    return dt


# %%
df = pd.read_csv('../from_Eddie/br.data', skiprows=12, skipinitialspace=True, header=None, sep=' ',
                 date_parser=dparse, parse_dates={'thedate': [0,1,2,3,4,5]})

# %%
df.shape

# %%
df.head()

# %%
df.tail()

# %%
# df.isnull()
df[df.isnull().any(axis=1)].shape[0]

# %%
df.info()

# %%
ff = open('../from_Eddie/br_data_tmin_trunc.txt', 'r')

ll = []

for lines in ff:
    ll.append(len(lines))

print max(ll)

# %%
# df = pd.read_csv('../from_Eddie/br_data_tmin.txt', skiprows=11, skipinitialspace=True, header=None, sep=' ',
#                  date_parser=dparse, parse_dates={'thedate': [0,1,2,3,4,5]})
df = pd.read_csv('../from_Eddie/br_data_tmax.txt', skiprows=11, skipinitialspace=True, header=None, sep=' ')


# %%
df.head()

# %%
df.to_csv('../from_Eddie/br_data_tmax_trunc.txt', float_format='%d', header=False, index=False)

# %%
