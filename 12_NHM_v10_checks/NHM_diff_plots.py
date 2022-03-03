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
# %matplotlib inline
import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as dates
from matplotlib.dates import DayLocator, MonthLocator, YearLocator

import datetime
import os
import xarray as xr
import numpy as np
import pandas as pd

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/parameters'


# %%
AB_df = xr.open_dataset(f'{workdir}/A_minus_B.nc', mask_and_scale=False)
AC_df = xr.open_dataset(f'{workdir}/A_minus_C.nc', mask_and_scale=False)
AD_df = xr.open_dataset(f'{workdir}/A_minus_D.nc', mask_and_scale=False)

BC_df = xr.open_dataset(f'{workdir}/B_minus_C.nc', mask_and_scale=False)
BD_df = xr.open_dataset(f'{workdir}/B_minus_D.nc', mask_and_scale=False)

CD_df = xr.open_dataset(f'{workdir}/C_minus_D.nc', mask_and_scale=False)

CLH_df = xr.open_dataset(f'{workdir}/C_minus_lauren_byHWmusk.nc', mask_and_scale=False)

BLH_df = xr.open_dataset(f'{workdir}/B_minus_lauren_byHRU.nc', mask_and_scale=False)

LH_test_df = xr.open_dataset(f'{workdir}/LHbyHRU_minus_test.nc', mask_and_scale=False)

# %%
AB_df

# %%
AB_df['wrain_intcp'].plot()


# %%
def plot_param_3df(varname, df1, df2, df3, df1_name, df2_name, df3_name):
#     varname = 'jh_coef'
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20, 4))
    ax = axes.flatten()

    df1[varname].plot(ax=ax[0])
    df2[varname].plot(ax=ax[1])
    df3[varname].plot(ax=ax[2])

    ax[0].set_title(df1_name, fontsize=12)
    ax[1].set_title(df2_name, fontsize=12)
    ax[2].set_title(df3_name, fontsize=12)


# %%
plot_param_3df('carea_max', AB_df, AC_df, BLH_df, 'master-byHRU', 'master-musk', 'byHRU-lauren_byHRU')

# %%
# a = np.ma.masked_where(a < 0.05, a)

# cmap = plt.cm.OrRd
# cmap.set_bad(color='black')

aa = np.ma.masked_where(AC_df['adjmix_rain'].values == 0.0, AC_df['adjmix_rain'].values)
aa['adjmix_rain'].plot()

# %%
plot_param_3df('fastcoef_lin', AB_df, AC_df, AD_df)

# %%
np.all(AB_df['jh_coef'].values == 0)


# %%
def get_non_zero(df):
    vars = list(df.keys())
    nz = []
    
    for vv in vars:
        all_zero = np.all(df[vv].values == 0)

        if not all_zero:
            nz.append(vv)
            #print(f'{vv}: {all_zero}')
    return nz


# %%
AB_nz = get_non_zero(AB_df)
AC_nz = get_non_zero(AC_df)
AD_nz = get_non_zero(AD_df)
BC_nz = get_non_zero(BC_df)
BD_nz = get_non_zero(BD_df)
CD_nz = get_non_zero(CD_df)
LH_test_nz = get_non_zero(LH_test_df)

# %%
# for xx, yy in AB_df.variables.items():
#     print(str(xx))
    
#     print(f'  {AB_df[xx].encoding["dtype"]}')
#     for zz in yy.attrs:
#         print(f'  {zz}')
    
#     for zz in yy.dims:
#         print(f'  {zz}')
aa = list(set(AB_nz).union(set(AC_nz).union(set(AD_nz).union(set(BC_nz).union(set(BD_nz).union(CD_nz))))))
aa.sort()
aa

# %%
AB_nz

# %%
for vv in AB_nz:
    print(vv)
    if vv not in ['poi_gage_id', 'K_coef', 'snarea_curve']:
        plot_param_3df(vv, AB_df, AC_df, AD_df, 'master-byHRU', 'master-musk', 'master-musk_obs')

# %%
CD_nz = get_non_zero(CD_df)

# %%
CD_nz

# %%
lim = 0.00000000001
CD_df['K_coef'].plot(ylim=[-lim, lim], label='K_coef')


# %%
def plot_param(varname, df):
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
#     ax = axes.flatten()

    df[varname].plot(ax=axes)
    axes.set_title(f'{varname}', fontsize=12)

for vv in CD_nz:
#     print(vv)
    if vv not in ['poi_gage_id', 'K_coef']:
        plot_param(vv, CD_df)


# %%
BC_nz = get_non_zero(BC_df)

for vv in BC_nz:
    if vv not in ['poi_gage_id', 'K_coef', 'snarea_curve']:
        plot_param_3df(vv, BC_df, BD_df, CD_df, 'byHRU-musk', 'byHRU-musk_obs', 'musk - musk_obs')

# %%
AB_df

# %%
AB_df['tmax_cbh_adj'].T.dims

# %%
