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
basedir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/sagehen'
workdir_prms5 = f'{basedir}/prms5/output/netcdf'
workdir_prms6 = f'{basedir}/prms6/output'

# %%
prms5_file = f'{workdir_prms5}/*.nc'
prms6_file = f'{workdir_prms6}/summary_daily.nc'

# %%
# ds = xr.open_mfdataset(f'{workdir}/*.nc', chunks={'hruid': 1040}, combine='by_coords', decode_cf=True,
#                        attrs_file=f'{workdir}/gm_climate_1980.nc')

# p5_df = xr.open_dataset(prms5_file)
p5_df = xr.open_mfdataset(prms5_file, combine='by_coords')
p5_df

# %%

# %%
p6_df = xr.open_dataset(prms6_file)
p6_df

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(17,11))
ax = axes.flatten()

p6_df['seg_outflow'].plot(ax=ax[0])
p5_df['seg_outflow'].plot(ax=ax[1])
(p6_df['seg_outflow'] - p5_df['seg_outflow']).plot(ax=ax[2], cmap='bwr', vmin=-0.5, vmax=0.5)

# ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim([1, 15])
ax[1].set_xlim([1, 15])
ax[2].set_xlim([1, 15])

ax[0].set_title('PRMS6', fontsize=12)
ax[1].set_title('PRMS5', fontsize=12)
ax[2].set_title('PRMS6 - PRMS5', fontsize=12)

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(17,11))
ax = axes.flatten()

p6_df['hru_outflow'].plot(ax=ax[0], cmap='terrain_r')
p5_df['hru_outflow'].plot(ax=ax[1], cmap='terrain_r')
(p6_df['hru_outflow'] - p5_df['hru_outflow']).plot(ax=ax[2], cmap='bwr', vmin=-0.01, vmax=0.01)

# ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim([1, 15])
ax[1].set_xlim([1, 15])
ax[2].set_xlim([1, 15])

ax[0].set_title('PRMS6', fontsize=12)
ax[1].set_title('PRMS5', fontsize=12)
ax[2].set_title('PRMS6 - PRMS5', fontsize=12)

# %%
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(17,11))
# ax = axes.flatten()

p5_df['seg_outflow'][2555:4015, 13].plot(ax=axes, color='#666666', linewidth=1, alpha=1)
p6_df['seg_outflow'][2555:4015, 13].plot(ax=axes, color='#cc3333', linewidth=0.5)

# axes.plot(p5_df['seg_outflow'][:, 4])
# axes.plot(p6_df['seg_outflow'][:, 4], linewidth=0.5, alpha=0.8)
# ax[0].plot(sim_mon_mn.index.to_pydatetime(), sim_mon_mn, linewidth=.25, color='#cc3333', label='MOCOM', alpha=0.5)

# ax[0].set_title('Monthly Streamflow', fontsize=12)
# plt.suptitle('basin: %s\nrunid: %s' % (basinid, runid), fontsize=15)

# ax[1].xaxis.set_major_locator(YearLocator(base=1, month=1, day=1))
# ax[1].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(17,11))
ax = axes.flatten()

p6_df['recharge'].plot(ax=ax[0], cmap='terrain_r')
p5_df['recharge'].plot(ax=ax[1], cmap='terrain_r')
(p6_df['recharge'] - p5_df['recharge']).plot(ax=ax[2], cmap='bwr', vmin=-0.001, vmax=0.001)

# ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim([1, 15])
ax[1].set_xlim([1, 15])
ax[2].set_xlim([1, 15])

ax[0].set_title('PRMS6', fontsize=12)
ax[1].set_title('PRMS5', fontsize=12)
ax[2].set_title('PRMS6 - PRMS5', fontsize=12)

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(17,11))
ax = axes.flatten()

p6_df['soil_moist'].plot(ax=ax[0], cmap='terrain_r')
p5_df['soil_moist'].plot(ax=ax[1], cmap='terrain_r')
(p6_df['soil_moist'] - p5_df['soil_moist']).plot(ax=ax[2], cmap='bwr', vmin=-0.01, vmax=0.01)

# ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim([1, 15])
ax[1].set_xlim([1, 15])
ax[2].set_xlim([1, 15])

ax[0].set_title('PRMS6', fontsize=12)
ax[1].set_title('PRMS5', fontsize=12)
ax[2].set_title('PRMS6 - PRMS5', fontsize=12)

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(17,11))
ax = axes.flatten()

p6_df['infil'].plot(ax=ax[0], cmap='terrain_r')
p5_df['infil'].plot(ax=ax[1], cmap='terrain_r')
(p6_df['infil'] - p5_df['infil']).plot(ax=ax[2], cmap='bwr', vmin=-0.01, vmax=0.01)

# ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim([1, 15])
ax[1].set_xlim([1, 15])
ax[2].set_xlim([1, 15])

ax[0].set_title('PRMS6', fontsize=12)
ax[1].set_title('PRMS5', fontsize=12)
ax[2].set_title('PRMS6 - PRMS5', fontsize=12)

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(17,11))
ax = axes.flatten()

p6_df['potet'].plot(ax=ax[0], cmap='terrain_r')
p5_df['potet'].plot(ax=ax[1], cmap='terrain_r')
(p6_df['potet'] - p5_df['potet']).plot(ax=ax[2], cmap='bwr', vmin=-0.001, vmax=0.001)

# ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim([1, 128])
ax[1].set_xlim([1, 128])
ax[2].set_xlim([1, 128])

ax[0].set_title('PRMS6', fontsize=12)
ax[1].set_title('PRMS5', fontsize=12)
ax[2].set_title('PRMS6 - PRMS5', fontsize=12)

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(17,11))
ax = axes.flatten()

p6_df['hru_storage'].plot(ax=ax[0], cmap='terrain_r')
p5_df['hru_storage'].plot(ax=ax[1], cmap='terrain_r')
(p6_df['hru_storage'] - p5_df['hru_storage']).plot(ax=ax[2], cmap='bwr', vmin=-0.1, vmax=0.1)

# ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim([1, 128])
ax[1].set_xlim([1, 128])
ax[2].set_xlim([1, 128])

ax[0].set_title('PRMS6', fontsize=12)
ax[1].set_title('PRMS5', fontsize=12)
ax[2].set_title('PRMS6 - PRMS5', fontsize=12)

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(17,11))
ax = axes.flatten()

p6_df['hru_streamflow_out'].plot(ax=ax[0], cmap='terrain_r')
p5_df['hru_streamflow_out'].plot(ax=ax[1], cmap='terrain_r')
(p6_df['hru_streamflow_out'] - p5_df['hru_streamflow_out']).plot(ax=ax[2], cmap='bwr', vmin=-0.01, vmax=0.01)

# ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim([1, 128])
ax[1].set_xlim([1, 128])
ax[2].set_xlim([1, 128])

ax[0].set_title('PRMS6', fontsize=12)
ax[1].set_title('PRMS5', fontsize=12)
ax[2].set_title('PRMS6 - PRMS5', fontsize=12)

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(17,11))
ax = axes.flatten()

p6_df['frac_swe'].plot(ax=ax[0], cmap='terrain_r')
p5_df['frac_swe'].plot(ax=ax[1], cmap='terrain_r')
(p6_df['frac_swe'] - p5_df['frac_swe']).plot(ax=ax[2], cmap='bwr', vmin=-0.01, vmax=0.01)

# ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim([1, 128])
ax[1].set_xlim([1, 128])
ax[2].set_xlim([1, 128])

ax[0].set_title('PRMS6', fontsize=12)
ax[1].set_title('PRMS5', fontsize=12)
ax[2].set_title('PRMS6 - PRMS5', fontsize=12)

# %%
p5_df['frac_swe'][:,87].plot()

# %%
p5_swe = p5_df['frac_swe'][:,:].values

print(f'v5 max frac_swe: {np.nanmax(p5_swe)}')
print(f'v5 num frac_swe > 1: {p5_swe[p5_swe > 1.0].size}')
# np.nanmax(bb[bb != np.inf])

# %%
p6_swe = p6_df['frac_swe'][:,:].values

print(f'v6 max frac_swe: {np.nanmax(p6_swe)}')
print(f'v6 num frac_swe > 1: {p6_swe[p6_swe > 1.0].size}')

# %%
p6_ai = p6_df['ai'].values
print(f'v6 max ai: {np.nanmax(p6_ai)}')

# %%
p6_pst = p6_df['pst'].values
print(f'v6 max pst: {np.nanmax(p6_pst)}')

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(17,11))
ax = axes.flatten()

p6_df['ai'].plot(ax=ax[0], cmap='terrain_r')
p6_df['pst'].plot(ax=ax[1], cmap='terrain_r')
(p6_df['ai'] - p6_df['pst']).plot(ax=ax[2], cmap='bwr', vmin=-0.01, vmax=0.01)

# ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim([1, 128])
ax[1].set_xlim([1, 128])
ax[2].set_xlim([1, 128])

ax[0].set_title('ai', fontsize=12)
ax[1].set_title('pst', fontsize=12)
ax[2].set_title('ai - pst', fontsize=12)

# %%
xlim = [1,128]
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(17,11))
ax = axes.flatten()

p5_df['ai'].plot(ax=ax[0], cmap='terrain_r')
p5_df['pst'].plot(ax=ax[1], cmap='terrain_r')
p5_df['pkwater_equiv'].plot(ax=ax[2], cmap='terrain_r')

# ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim(xlim)
ax[1].set_xlim(xlim)
ax[2].set_xlim(xlim)

ax[0].set_title('ai', fontsize=12)
ax[1].set_title('pst', fontsize=12)
ax[2].set_title('pkwater_equiv', fontsize=12)

# %%
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20, 6))
ax = axes.flatten()

hru = 70
p5_df['ai'][:, hru].plot(ax=ax[0])
p5_df['pst'][:, hru].plot(ax=ax[1])
p5_df['pkwater_equiv'][:, hru].plot(ax=ax[2])


# %%
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(20, 6))
# ax = axes.flatten()

hru = 95

p5_df['pkwater_equiv'][:, hru].plot(ax=axes, label='pkwater_equiv', color='gray', linewidth=.9)
p5_df['pst'][:, hru].plot(ax=axes, label='pst', color='goldenrod', linewidth=6.5, alpha=.7)
p5_df['ai'][:, hru].plot(ax=axes, label='ai', color='dodgerblue', linewidth=1.5, alpha=1.0)
plt.legend(loc="upper left")

# %%

# %%
