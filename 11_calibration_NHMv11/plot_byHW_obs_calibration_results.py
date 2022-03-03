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
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

# %%

# %%

# %%

# %%

# %%

# %%

# %%
col_names = ['HW', 'run_num', 'poi_id', 'of_prms', 'NS', 'NSlog', 'weight', 'weightNUM', 'sim_da_ratio']

# NOTE: Must use Int64 dtype for nullable-integer fields
col_types = [str, int, str, float, float, float, float, float, float]

cols = dict(zip(col_names, col_types))


# %%

# %%
headwater = '0833'
hw_suffix = ''
workdir = f'/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHW_obs_sample/hw_{headwater}{hw_suffix}/RESULTS'
ofs_file = f'{workdir}/objfun_{headwater}'

df = pd.read_csv(ofs_file, sep='\s+', skipinitialspace=True, dtype=cols)        

x_vars = df.columns.tolist()[2:]

poi_id = '02157470'

# %%

# %%
pois = df.loc[:, 'poi_id'].unique().tolist()

xmax = round(df.NS.max() + 0.4)
xmin = round(df.NS.min() - 0.4)

# Layout info at: https://matplotlib.org/stable/tutorials/intermediate/constrainedlayout_guide.html
ncols = 4
numrows = int(round(len(pois) / float(ncols) + 0.4))

fig, axes = plt.subplots(nrows=numrows, ncols=ncols, figsize=(20, 10), constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0.1, wspace=0.2)
ax = axes.flatten()

for ii, poi in enumerate(pois):
    ax[ii].set_title(f'{poi}')
    tmp_df = df[df.poi_id == poi]
    tmp_df.plot(ax=ax[ii], kind='scatter', x='NS', y='of_prms', color='red', alpha=0.2)
    
    df_final = tmp_df.iloc[[-1]]
    df_final.plot(ax=ax[ii], kind='scatter', x='NS', y='of_prms', color='black')
    ax[ii].set_xlim([xmin, xmax])

# %%
xmax = round(df.NSlog.max() + 0.4)
xmin = round(df.NSlog.min() - 0.4)

# Layout info at: https://matplotlib.org/stable/tutorials/intermediate/constrainedlayout_guide.html
ncols = 4
numrows = int(round(len(pois) / float(ncols) + 0.4))

fig, axes = plt.subplots(nrows=numrows, ncols=ncols, figsize=(20, 10), constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0.1, wspace=0.2)
ax = axes.flatten()

for ii, poi in enumerate(pois):
    ax[ii].set_title(f'{poi}')
    tmp_df = df[df.poi_id == poi]
    tmp_df.plot(ax=ax[ii], kind='scatter', x='NSlog', y='of_prms', color='red', alpha=0.2)
    
    df_final = tmp_df.iloc[[-1]]
    df_final.plot(ax=ax[ii], kind='scatter', x='NSlog', y='of_prms', color='black')
    ax[ii].set_xlim([xmin, xmax])

# %% [markdown]
# ### Plot OFS from the original byHRU calibration

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHRU_sample/HRU3505_run1/RESULTS'
ofs_file = f'{workdir}/OFS_HRU3505'
df = pd.read_csv(ofs_file, sep='\s+', skipinitialspace=True, header=0)     

# df.plot(kind='scatter',x='num_children',y='num_pets',color='red')
ax = plt.gca()

df.plot(kind='scatter', x='ofRUN', y='prmsOF', ax=ax, color='red', alpha=0.2)
df.plot(kind='scatter', x='ofAET', y='prmsOF', ax=ax, color='green', alpha=0.2)
df.plot(kind='scatter', x='ofSCA', y='prmsOF', ax=ax, color='blue', alpha=0.2)
df.plot(kind='scatter', x='ofRCH', y='prmsOF', ax=ax, color='yellow', alpha=0.2)
df.plot(kind='scatter', x='ofSOM', y='prmsOF', ax=ax, color='purple', alpha=0.2)

df_final = df.iloc[[-1]]
df_final.plot(kind='scatter', x='ofRUN', y='prmsOF', ax=ax, color='black')
df_final.plot(kind='scatter', x='ofAET', y='prmsOF', ax=ax, color='black')
df_final.plot(kind='scatter', x='ofSCA', y='prmsOF', ax=ax, color='black')
df_final.plot(kind='scatter', x='ofRCH', y='prmsOF', ax=ax, color='black')
df_final.plot(kind='scatter', x='ofSOM', y='prmsOF', ax=ax, color='black')

# %%

# %% [markdown]
# ### Plot params

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHRU_sample/HRU3505_run2/RESULTS'
ofs_file = f'{workdir}/PARAMS_HRU3505'

df = pd.read_csv(ofs_file, sep='\s+', skipinitialspace=True, header=0)        

ax = plt.gca()

df.plot(kind='scatter', x='carea_max', y='RUN', ax=ax, color='red', alpha=0.2)
df.plot(kind='scatter', x='fastcoef_lin', y='RUN', ax=ax, color='green', alpha=0.2)
df.plot(kind='scatter', x='freeh2o_cap', y='RUN', ax=ax, color='blue', alpha=0.2)
df.plot(kind='scatter', x='gwflow_coef', y='RUN', ax=ax, color='yellow', alpha=0.2)
df.plot(kind='scatter', x='jh_coef', y='RUN', ax=ax, color='purple', alpha=0.2)

df_final = df.iloc[[-1]]
df_final.plot(kind='scatter', x='carea_max', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x='fastcoef_lin', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x='freeh2o_cap', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x='gwflow_coef', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x='jh_coef', y='RUN', ax=ax, color='black')

# %% [markdown]
# ### Plot params from original calibration

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHRU_sample/HRU3505_run1/RESULTS'
ofs_file = f'{workdir}/PARAMS_HRU3505'

df = pd.read_csv(ofs_file, sep='\s+', skipinitialspace=True, header=0)        

ax = plt.gca()

df.plot(kind='scatter', x='carea_max', y='RUN', ax=ax, color='red', alpha=0.2)
df.plot(kind='scatter', x='fastcoef_lin', y='RUN', ax=ax, color='green', alpha=0.2)
df.plot(kind='scatter', x='freeh2o_cap', y='RUN', ax=ax, color='blue', alpha=0.2)
df.plot(kind='scatter', x='gwflow_coef', y='RUN', ax=ax, color='yellow', alpha=0.2)
df.plot(kind='scatter', x='jh_coef', y='RUN', ax=ax, color='purple', alpha=0.2)

df_final = df.iloc[[-1]]
df_final.plot(kind='scatter', x='carea_max', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x='fastcoef_lin', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x='freeh2o_cap', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x='gwflow_coef', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x='jh_coef', y='RUN', ax=ax, color='black')

# %%
ax = plt.gca()

df.plot(kind='scatter', x='carea_max', y='RUN', ax=ax, color='red', alpha=0.2)
df_final = df.iloc[[-1]]
df_final.plot(kind='scatter', x='carea_max', y='RUN', ax=ax, color='black')

# %%
df.columns

# %%

# %%

# %%
var = 'tmin_cbh_adj'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHRU_sample/HRU3505_run2/RESULTS'
ofs_file = f'{workdir}/PARAMS_HRU3505'

df = pd.read_csv(ofs_file, sep='\s+', skipinitialspace=True, header=0)        

ax = plt.gca()

df.plot(kind='scatter', x=f'{var}', y='RUN', ax=ax, color='red', alpha=0.2)
df.plot(kind='scatter', x=f'{var}.1', y='RUN', ax=ax, color='green', alpha=0.2)
df.plot(kind='scatter', x=f'{var}.2', y='RUN', ax=ax, color='blue', alpha=0.2)
df.plot(kind='scatter', x=f'{var}.3', y='RUN', ax=ax, color='purple', alpha=0.2)

df_final = df.iloc[[-1]]
df_final.plot(kind='scatter', x=f'{var}', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x=f'{var}.1', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x=f'{var}.2', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x=f'{var}.3', y='RUN', ax=ax, color='black')

# %%

# %%
var = 'tmin_cbh_adj'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHRU_sample/HRU3505_run1/RESULTS'
ofs_file = f'{workdir}/PARAMS_HRU3505'

df = pd.read_csv(ofs_file, sep='\s+', skipinitialspace=True, header=0)        

ax = plt.gca()

df.plot(kind='scatter', x=f'{var}', y='RUN', ax=ax, color='red', alpha=0.2)
df.plot(kind='scatter', x=f'{var}.1', y='RUN', ax=ax, color='green', alpha=0.2)
df.plot(kind='scatter', x=f'{var}.2', y='RUN', ax=ax, color='blue', alpha=0.2)
df.plot(kind='scatter', x=f'{var}.3', y='RUN', ax=ax, color='purple', alpha=0.2)

df_final = df.iloc[[-1]]
df_final.plot(kind='scatter', x=f'{var}', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x=f'{var}.1', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x=f'{var}.2', y='RUN', ax=ax, color='black')
df_final.plot(kind='scatter', x=f'{var}.3', y='RUN', ax=ax, color='black')

# %%
cnt = 0

for xx in range(1981, 2010):
    if xx%2 == 1:
        print(xx)
        cnt += 1
print(cnt)

# %%
