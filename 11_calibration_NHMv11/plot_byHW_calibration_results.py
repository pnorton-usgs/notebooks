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
headwater = '0259'
hw_suffix = ''
workdir = f'/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHW_sample/HW{headwater}{hw_suffix}/RESULTS'
ofs_file = f'{workdir}/objfun_{headwater}'

df = pd.read_csv(ofs_file, sep='\s+', skipinitialspace=True)        

x_vars = df.columns.tolist()[3:]

ncols = 3
numrows = int(round(len(x_vars) / float(ncols) + 0.5))

cstep = 4
# of_var = 'of_som'

# Layout info at: https://matplotlib.org/stable/tutorials/intermediate/constrainedlayout_guide.html
fig, axes = plt.subplots(nrows=numrows, ncols=ncols, figsize=(10, 10), constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0.1, wspace=0.2)
ax = axes.flatten()

for ii,of in enumerate(x_vars):
    ax[ii].set_title(f'of_prms vs {of}')
    step_df = df[df.step == cstep]
    step_df.plot(ax=ax[ii], kind='scatter', x=of, y='of_prms', color='red', alpha=0.2)
    
    df_final = step_df.iloc[[-1]]
    df_final.plot(ax=ax[ii], kind='scatter', x=of, y='of_prms', color='black')
    
#     precal_ns_ref_df.plot(ax=ax[0], x='OF', y=precal_ns_ref_df.columns[1], ylim=(0.0, 1.0), color=calib_color,
#                           label='PRECAL-ref')


# ax = plt.gca()
# step_df = df[df.step == cstep]
# df_final = step_df.iloc[[-1]]
# step_df.plot(kind='scatter', x=of_var, y='of_prms', ax=ax, color='red', alpha=0.2)
# df_final.plot(kind='scatter', x=of_var, y='of_prms', ax=ax, color='black')




# step_two = df[df.step == 2]
# step_two.plot(kind='scatter', x=of_var, y='of_prms', ax=ax, color='green', alpha=0.2)

# step_three = df[df.step == 3]
# step_three.plot(kind='scatter', x=of_var, y='of_prms', ax=ax, color='blue', alpha=0.2)

# step_four = df[df.step == 4]
# step_four.plot(kind='scatter', x=of_var, y='of_prms', ax=ax, color='yellow', alpha=0.2)

# df_final = step_one.iloc[[-1]]
# df_final.plot(kind='scatter', x='ofRUN', y='prmsOF', ax=ax, color='black')
# df_final.plot(kind='scatter', x='ofAET', y='prmsOF', ax=ax, color='black')
# df_final.plot(kind='scatter', x='ofSCA', y='prmsOF', ax=ax, color='black')
# df_final.plot(kind='scatter', x='ofRCH', y='prmsOF', ax=ax, color='black')
# df_final.plot(kind='scatter', x='ofSOM', y='prmsOF', ax=ax, color='black')


# %%
len(df.columns.tolist()[2:])

# %%
colors = ['red', 'green', 'blue', 'yellow']

ncols = 3
numrows = int(round(len(x_vars) / float(ncols) + 0.5))

rnd = 3
# of_var = 'of_som'
df = df[df.loc[:, 'round'] == rnd]

# Layout info at: https://matplotlib.org/stable/tutorials/intermediate/constrainedlayout_guide.html
fig, axes = plt.subplots(nrows=numrows, ncols=ncols, figsize=(15, 15), constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0.1, wspace=0.2)
ax = axes.flatten()

for ii,of in enumerate(x_vars):
    ax[ii].set_title(f'of_prms vs {of}')
    
    for xx in range(1, 5):
        p_df = df[df.step == xx]
        p_df.plot(ax=ax[ii], kind='scatter', x=of, y='of_prms', color=colors[xx-1], alpha=0.2)
    
    df_final = p_df.iloc[[-1]]
    df_final.plot(ax=ax[ii], kind='scatter', x=of, y='of_prms', color='black')

# %%
df[df.loc[:, 'round'] == 1]

# %%
df.head()

# %%
df.info()

# %%

# %%

# %%

# %%

# %%
x_vars

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
