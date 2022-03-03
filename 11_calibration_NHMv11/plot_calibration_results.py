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
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHRU_sample/HRU3505_run2/RESULTS'
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

# %%

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
