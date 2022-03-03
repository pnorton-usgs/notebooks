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

# %%

# %%
import pandas as pd

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/slurms_redo'

notesdir = f'{workdir}/notes'


# %%
# HRUs missing the *PARAMS* file did not finish calibration to completion
df_miss_params = pd.read_csv(f'{notesdir}/missing_PARAMS_file.txt', sep=',', header=0)
df_miss_params.head()

miss_params = df_miss_params.loc[:,'file'].tolist()
hru_miss_params = df_miss_params.loc[:, 'hru'].tolist()

# %%
print(len(set(miss_params)))
print(len(set(hru_miss_params)))
print(len(hru_miss_params))

# %%
# HRUs which generate a bad index error most likely are missing usable baseline data
df_bad_index = pd.read_csv(f'{notesdir}/bad_index.txt', header=0)
df_bad_index.head()

bad_index = df_bad_index.loc[:, 'file'].tolist()

# %%
# Missing PARAMS file and had a 'bad index' error

set(miss_params) & set(bad_index)

# %%
# Calibrations that tripped an IEEE_INVALID_FLAG while running
df_ieee = pd.read_csv(f'{notesdir}/IEEE_INVALID_FLAG.txt', header=0)
df_ieee.head()

ieee_flag = df_ieee.loc[:, 'file'].tolist()

# %%
# Incomplete calibrations that tripped the IEEE_INVALID_FLAG
set(miss_params) & set(ieee_flag)

# %%

# %%
df_completed = pd.read_csv(f'{notesdir}/finished_hrus.txt', sep=',', header=0)
df_completed.head()

completed = df_completed.loc[:, 'file'].tolist()

# %%
# Jobs that finished calibration for at least some of the HRUs but raised the IEEE_INVALID_FLAG
set(completed) & set(ieee_flag)

# %%
# Jobs that finished calibration but had one or more incomplete HRU
set(completed) & set(miss_params)

# %%
set(completed) & set(miss_params) & set(ieee_flag)

# %%
df_bad_RCH = pd.read_csv(f'{notesdir}/hru_bad_RCH_baseline.txt', header=0)
df_bad_RCH.head()

bad_RCH = df_bad_RCH.loc[:,'hru'].tolist()

# %%
df_hru_bad_index = pd.read_csv(f'{notesdir}/hru_bad_index.txt', sep=',', header=0)

hru_bad_index = df_hru_bad_index.loc[:, 'hru'].tolist()

# %%
# HRUs that generated the bad index error but did not have bad RCH values
set(hru_bad_index) - set(bad_RCH)

# %%
# HRUs have bad RCH values are not in the list of HRUs which generated the bad index error
set(bad_RCH) - set(hru_bad_index)

# %%
len(set(hru_bad_index) - set(hru_miss_params))

# %%
len(hru_miss_params)

# %%
df_bad_proc_wgt = pd.read_csv(f'{notesdir}/bad_process_weights.txt', sep=',', header=0)
df_bad_proc_wgt.head()

hru_bad_proc_wgt = df_bad_proc_wgt.loc[:, 'hru'].tolist()

# %%
len(set(hru_bad_proc_wgt))

# %%
set(hru_bad_proc_wgt) - set(hru_miss_params)

# %%
set(hru_miss_params) - set(hru_bad_proc_wgt)

# %%
