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
#     display_name: Python [conda env:bandit]
#     language: python
#     name: conda-env-bandit-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import glob
import os

# %%
workdir_tb = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/tbparamdb'
workdir_v1 = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/v1_paramdb'

# %%
tb_files = []

file_it = glob.glob('{}/*.csv'.format(workdir_tb))
for kk in file_it:
    tb_files.append(os.path.basename(os.path.splitext(kk)[0]))

# %%
v1_files = []

file_it = glob.glob('{}/*.csv'.format(workdir_v1))
for kk in file_it:
    v1_files.append(os.path.basename(os.path.splitext(kk)[0]))

# %%
# Files in v1 or TB but not in both
set(tb_files) ^ set(v1_files)

# %%
