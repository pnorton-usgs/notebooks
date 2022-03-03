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
import numpy as np
from collections import OrderedDict
from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParamDbRegion import ParamDbRegion
from pyPRMS.ParameterSet import ParameterSet

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/lauren_repo/Daymet'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/paramdb_daymet_CONUS'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_lauren_daymet_byHRU'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/nhmparamdb_lhay_daymet_byHRU'

pdb = ParamDbRegion(paramdb_dir=workdir, verbose=True, verify=True)

# %%
param_name = 'tmax_cbh_adj'

# %%
pdb.parameters[param_name].data.shape

# %%

for xx in pdb.parameters[param_name].data[2118,:]:
    print(xx)

# %%

# %%

# %% [markdown]
# ## By-region parameter database

# %%
# Read from a by-region paramdb file
region_workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/lauren_repo/Maurer'
# region_workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/lauren_repo/Daymet'
# region_workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb'

# %%
pdb_r = ParamDbRegion(region_workdir, verbose=True, verify=True)

# %%
pdb_r.parameters[param_name].data[2118]

# %%
pdb_r.parameters[param_name].data[2118,:]

# %%
pdb_r.parameters[param_name].data[1962,:]

# %%
pdb_r.parameters[param_name].data[0,:]

# %%
pdb_r.dimensions['nhru'].size

# %%
pdb_r.parameters.check()

# %%
pdb_r.parameters[param_name].stats()

# %%
pdb_r.parameters[param_name].maximum

# %%
