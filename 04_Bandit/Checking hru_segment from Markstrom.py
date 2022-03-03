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
#     display_name: Python [conda env:bandit_py3]
#     language: python
#     name: conda-env-bandit_py3-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import pandas as pd
import numpy as np
import pydot
import networkx as nx

from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ParamDb import ParamDb

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190208_DelawareRiver'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190913_Delware_streamtemp'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190916_Delaware_CONUS'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/streamtemp_paramdb/source'
filename = '{}/myparam.param'.format(workdir)

# Load parameter file
pfile = ParameterFile(filename)

# %%
workdir_CONUS = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/streamtemp_paramdb'
# workdir_CONUS = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_CONUS'

pdb = ParamDb(paramdb_dir=workdir_CONUS, verbose=True, verify=True)

# %%
master = pdb.parameters['hru_segment_nhm'].data

# %%
comp = pfile.parameters['hru_segment_nhm'].data

# %%
for ii, vv in enumerate(master):
#     if vv != comp[ii]:
    if vv == 0 or comp[ii] == 0:
        print('index: {}; master = {}, src = {}'.format(ii+1, vv, comp[ii]))

# %%
for ii, vv in enumerate(comp):
    if vv == 0:
        print('index: {}; src = {}'.format(ii+1, vv))

# %%
nr_master = 0
for ii, vv in enumerate(master):
    if vv == 0:
        nr_master += 1
        
print('Non-routed HRUs: {}'.format(nr_master))

# %%
nr_comp = 0
for ii, vv in enumerate(comp):
    if vv == 0:
        nr_comp += 1
        
print('Non-routed HRUs: {}'.format(nr_comp))

# %%
comp_nhm_id = pfile.parameters['nhm_id'].data


# %%
aa = set(comp_nhm_id)

# %%
len(aa)

# %%
nhm_id = pfile.parameters['nhm_id'].data

# %%
min(nhm_id)

# %%
max(nhm_id)

# %%
pfile.degenerate_parameters()

# %%
aa = 'seginc_gwflow'
aa[0:3] == 'seg'

# %%
