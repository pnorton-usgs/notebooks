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
from future.utils import iteritems

import numpy as np
from pyPRMS.ControlFile import ControlFile
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ValidParams import ValidParams

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/conus'
filename = '{}/myparam.param'.format(workdir)

outfile = '{}/myparam.nc'.format(workdir)

# %%
pfile = ParameterFile(filename)

# %%

# %%
# Read the control file and get a list of requested modules
ctl = ControlFile('{}/control.netcdf'.format(workdir))

pfile.reduce_by_modules(control=ctl)

# %%
# Write the parameters to netcdf format
pfile.write_netcdf(filename=outfile)
