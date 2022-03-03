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
from future.utils import iteritems

import pyPRMS.ParameterFile as pf
reload(pf)

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/prms5/test4'
filename = '{}/myparam.param'.format(workdir)

# filename = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/regions/subset_testing/crap.param'

pfile = pf.ParameterFile(filename)

# Set pointers to the parameters and dimensions within the parameterSet
# This simplifies code for access to the members of parameters and dimensions
myparams = pfile.parameters
mydims = pfile.dimensions

# %%
curr_param = 'dday_slope'

# Here's one way to get parameter information (three approaches)
# print(pfile.parameterset.parameters.get('gwflow_coef'))
print(myparams.get(curr_param))

# %%

# %%
# Convert relative temperatures (or temperature differences) from Fahrenheit to Celsius
temp_diffs = ['tmax_cbh_adj', 'tmin_cbh_adj', 'dday_slope', 'radadj_slope', 'tmax_allrain_offset']
ftoc_rel = 5.0 / 9.0

for xx in temp_diffs:
    aa = myparams.get(xx).tolist()
    
    myparams.get(xx).data = [bb * ftoc_rel for bb in aa]
    

# %%
# Convert absolute temperatures from Fahrenheit to Celsius
temp_abs = ['tmax_index', 'transp_tmax', 'tmax_allsnow']

def f_to_c(temp_f):
    return (temp_f - 32.) / 1.8

for xx in temp_abs:
    aa = myparams.get(xx).tolist()
    
    myparams.get(xx).data = [f_to_c(bb) for bb in aa]

# %%
myparams.get('dday_intcp').data

# %%
