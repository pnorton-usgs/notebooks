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
import  pyPRMS.ValidParams_v2 as vparm
from pyPRMS.ControlFile import ControlFile
reload(vparm)

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/fortran/prms6/src/xml'
workdir = '/Users/pnorton/tmp'

# Load valid parameters from xml
vpdb = vparm.ValidParams_v2()
# vpdb = vparm.ValidParams_v2('{}/parameters.xml'.format(workdir))

# %%
# Read the control file and get a list of requested modules
ctl = ControlFile('/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/control.default')

modules_used = []

for xx in ctl.modules.keys():
    if xx == 'precip_module':
        if ctl.modules[xx] == 'climate_hru':
            modules_used.append('precipitation_hru')
    elif xx == 'temp_module':
        if ctl.modules[xx] == 'climate_hru':
            modules_used.append('temperature_hru')
    else:
        modules_used.append(ctl.modules[xx])
print(modules_used)     

# %%
# Get the list of parameters that are required for the selected modules
required_params = vpdb.get_params_for_modules(modules=modules_used)
print('Number of required parameters: {}'.format(len(required_params)))

# %%
# Get the list of parameters that can be removed from the parameter files
remove_list = set(vpdb.parameters.keys()).difference(required_params)
print(remove_list)

# %%
# Build dictionary of parameters by module
params_by_module = {}

for xx in vpdb.parameters.values():
    for mm in xx.modules:
        if mm not in params_by_module:
            params_by_module[mm] = []
        params_by_module[mm].append(xx.name)

# %%
# Build list of parameters required for all modules
req_params = []
for xx, yy in params_by_module.iteritems():
    for pp in yy:
        req_params.append(pp)
        
print(len(req_params))

# Unique parameter names
print(len(set(req_params)))

# %%
# for xx in params_by_module.keys():
#     print('{}: {}'.format(xx, params_by_module[xx]))

# %%
import pyPRMS.NhmParamDb as nhm
import pandas as pd

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb'

pdb = nhm.NhmParamDb(workdir)

# %%
pdb.parameters.exists('print_type')

# %%
# curr_modules = ['basin', 'potet_jh', 'precipitation_hru', 'temperature_hru', 'soilzone', 
#                 'ddsolrad', 'srunoff_smidx', 'muskingum', 'transp_tindex', 'intcp', 'snowcomp',
#                 'gwflow']

for xx in params_by_module.keys():
    if xx in modules_used:
        print(xx)
        
        for yy in params_by_module[xx]:
            if not pdb.parameters.exists(yy):
                if yy in ['basin_solsta', 'hru_solsta', 'rad_conv']:
                    try:
                        if pdb.dimensions.get('nsol') > 0:
                            print('\tMissing: {}'.format(yy))
                    except ValueError:
                            pass
                elif yy == 'irr_type':
                    try:
                        if pdb.dimensions.get('nwateruse') == 1:
                            print('\tMissing: {}'.format(yy))
                    except ValueError:
                        pass
                elif yy == 'gvr_hru_id':
                    try:
                        if ctl.get_var('mapOutON_OFF') and ctl.get_values('mapOutON_OFF') == 1:
                            print('\tMissing: {}'.format(yy))
                    except ValueError:
                        pass
                else:
                    print('\tMissing: {}'.format(yy))

# %% [markdown]
# ### Read modules from control file

# %%
ctl = ControlFile('/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/control.default')

# %%
ctl.modules

# %%
modules_used = []

for xx in ctl.modules.keys():
    if xx == 'precip_module':
        if ctl.modules[xx] == 'climate_hru':
            modules_used.append('precipitation_hru')
    elif xx == 'temp_module':
        if ctl.modules[xx] == 'climate_hru':
            modules_used.append('temperature_hru')
    else:
        modules_used.append(ctl.modules[xx])
print(modules_used)     

# %%
remove_list = set(params_by_module.keys()).difference(set(modules_used))
print(remove_list)

# %%
# Create list of modules not used by current model
for cmod in remove_list:
    del params_by_module[cmod]
    
# Build list of parameters required for all modules
req_params = []
for xx, yy in params_by_module.iteritems():
    for pp in yy:
        req_params.append(pp)

# Unique parameter names
print(len(set(req_params)))
print(req_params)

# %%
print(pdb.dimensions)

# %%
pdb.dimensions.exists('nsol')

# %%
pdb.dimensions.exists('nsol')

# %%
pdb.dimensions.get('nsol') 

# %%
ctl

# %%
ctl.get('tmax_day').values

# %%
ctl.get('mapOutON_OFF').values

# %%
