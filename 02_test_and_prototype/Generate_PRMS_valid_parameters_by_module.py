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
#     display_name: Python [default]
#     language: python
#     name: python2
# ---

# %%

# %%
workdir = '/media/scratch/PRMS/datasets/prms_param_work'
cmd_prms = '/media/scratch/PRMS/bin/prms'

# %% [markdown]
# ### Build dictionary of modules by control file module variable

# %%

# Maps process name to the control file module variables
mod_map = {'Temperature Dist': 'temp_module', 
           'Precip Dist': 'precip_module',
           'Temp & Precip Dist': ['precip_module', 'temp_module'],
           'Solar Rad Dist': 'solrad_module', 
           'Transpiration Dist': 'transp_module',
           'Potential ET': 'et_module', 
           'Surface Runoff': 'srunoff_module',
           'Streamflow Routing': 
           'strmflow_module'}

valid_procs = ['Basin Definition', 'Cascading Flow', 'Time Series Data', 'Potet Solar Rad',
               'Temperature Dist', 'Precip Dist', 'Temp & Precip Dist', 'Solar Rad Dist',
               'Transpiration Dist', 'Potential ET', 'Interception', 'Snow Dynamics',
               'Surface Runoff', 'Soil Zone', 'Groundwater', 'Streamflow Routing',
               'Output Summary', 'Preprocessing']

# Maps the control file variables to the available modules
modvar_modules = {}

# TODO: Run prms -print

# Read model.out
# Build dictionary of modules by control file module variable
fhdl = open('{}/model.out'.format(workdir), 'r')

# Look for line before the modules are listed
for rec in fhdl:
    if rec[0] == '-':
        break

for rec in fhdl:
    if rec[0] == '-':
        break
    flds = rec.strip().split(':')
    if flds[0] in valid_procs:
        mods = flds[1].strip().rstrip(',').split(', ')
        
        if flds[0] in mod_map.keys():
            if isinstance(mod_map[flds[0]], list):
                for ff in mod_map[flds[0]]:
                    if ff not in modvar_modules:
                        modvar_modules[mod_map[ff]] = []
                                   
                    for mm in mods:
                        modvar_modules[ff].append(mm)
            else:
                modvar_modules[mod_map[flds[0]]] = []

                for mm in mods:
                    modvar_modules[mod_map[flds[0]]].append(mm)
    else:
        mods = flds[0].strip().split(',')

        if prior in modvar_modules:
            for mm in mods:
                modvar_modules[mod_map[prior]].append(mm)

    prior = flds[0]
#     print(rec.strip().split(':'))

print(modvar_modules)

fhdl.close()

# %% [markdown]
# ### Run 'prms -print' for each module

# %%

import os
import subprocess
import sys
import traceback

import prms_lib as prms

# %%
# Additional dimensions are required by certain modules
mod_dims = {'precip_dist2': {'nrain': 2}, 'potet_pan': {'nevap': 1}, 'temp_dist2': {'ntemp': 2}}

# %%
src_ctl_file = 'control.default'
wk_ctl_file = 'control.run'
src_param_file = 'myparams.param.orig'
wk_param_file = 'myparams.param'

output_dir = 'all_mods'
num = 1

# change working directory to prms model area
orig_dir = os.getcwd()

try:
    os.chdir(workdir)
    print(os.getcwd())

    # loop through each module variable. For each module:
    for vv, mods in modvar_modules.iteritems():
        print('Module variable: {}'.format(vv))
        
        for mm in mods:
            print('\t{}'.format(mm))
            #   1) modify the control file
            ctl = prms.control(src_ctl_file)

            ctl.replace_values(vv, mm)
            ctl.write_control_file(wk_ctl_file)

            # Modify the parameter file if necessary
            if mm in mod_dims:
                prm = prms.parameters(src_param_file)
                
                for dd, dv in mod_dims[mm].iteritems():
                    if prm.get_dim(dd) is None:
                        prm.add_dimension(dd, dv)
                prm.write_param_file(wk_param_file)
            
            
            #   2) run 'prms -print'
            cmd_opts = ' -C{} -print'.format(wk_ctl_file)
            result = subprocess.call(cmd_prms + cmd_opts, shell=True)

            #   3) save *.par_name to a numbered file in the output directory
            os.rename('{}.par_name'.format(wk_ctl_file), '{0:2d}.par_name'.format(num))
            num += 1
except NameError as err:
    exc_type, exc_value, exc_traceback = sys.exc_info()
    traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
    print('Blah!')
    print(err)
#     print(sys.exc_info()[0])
finally:
    os.chdir(orig_dir)


# %% [markdown]
# ### Show every parameter and the modules they are used in

# %%
pdb = prms.param_db('{}/all_mods'.format(workdir))

# Output every parameter and the modules that use it
for kk, vv in pdb.paramdb.iteritems():
    print('{}:'.format(kk)),
    
    for mm in vv['Module']:
        print('{}'.format(mm)),
    print('')

# %%
print(pdb.module_params({'climate_hru': 'temp_module'}))

# %%
