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

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
from future.utils import iteritems

import os
import subprocess
import sys
import traceback

import pyPRMS.ParameterFile as prms
import pyPRMS.Control as prms_ctl
import  pyPRMS.ValidParams_v2 as vparm
# reload(prms)
# reload(prms_ctl)

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/prms'
cmd_prms = '/Users/pnorton/Projects/National_Hydrology_Model/code/PRMS/20180314_sntemp/bin/prms'

src_ctl_file = 'control.default'
wk_ctl_file = 'control.run'
src_param_file = 'myparams.param.orig'
wk_param_file = 'myparams.param'

output_dir = 'all_modules'

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
# Additional dimensions are required by certain modules
mod_dims = {'precip_dist2': {'nrain': 2}, 'potet_pan': {'nevap': 1}, 'temp_dist2': {'ntemp': 2}}

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
            ctl = prms_ctl.Control(src_ctl_file)

            ctl.replace_values(vv, mm)
            ctl.write_control_file(wk_ctl_file)

            # Modify the parameter file if necessary
            if mm in mod_dims:
                prm = prms.ParameterFile(src_param_file)
#                 myparams = pfile.parameters
                mydims = prm.dimensions
                
                for dd, dv in mod_dims[mm].iteritems():
                    if dd not in mydims.keys():
                        mydims.add(dd, dv)
#                     if prm.get_dim(dd) is None:
#                         prm.add_dimension(dd, dv)
                prm.write(wk_param_file)
            
            
            #   2) run 'prms -print'
            cmd_opts = ' -C{} -print'.format(wk_ctl_file)
            result = subprocess.call(cmd_prms + cmd_opts, shell=True)

            #   3) save *.par_name to a numbered file in the output directory
            os.rename('{}.par_name'.format(wk_ctl_file), '{:02d}.par_name'.format(num))
            num += 1
except NameError as err:
    exc_type, exc_value, exc_traceback = sys.exc_info()
    traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
    print('Blah!')
    print(err)
#     print(sys.exc_info()[0])
finally:
    os.chdir(orig_dir)


# %%
pdb = vparm.ValidParams_v2('{}/{}'.format(workdir, output_dir))

myparams = pdb.parameters
mydims = pdb.dimensions

# %%
print(len(myparams.keys()))

# %%
print(mydims)

# %%
print(myparams.get('poi_gage_id'))

# %%
print(myparams.get('poi_gage_id').minimum)
print(mydims.get('nwateruse').description)

# %% [markdown]
# ### Get global dimensions in xml format

# %%
import xml.dom.minidom as minidom
import xml.etree.ElementTree as xmlET

xmlstr = minidom.parseString(xmlET.tostring(pdb.xml_global_dimensions)).toprettyxml(indent='    ')
print(xmlstr)

pdb.write_dimensions_xml('/Users/pnorton/tmp')
# with open('crap.xml', 'w') as f:
#     f.write(xmlstr.encode('utf-8'))

# %% [markdown]
# ### Get parameters in xml format

# %%
xmlstr = minidom.parseString(xmlET.tostring(pdb.xml_global_parameters)).toprettyxml(indent='    ')
print(xmlstr)

pdb.write_parameters_xml('/Users/pnorton/tmp')

# %%

# %%
# Loop through the parameter objects in the ParameterSet and display name, datatype
# Create dictionary of parameter names to datatypes
param_to_datatype = {}

for vv in myparams.values():
    param_to_datatype[vv.name] = vv.datatype
    print('{}: datatype = {}'.format(vv.name, vv.datatype))


# %%
param_to_datatype

# %%
for ii, jj in iteritems(mydims):
    print ii, jj.size

# %%
