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
#     display_name: Python 2
#     language: python
#     name: python2
# ---

# %%
import re
import glob

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/code/PRMS/20151229_rev7844/build/prms'
filename = 'call_modules.f90'

defaults = {'data_file': 'prms.data', 'model_mode': 'PRMS', 'model_output_file': 'prms.out', 
            'param_file': 'prms.params', 'start_time': [2000,10,1,0,0,0], 'end_time': [2001,9,30,0,0,0],
            'et_module': 'potet_jh', 'precip_module': 'climate_hru', 'soilzone_module': 'soilzone',
            'solrad_module': 'ddsolrad', 'srunoff_module': 'srunoff_smidx',
            'strmflow_module': 'muskingum', 'temp_module': 'climate_hru', 'transp_module': 'transp_tindex',
            'humidity_day': 'humidity.day', 'potet_day': 'potet.day', 'precip_day': 'precip.day',
            'swrad_day': 'swrad.day', 'tmax_day': 'tmax.day', 'tmin_day': 'tmin.day', 
            'transp_day': 'transp.day', 'windspeed_day': 'windspeed.day', 'stat_var_file': 'statvar.out',
            'var_init_file': 'prms_ic.in', 'var_save_file': 'prms_ic.out', 
            'ani_output_file': 'animation.out', 'executable_desc': 'MOWS executable',
            'executable_model': 'prmsIV', 'initial_deltat': 24.0, 'dispGraphsBuffSize': 50}

# %%
filelist = [el for el in glob.glob('%s/*.f*' % workdir)]

# %%
# Get integer control fields
fhdl = open('%s/%s' % (workdir, filename), 'r')

ctl_dict = {}
module_dict = {}

for ff in filelist:
    fhdl = open(ff, 'r')
    
    for line in fhdl:
        ctl_type = 0
        flds = None
        if "control_integer(" in line:
            flds = re.split(",|'", line)[2]
            ctl_type = 1
        elif "control_string(" in line:
            flds = re.split(",|'", line)[2]
            ctl_type = 4
        elif "control_string_array(" in line:
            flds = re.split(",|'", line)[2]
            ctl_type = 4
        
        if flds is not None:
            if flds in ctl_dict:
                print "Duplicate control file variable - skipping", flds
            else:
                ctl_size = 1

                if ctl_type == 1:
                    # For integers start with zero
                    ctl_values = 0
                    if '_flag' in flds:
                        # Deal with the flag variables
                        ctl_values = 0
                        if flds == 'cbh_check_flag':
                            ctl_values = 1
                        elif flds == 'parameter_check_flag':
                            ctl_values = 1
                    elif 'ON_OFF' in flds:
                        # Deal with variables that turn things on or off
                        ctl_values = 0
                elif ctl_type == 4:
                    # For strings start with an empty string
                    ctl_values = ''
            
                    if flds in defaults:
                        # Variables that have a defined default 
                        ctl_values = defaults[flds]
                        
                    if '_module' in flds:
                        # Add to the module dictionary
                        module_dict[flds] = defaults[flds]
                        
                ctl_dict[flds] = {'Type': ctl_type, 'Size': ctl_size, 'Values': ctl_values}

for kk, vv in module_dict.iteritems():
    # Strip out unneeded control variables
    if kk == 'et_module':
        if vv != 'potet_pm':
            del ctl_dict['humidity_day']
            del ctl_dict['windspeed_day']
        if vv != 'climate_hru':
            del ctl_dict['potet_day']
    elif kk == 'solrad_module':
        if vv != 'climate_hru':
            del ctl_dict['orad_flag']
            del ctl_dict['swrad_day']
    elif kk == 'precip_module':
        if vv != 'climate_hru':
            del ctl_dict['precip_day']
    elif kk == 'temp_module':
        if vv != 'climate_hru':
            del ctl_dict['tmax_day']
            del ctl_dict['tmin_day']
    elif kk == 'transp_module':
            del ctl_dict['transp_day']
            
outfile = 'crap.csv'
outhdl = open(outfile, 'w')
outhdl.write('control_var, control_type, control_value\n')
for kk, vv in ctl_dict.iteritems():
    if vv['Type'] == 1:
        outhdl.write('%s, %d, %d\n' % (kk, vv['Type'], vv['Values']))
    elif vv['Type'] == 4:
        outhdl.write('%s, %d, %s\n' % (kk, vv['Type'], vv['Values']))
    print kk, vv
outhdl.close()
fhdl.close()

# %%
ctl_order = []

for kk, vv in ctl_dict.iteritems():
    ctl_order.append(kk)
    
ctl_order.sort()
print ctl_order

# %%
