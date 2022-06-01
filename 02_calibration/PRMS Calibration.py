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
import collections
import ConfigParser
import datetime
import os
import re
import shutil
# import subprocess
import numpy as np
import prms_lib as prms
# reload(prms)


# %%
def to_datetime(date_str):
    """Takes a date string of the form 'YYYY-MM-DD HH:mm:ss' (and variations thereof)
       and converts it to a datetime"""
    return datetime.datetime(*[int(x) for x in re.split('-| |:', date_str)])


# %%
#configfile = 'basin.cfg'
configfile = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t4/basin.cfg'
basinConfigFile = 'basin.cfg'

config = ConfigParser.SafeConfigParser()
config.read(configfile)

print(os.getcwd())
#print config.options('Paths')

# %% [markdown]
# ### Process the configuration file

# %%
# Load path information
base_dir = config.get('DEFAULT', 'base_dir', 0)
base_calib_dir = config.get('DEFAULT', 'base_calib_dir', 0)
template_dir = config.get('DEFAULT', 'template_dir', 0)


# Get the names of the PRMS subdirectories for control, input, and output areas
prms_control_sub = config.get('DEFAULT', 'prms_control_sub', 1)
prms_input_sub = config.get('DEFAULT', 'prms_input_sub', 0)
prms_output_sub = config.get('DEFAULT', 'prms_output_sub', 0)

# Get dates for statistics
start_date_model = config.get('General', 'start_date_model', 0)
start_date = config.get('General', 'start_date', 0)
end_date = config.get('General', 'end_date', 0)
obs_var = config.get('General', 'obs_var', 0)
sim_var = config.get('General', 'sim_var', 0)

# Get filename containing basins to calibrate
basins_xref_file = config.get('Files', 'basins_xref_file', 0)
basins_file = config.get('Files', 'basins_file', 0)

# Parameter limits file
param_range_file = config.get('Files', 'param_range_file', 1)

# margins file
test_func_margin_file = config.get('Files', 'test_func_margin_file', 1)

# PRMS files
prms_control_file = config.get('Files', 'prms_control_file')
prms_input_file = config.get('Files', 'prms_input_file')

# Default parameters ranges for PRMS
param_list_file = config.get('Files', 'param_list_file')

# MOCOM program
mocom_cmd = config.get('Commands', 'mocom')

# PRMS program
prms_cmd = config.get('Commands', 'prms', 0)

# Calibration script
calib_run = config.get('Commands', 'calib_run', 0)

stats_script = config.get('Commands', 'stats_script', 0)
plot_script = config.get('Commands', 'plot_script', 0)
cp_cmd = config.get('Commands', 'cp', 0)

# Log file for mocom runs
log_file = config.get('Files', 'log_file', 1)

# Calibration settings
nstart = config.getint('Calibration', 'nstart')
nsets = config.getint('Calibration', 'nsets')
nparam = config.getint('Calibration', 'nparam')
ntests = config.getint('Calibration', 'ntests')

# Objective function information
ofs = config.get('ObjFcns', 'fcn').split(',')
interval = config.get('ObjFcns', 'interval').split(',')

# %%
# Could do testing of config parameters to make sure everything exists
print('\n', '-' * 5, 'Dates', '-' * 5)
print('  start_date_model:', start_date_model)
print('start_date (calib):', start_date)
print('          end_date:', end_date)
print('\n', '-' * 5, 'Paths', '-' * 5)
print('      base_dir:', base_dir)
print('base_calib_dir:', base_calib_dir)
print('  template_dir:', template_dir)
print('PRMS control subdirectory:', prms_control_sub)
print('  PRMS input subdirectory:', prms_input_sub)
print(' PRMS output subdirectory:', prms_output_sub)
print('\n', '-' * 5, 'Files', '-' * 5)
print('  input parameter file:', prms_input_file)
print('control parameter file:', prms_control_file)
print('           basins_file:', basins_file)
print('      param_range_file:', param_range_file)
print(' test_func_margin_file:', test_func_margin_file)
print('  calibration log file:', log_file)
print('\n', '-' * 5, 'Scripts', '-' * 5)
print('      mocom executable:', mocom_cmd)
print('       prms executable:', prms_cmd)
print('calibration run script:', calib_run)
print('\n', '-' * 5, 'Calibration settings', '-' * 5)
print('nstart:', nstart)
print(' nsets:', nsets)
print('nparam:', nparam)
print('ntests:', ntests)

# %% [markdown]
# <b>If necessary create basins_file from special lookup file</b>

# %%
# On yeti we'll use a pre-defined cross-reference file to create a default basins_file
# This cross references between station id and the lumen-derived directory.
xref = pd.read_csv(basins_xref_file)
xref['lumen_name'].to_csv(basins_file, index=False)

# %% [markdown]
# <b>Open basins_file and set the runid</b>

# %%
# Read the basins_file. Single basin per line
bfile = open(basins_file, "r")
basins = bfile.read().splitlines()
bfile.close()
print basins

# Define a run id
runid = datetime.datetime.now().strftime('%Y-%m-%d_%H%M')
print 'runid:', runid

# %%
# Create parameter default ranges file from from PRMS -print results

#infile = open('%s/default.control.par_name' % template_dir, 'r')
infile = open(param_list_file, 'r')

rawdata = infile.read().splitlines()
infile.close()

it = iter(rawdata)

for line in it:
    if line == '--------------- PARAMETERS ---------------':
        break
param_dict = {}

for line in it:
    flds = line.split(':')
    
    if len(flds) < 2:
        continue
    
    key = flds[0].strip()
    val = flds[1].strip()
    
    if key == 'Name':
        cparam = val
        param_dict[cparam] = {}
    else:
        param_dict[cparam][key] = val  
    
#outfile = open(param_list_file, 'w')
outfile = open('param_ranges.txt', 'w')
outfile.write('parameter max min\n')

for kk, vv in param_dict.iteritems():
    outfile.write('%s %f %f\n' % (kk, float(vv['Max']), float(vv['Min'])))
    
outfile.close()  


# %%
filelist = [test_func_margin_file]

# For each basin create the directory structure and populate it
for bb in basins:
    print '-'*60
    print 'Basin: %s' % bb
    cbasin = '%s/%s' % (base_calib_dir, bb)  # directory for the calibration run of this basin
    csrc = '%s/%s' % (template_dir, bb)   # source files directory for this basin
    cdirs = [cbasin, '%s/runs' % cbasin, '%s/pdf' % cbasin]
    
    for dd in cdirs:
        if not os.path.exists(dd):
            # Make the basin directory (bb)
            try:
                #print "\tCreating %s" % dd
                os.makedirs(dd)
            except OSError, err:
                print "\tError creating directory: %s" % err
                exit(1)

    # Get a list of input files from the control file.
    # NOTE: Assumes sane paths to input files in the control file
    ctlfile = '%s/%s' % (csrc, prms_control_file) #currently doesn't handle subdir for ctl file
    src_config = ctlfile+'\n'    
    
    paramfile = '%s/%s' % (csrc, prms_input_file)
    src_config += paramfile+'\n'
    
    ctl = prms.control(ctlfile)
    input_file_list = ['tmax_day', 'tmin_day', 'precip_day', 'data_file']
    
    src_files = ''
    for vv in input_file_list:
        #print vv, ctl.get_var(vv)
        src_files += '%s/%s\n' % (csrc, ctl.get_var(vv)[0]['values'][0])
    
    print '-'*20
    print src_files
    print '-'*20                
    
    # ==================================================================    
    # Verify the date ranges for the basin from the streamflow.data file
    streamflow_file = '%s/%s' % (csrc, ctl.get_var('data_file')[0]['values'][0])
    streamflow = prms.streamflow(streamflow_file).data
    streamflow.dropna(axis=0, how='any', inplace=True)
    first_date = streamflow[streamflow.notnull()].index.min()
    last_date = streamflow[streamflow.notnull()].index.max()

    if first_date > to_datetime(start_date_model):
        print "\t* Adjusting start_date_model and start_date to reflect available streamflow data"
        start_date_model = first_date.strftime('%Y-%m-%d')
        yr = int(first_date.strftime('%Y'))
        mo = int(first_date.strftime('%m'))
        dy = int(first_date.strftime('%d'))
        
        # Set calibration start date to 2 years after model start date
        start_date = '%d-%d-%d' % (yr+2 ,mo, dy)
    if last_date < to_datetime(end_date):
        print "\t* Adjusting end_date to reflect available streamflow data"
        end_date = last_date.strftime('%Y-%m-%d')
        
    print "\tstart_date_model:", start_date_model
    print "\t      start_date:", start_date
    print "\t        end_date:", end_date
    # ------------------------------------------------------------------
    
    
    # ==================================================================
    # Create the param_range_file from the sensitive parameters file and the default ranges
    
    # Some parameters should always be excluded even if they're sensitive
    exclude_params = ['dday_slope', 'dday_intcp', 'jh_coef_hru', 'jh_coef']
    
    # Some parameter should always be included
    include_params = ['tmin_cbh_adj', 'tmax_cbh_adj', 'rain_cbh_adj', 'snow_cbh_adj']

    try:
        sensparams_file = open('%s/hruSens.csv' % csrc, 'r')
    except IOError:
        print "\tERROR: Missing hruSens.csv file for %s... skipping" % bb
        continue
        
    default_rng_file = open('%s/param_ranges.txt' % base_dir, 'r')
    
    # Read in the default parameter ranges
    raw_range = default_rng_file.read().splitlines()
    default_rng_file.close()
    it = iter(raw_range)
    
    def_ranges = {}
    it.next()
    for line in it:
        flds = line.split(' ')
        def_ranges[flds[0]] = {'max': float(flds[1]), 'min': float(flds[2])}
        
    # Read in the sensitive parameters
    rawdata = sensparams_file.read().splitlines()
    sensparams_file.close()
    it = iter(rawdata)

    counts = {}
    for line in it:
        flds = line.split(',')
        for ff in flds:
            try: 
                int(ff)
            except:
                if ff not in exclude_params:
                    if ff not in counts:
                        counts[ff] = 0
                    counts[ff] += 1
    
    # Add in the include_params if they are missing from the sensitive parameter list
    for pp in include_params:
        if pp not in counts:
            counts[pp] = 0
    
    # ---------------------------------------------------------------
    # Adjust the min/max values for each parameter based on the 
    # initial values from the .params file
    src_params = prms.parameters(paramfile)
                     
    # Write the param_list file
    outfile = open('%s/%s' % (cbasin, param_range_file), 'w')
    for kk, vv in counts.iteritems():
        # Grab the current param (kk) from the .params file and verify the
        # upper and lower bounds. Modify them if necessary.
        src_vals = src_params.get_var(kk)['values']
        src_mean = np.mean(src_vals)
        src_min = np.min(src_vals)
        src_max = np.max(src_vals)
        
        # Set upper and lower bounds
        user_min = def_ranges[kk]['min']
        user_max = def_ranges[kk]['max']
        if user_min > src_min:
            user_min = src_min
        if user_max < src_max:
            user_max = src_max
        
        C = abs(user_min) + 10.
        
        adjMin = ((user_min + C) * (src_mean + C) / (src_min + C)) - C
        adjMax = ((user_max + C) * (src_mean + C) / (src_max + C)) - C
                
        if round(adjMin, 5) != round(def_ranges[kk]['min'], 5):
            print '\t%s: lower bound adjusted (%f to %f)' % (kk, def_ranges[kk]['min'], adjMin)
        if round(adjMax, 5) != round(def_ranges[kk]['max'], 5):
            print '\t%s: upper bound adjusted (%f to %f)' % (kk, def_ranges[kk]['max'], adjMax)
        
        outfile.write('%s %f %f\n' % (kk, adjMax, adjMin))
        #outfile.write('%s %f %f\n' % (kk, def_ranges[kk]['max'], def_ranges[kk]['min']))
    outfile.close()
    # ---------------------------------------------------------------
    
    
    # Update the nparam and nsets values for the number of parameters we are calibrating
    nparam = len(counts)
    nsets = int(nparam * 2.5)    # Can't have the number of sets too small or too big
    ntests = len(ofs)            # Set the number of tests according the size of ofs                
                                

    # Copy files to basin directory
    #     param_range_file
    #     test_func_margin_file
    for ff in filelist:
        #print "\tCopying %s" % ff
        shutil.copyfile('%s/%s' % (base_calib_dir, ff), '%s/%s' % (cbasin, ff))
            
        
    # Change into the basin directory
    orig_path = os.getcwd()
    os.chdir(cbasin)

    # ***************************************************************************
    # Write basin configuration file

    # Create a config file for the run
    newconfig = ConfigParser.SafeConfigParser()
    newconfig.set('DEFAULT', 'base_dir', base_dir)
    newconfig.set('DEFAULT', 'base_calib_dir', base_calib_dir)
    newconfig.set('DEFAULT', 'template_dir', template_dir)
    newconfig.set('DEFAULT', 'basin', bb)

    newconfig.add_section('General')
    newconfig.set('General', 'start_date_model', start_date_model)
    newconfig.set('General', 'start_date', start_date)
    newconfig.set('General', 'end_date', end_date)
    newconfig.set('General', 'obs_var', obs_var)
    newconfig.set('General', 'sim_var', sim_var)
    
    newconfig.add_section('Paths')
    newconfig.set('Paths', 'runs_sub', 'runs')
    newconfig.set('Paths', 'prms_control_sub', prms_control_sub)
    newconfig.set('Paths', 'prms_input_sub', prms_input_sub)
    newconfig.set('Paths', 'prms_output_sub', prms_output_sub)
    
    newconfig.add_section('Sources')
    newconfig.set('Sources', 'source_data', src_files)
    newconfig.set('Sources', 'source_config', src_config)
    
    newconfig.add_section('Files')
    newconfig.set('Files', 'param_range_file', param_range_file)
    newconfig.set('Files', 'test_func_margin_file', test_func_margin_file)
    newconfig.set('Files', 'prms_control_file', prms_control_file)
    newconfig.set('Files', 'prms_input_file', prms_input_file)

    newconfig.add_section('Commands')
    newconfig.set('Commands', 'prms', prms_cmd)
    newconfig.set('Commands', 'stats_script', stats_script)
    newconfig.set('Commands', 'plot_script', plot_script)
    
    newconfig.set('Commands', 'cp', cp_cmd)
    
    newconfig.add_section('ObjFcns')
    newconfig.set('ObjFcns', 'fcn', ','.join(ofs))
    newconfig.set('ObjFcns', 'interval', ','.join(interval))
    
    with open(basinConfigFile, 'w') as outfile:
        newconfig.write(outfile)

    # Run MOCOM with the following arguments:
    #     nstart
    #     nsets
    #     nparam
    #     ntests
    #     calib_run
    #     run_id
    #     log_file
    #     param_range_file
    #     test_func_margin_file
    #print '\tRunning MOCOM'
    cmd_opts = ' %d %d %d %d %s %s %s %s %s' % (nstart, nsets, nparam, ntests, 
                                                calib_run, runid, log_file, 
                                                param_range_file, test_func_margin_file)
    print "\tRunning MOCOM..."
    print mocom_cmd + cmd_opts
    #subprocess.call(mocom + cmd_opts, shell=True)
    
    os.chdir(orig_path)


# %%
os.chdir(orig_path)

# %%
os.chdir('/Users/pnorton/Projects/National_Hydrology_Model/src/prms_calibration')
os.getcwd()

# %%
os.getcwd()

# %%
type(cdirs)

# %%
