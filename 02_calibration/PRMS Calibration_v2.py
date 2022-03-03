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
import collections
import ConfigParser
import datetime
import os
import re
import shutil
import subprocess
import numpy as np
import prms_lib as prms
import prms_cfg_yaml as prms_cfg
from prms_calib_helpers import read_default_params, read_sens_params, adjust_param_ranges
reload(prms)
reload(prms_cfg)

# %%

# %%
#configfile = 'basin.cfg'
configfile = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t6/basin.cfg'
basinConfigFile = 'basin.cfg'

cfg = prms_cfg.cfg(configfile)
print os.getcwd()

# %% [markdown]
# ### Process the configuration file

# %%
#cfg.list

# %% [markdown]
# <b>If necessary create basins_file from special lookup file</b>

# %%
# On yeti we'll use a pre-defined cross-reference file to create a default basins_file
# This cross references between station id and the lumen-derived directory.
xref = pd.read_csv(cfg['basins_xref_file']['value'])
xref['lumen_name'].to_csv(basins_file, index=False)

# %% [markdown]
# <b>Open basins_file and set the runid</b>

# %%
# Read the basins_file. Single basin per line
bfile = open(cfg.get_value('basins_file'), "r")
basins = bfile.read().splitlines()
bfile.close()
print basins

# Define a run id
runid = datetime.datetime.now().strftime('%Y-%m-%d_%H%M')
print 'runid:', runid

# %%
filelist = [cfg.get_value('test_func_margin_file')]

# For each basin create the directory structure and populate it
for bb in basins:
    print '-'*60
    print 'Basin: %s' % bb
    
    # Add the current basin to the configuration dictionary
    cfg.add_var('basin', bb)
    
    cbasin = '%s/%s' % (cfg.get_value('base_calib_dir'), bb)  # directory for the calibration run of this basin
    csrc = '%s/%s' % (cfg.get_value('template_dir'), bb)   # source files directory for this basin
    cdirs = [cbasin, '%s/runs' % cbasin, '%s/pdf' % cbasin]
    
    for dd in cdirs:
        if not os.path.exists(dd):
            # Make the basin directory (bb)
            try:
                os.makedirs(dd)
            except OSError, err:
                print "\tError creating directory: %s" % err
                exit(1)

    # Get a list of input files from the control file.
    # NOTE: Assumes sane paths to input files in the control file
    src_config = []
    ctlfile = '%s/%s' % (csrc, cfg.get_value('prms_control_file')) #currently doesn't handle subdir for ctl file
    src_config.append(ctlfile)
    
    paramfile = '%s/%s' % (csrc, cfg.get_value('prms_input_file'))
    src_config.append(paramfile)
    
    # Read the PRMS control file
    ctl = prms.control(ctlfile)
    input_file_list = ['tmax_day', 'tmin_day', 'precip_day', 'data_file']
    
    src_files = []
    # Build list of data files from what is found in the control file
    # These files will be linked into the run directory
    for vv in input_file_list:
        src_files.append('%s/%s' % (csrc, ctl.get_var(vv)['values'][0]))
    
    # Add any subdivide files and observed range/value files that are specified for the objective functions
    of_dict = cfg.get_value('objfcn')
    for kk, vv in of_dict.iteritems():
        for k2, v2 in vv.iteritems():
            # For each objective function look for sd_file and obs_file
            if k2 in ['sd_file', 'obs_file']:
                cfile = '%s/%s' % (csrc, v2)
                if v2 is not None and cfile not in src_files:
                    src_files.append(cfile)
    
    
    cfg.add_var('source_config', src_config)
    cfg.add_var('source_data', src_files)
    
    print '-'*20
    print src_files
    print '-'*20                
    
    # ==================================================================    
    # Verify the date ranges for the basin from the streamflow.data file
    streamflow_file = '%s/%s' % (csrc, ctl.get_var('data_file')['values'][0])
    first_date, last_date = prms.streamflow(streamflow_file).date_range

    if first_date > prms.to_datetime(cfg.get_value('start_date_model')):
        print "\t* Adjusting start_date_model and start_date to reflect available streamflow data"
        cfg.update_value('start_date_model', first_date.strftime('%Y-%m-%d'))
        yr = int(first_date.strftime('%Y'))
        mo = int(first_date.strftime('%m'))
        dy = int(first_date.strftime('%d'))
        
        # Set calibration start date to 2 years after model start date
        cfg.update_value('start_date', '%d-%d-%d' % (yr+2 ,mo, dy))
        
    if last_date < prms.to_datetime(cfg.get_value('end_date')):
        print "\t* Adjusting end_date to reflect available streamflow data"
        cfg.update_value('end_date', last_date.strftime('%Y-%m-%d'))
        
    print "\tstart_date_model:", cfg.get_value('start_date_model')
    print "\t      start_date:", cfg.get_value('start_date')
    print "\t        end_date:", cfg.get_value('end_date')
    # ------------------------------------------------------------------
    
    
    # ==================================================================
    # Create the param_range_file from the sensitive parameters file 
    # and the default ranges
    
    # Some parameters should always be excluded even if they're sensitive
    exclude_params = ['dday_slope', 'dday_intcp', 'jh_coef_hru', 'jh_coef']
    
    # Some parameter should always be included
    #include_params = ['tmin_cbh_adj', 'tmax_cbh_adj', 'rain_cbh_adj', 'snow_cbh_adj']
    include_params = []
    
    default_rng_filename = '%s/%s' % (cbasin, cfg.get_value('default_param_list_file'))
    if not os.path.isfile(default_rng_filename):
        # We don't have a list of default ranges yet, so create it
        if not prms.create_default_range_file('%s.par_name' % ctlfile, default_rng_filename):
            # Only happens when input file can't be read or output file can't be written
            # Assume the input file is not created yet.
            hld = os.getcwd()
            os.chdir(csrc)
            cmd_opts = " -C%s -set param_file %s -print" % (cfg.get_value('prms_control_file'),
                                                            cfg.get_value('prms_input_file'))
            subprocess.call(cfg.get_value('cmd_prms') + cmd_opts, shell=True)
            os.chdir(hld)
            
            # Try to get the default ranges one more time
            if not prms.create_default_range_file('%s.par_name' % ctlfile, default_rng_filename):
                raise IOError
    
    # Read the default parameter ranges in
    def_ranges = read_default_params(default_rng_filename)
    
    # Read the sensitive parameters in
    sens_params = read_sens_params('%s/hruSens.csv' % csrc, include_params, exclude_params)
    # ---------------------------------------------------------------    
    
    # ==================================================================
    # Adjust the min/max values for each parameter based on the 
    # initial values from the input parameter file
    adjust_param_ranges(paramfile, sens_params, def_ranges, 
                        '%s/%s' % (cbasin, cfg.get_value('param_range_file')))
    
    
    # Update the nparam and nsets values for the number of parameters we are calibrating
    cfg.update_value('nparam', len(sens_params))
    cfg.update_value('nsets', int(cfg.get_value('nparam') * 2.5)) # Can't have the number of sets too small or too big
    cfg.update_value('ntests', len(cfg.get_value('of_link')))  # Set the number of tests according the size of ofs               

    # Copy files to basin directory
    print 'filelist:', filelist
    for ff in filelist:
        shutil.copyfile('%s/%s' % (cfg.get_value('base_calib_dir'), ff), '%s/%s' % (cbasin, ff))
                
    # Change into the basin directory
    orig_path = os.getcwd()
    os.chdir(cbasin)

    # ***************************************************************************
    # Write basin configuration file for run
    cfg.write_config(basinConfigFile)
    # ---------------------------------------------------------------------------
    
    # Run MOCOM with the following arguments:
    #     nstart, nsets, nparam, ntests, calib_run, run_id, log_file, param_range_file, test_func_margin_file
    cmd_opts = ' %d %d %d %d %s %s %s %s %s' % (cfg.get_value('nstart'), cfg.get_value('nsets'), cfg.get_value('nparam'),
                                                cfg.get_value('ntests'), cfg.get_value('calib_run'), runid, 
                                                cfg.get_value('log_file'), cfg.get_value('param_range_file'), 
                                                cfg.get_value('test_func_margin_file'))
    print "\tRunning MOCOM..."
    print cfg.get_value('cmd_mocom') + cmd_opts
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
