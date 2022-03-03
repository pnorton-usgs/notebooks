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
import prms_lib as prms
import prms_objfcn as objfcn
import prms_cfg
import ConfigParser
import re
import pandas as pd
import numpy as np
import datetime


# %%
def to_datetime(date_str):
    """Takes a date string of the form 'YYYY-MM-DD HH:mm:ss' (and variations thereof)
       and converts it to a datetime"""
    return datetime.datetime(*[int(x) for x in re.split('-| |:', date_str)])


# %%
# Setup model run information
#basinids = ['1.11_10173450', '1.141_09223000', '1.181_01193500', '1.183_06469400',
#            '1.20_03574500', '1.20_08269000', '1.9_11468000', '2.121_14301000',
#            '2.140_02016000', '2.18_05131500', '2.3_09510200', '2.67_04185000',
#            '2.74_05584500', '3.45_07249985', '4.39_03410500', '4.9_06191500']
basinids = ['4.39_03410500']
runid = '2015-04-29_1424'

basedir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t4'
templatedir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/gage_master'
statvar_file = 'daymet.statvar'

# %%
for basinid in basinids:
    # Read the config file
    configfile = '%s/%s/basin.cfg' % (basedir, basinid)
    
    config = ConfigParser.SafeConfigParser()
    config.read(configfile)
    cfg = prms_cfg.cfg_get_vars(config)
    
    templatedir = cfg['template_dir']['value']
    runs_sub = cfg['runs_sub']['value']
    workdir = '%s/%s/%s/%s' % (basedir, basinid, runs_sub, runid) 
    modeldir = '%s/%s' % (templatedir, basinid)
    
    sim_var = cfg['sim_var']['value']
    obs_var = cfg['obs_var']['value']
    st = to_datetime(cfg['start_date']['value'])
    en = to_datetime(cfg['end_date']['value'])
    
    mocom_file = '%s/optim_fixed.log' % workdir

    # Read in mocom file and use regex to specify variable white space between fields
    mocom = pd.read_csv(mocom_file, sep=',')
    maxgen = max(mocom['gennum'])    # Get the number of the last generation
    #print "Last generation = %d" % maxgen
    
    # Get list of solutions from the last generation in the optimization log
    modelrunids = mocom['soln_num'].loc[mocom['gennum'] == maxgen].tolist()
    
    results = []
    
    for rr in modelrunids:
        cfile = '%s/%05d/%s' % (workdir, rr, statvar_file)
        
        sv = prms.statvar(cfile)
        tmp_df = sv.data[st:en]
    
        ns = objfcn.compute_objfcn('NS', 'daily', tmp_df, obs_var, sim_var) * -1.
        results.append(ns)
        #print 'Set: % 4d\tNS = %0.5f' % (rr, ns)
    
    outdata = {'set': modelrunids, 
               'NS': results}
    df = pd.DataFrame(outdata)
    
    best_modelrunid = '%05d' % df.sort(['NS','set']).iloc[-1].set
    best_modelrunNS = df.sort(['NS', 'set']).iloc[-1].NS
    
    # Output the top NS value
    print '%s: %0.4f (%s)' % (basinid, best_modelrunNS, best_modelrunid)

# %%
df.sort(['NS','set']).iloc[-1].NS

# %%
df.sort(['NS','set'])

# %%
