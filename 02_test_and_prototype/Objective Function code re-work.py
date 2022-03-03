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
import pandas as pd
import numpy as np
import collections
import datetime
import re
import os

import prms_lib as prms
import prms_cfg
reload(prms_cfg)
reload(prms)


# %%
def dparse(*args):
#def dparse(yr, mo, dy, hr, minute, sec):
    # Date parser for working with the date format from PRMS files
    print args
    return datetime.datetime(*[int(x) for x in args])

    #dt = datetime.datetime(yr, mo, dy, hr, minute, sec)
    #return dt

def to_datetime(date_str):
    """Takes a date string of the form 'YYYY-MM-DD HH:mm:ss' (and variations thereof)
       and converts it to a datetime"""
    return datetime.datetime(*[int(x) for x in re.split('-| |:', date_str)])


# %%
# Setup model run information
basedir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t4'
basinid = '4.39_03410500'
configfile = '%s/%s/basin.cfg' % (basedir, basinid)

runid = '2015-04-29_1424'

cfg = prms_cfg.cfg(configfile)

tst_master_dir = '%s/%s/from_lauren' % (cfg.get_value('template_dir'), basinid)
tmp_obsfile = '%s/NWIS.error' % (tst_master_dir)
tmp_subfile = '%s/EFC.high6low7' % (tst_master_dir)

print 'testing master directory:', tst_master_dir
print 'testing observation ranges:', tmp_obsfile

objfcn_dict = collections.OrderedDict(OF1= dict(objfcn='NRMSE', of_intv ='daily',
                                                obstype='range', obsfile=tmp_obsfile, obs_intv='daily',
                                                sdval=6, sdfile=tmp_subfile,  
                                                obsvar='runoff', simvar='basin_cfs'))
orig_dir = os.getcwd()
os.chdir('%s/%s/runs/%s/00540' % (basedir, basinid, runid))
ctl = prms.control(cfg.get_value('prms_control_file'))

# Load the statvar file
df_sv = pd.DataFrame(prms.statvar(ctl.get_var('stat_var_file')['values'][0]).data['basin_cfs'])
print df_sv.head()

# %% [markdown]
# <b>Read the observation values or ranges from the obsfile</b>

# %%
# Equate objfcn values to number of columns expected
#collookup = {'range': 2, 'value': 1, 'daily': 3, 'monthly': 2, 'annual': 1, 'mnmonth': 1}
colnm_lookup = {'range': ['obs_lower', 'obs_upper'],
                'value': ['obs_val'],
                'daily': ['year', 'month', 'day'],
                'monthly': ['year', 'month'],
                'annual': ['year'],
                'mnmonth': ['month']}

# Range files from Lauren use -99.0 as missing, other files use -999.0
missing = [-999.0, -99.0]

for kk, vv in objfcn_dict.iteritems():
    print "Objective Function:", kk
    
    # Get the total number of columns for the dtype and obs_intv
    # and build the names to use for the dataframe.
    thecols = []
    #tcol = collookup[vv['obstype']] + collookup[vv['obs_intv']]
    thecols.extend(colnm_lookup[vv['obs_intv']])
    thecols.extend(colnm_lookup[vv['obstype']])
        
    # Read in the observation values/ranges
    if vv['obs_intv'] == 'mnmonth':
        # The index won't be a datetime, instead it's a month value
        df1 = pd.read_csv(vv['obsfile'], sep=r"\s*", engine='python', usecols=range(0,len(thecols)), 
                          header=None, na_values=missing, names=thecols, index_col=0)
    else:
        # NOTE: When parsing year-month dates pandas defaults to the 21st of each month. I'm not sure yet
        #       if this will cause a problem.
        #       Annual dates are parsed as Jan-1 of the given year.
        # TODO: if 'obsfile' == statvar then read the observed values in from the statvar file
        df1 = pd.read_csv(vv['obsfile'], sep=r"\s*", engine='python', usecols=range(0,len(thecols)), 
                          header=None, na_values=missing,
                          names=thecols, parse_dates={'thedate': colnm_lookup[vv['obs_intv']]}, index_col='thedate')
    print '----- Observed data -----'
    print df1.head()
    print '-'*30
    
    # Read in the subdivide data, if specified
    if len(vv['sdfile']) > 0:
        # The subdivide file must be a daily timestep
        thecols = ['year', 'month', 'day', 'sdval']
        
        df2 = pd.read_csv(vv['sdfile'], sep=r"\s*", engine='python', usecols=range(0,len(thecols)), 
                          header=None, na_values=missing,
                          names=thecols, parse_dates={'thedate': ['year', 'month', 'day']}, index_col='thedate')
        
        print '----- Subdivide data -----'
        print df2.head()
        print '-'*30
    
        # Merge the subdivide data with the observed data
        if vv['obs_intv'] != 'daily':
            # The observed data is not a daily timestep (subdivide data is daily) so raise an error. 
            print 'ERROR: observed data must be daily timestep when using subdivide data'
            pass
        
        df1_j1 = df1.join(df2, how='left')
        
        print '----- Observed data merged with subdivide data -----'
        print df1_j1.head()
        print '-'*30
        
        df_final = df_sv.join(df1_j1, how='left')
        print '----- statvar merged with observed data -----'
        print df_final.head()
        print '-'*30
    
        # Subset to only include values which match 'sdval'
        df_final = df_final[df_final['sdval'] == vv['sdval']]
        print '----- subset final to sdval -----'
        print df_final.head()
        print '-'*30
    
    # Now resample to specified of_intv
    if vv['of_intv'] == 'monthly':
        df_final = df_final.resample('M', how='mean')
    elif vv['of_intv'] == 'annual':
        # TODO: For now the annual interval is assumed to be water-year based
        df_final = df_final.resample('A-SEP', how='mean')
    elif vv['of_intv'] == 'mnmonth':
        monthly = df_final.resample('M', how='mean')
        df_final = monthly.resample('M', how='mean').groupby(monthly.index.month).mean()

    # TODO: strip rows with NaN observations out of dataframe
    df_final.dropna(axis=0, how='any', thresh=None, subset=colnm_lookup[vv['obstype']], inplace=True)

# %%
df_final.tail()

# %%
# Check the time interval of the obsfile using the first two dates
tst = df1.index.tolist()
print (tst[1] - tst[0]).days

# %%
tst_dt = [1980,12,14]
dparse(1980,12,14,11,59,03, 19)

# %%
