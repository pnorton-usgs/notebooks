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
import prms_lib as prms
reload(prms)
import datetime
import numpy as np
import pandas as pd

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/tmp_data'
filename = '%s/test1_luca.statvar' % workdir

# %%
# Load the data
data = prms.statvar(filename).data
data.drop(['rec'], axis=1, inplace=True)

data.head()

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/06191500'
filename = '%s/06191500_hrus.csv' % workdir

infile = open(filename, 'r')
rawdata = infile.read().splitlines()
infile.close()
it = iter(rawdata)

counts = {}
for line in it:
    flds = line.split(',')
    for ff in flds:
        try: 
            int(ff)
        except:
            if ff not in counts:
                counts[ff] = 0
            counts[ff] += 1

for kk, vv in counts.iteritems():
    print '%s: %02d' % (kk, vv)

# %% [markdown]
# ### Get min and max parameter values from final generation of a MOCOM calibration

# %%
# Setup model run information
basedir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t3'
basinid = '06191500'
runid = '2015-03-16_1121'

workdir = '%s/%s/runs/%s' % (basedir, basinid, runid)
mocom_file = '%s/optim_fixed.log' % workdir

# Read in mocom file and use regex to specify variable whitespace between fields
mocom = pd.read_csv(mocom_file, sep=',')
objfcns = [col for col in mocom.columns if 'test' in col]
of_names = {'test0':'NRMSE(monthly)',
            'test1':'NRMSE(annual)',
            'test2':'NRMSE(mean monthly)'}

# --------------------------------------------------------------
# Create an array of colors based on the generation number field
maxgen = max(mocom['gennum'])
mingen = min(mocom['gennum'])


# %%
# Write new parameter file with min and max values derived from
# the last generation in the given optimization run
lastgen = mocom[mocom['gennum'] == maxgen]

pmins = lastgen.min()
pmaxs = lastgen.max()

comb = pd.concat({'max': pmaxs, 'min': pmins}, axis=1)
comb.drop(['setnum', 'soln_num', 'gennum', 'rank'], inplace=True)
comb.drop([col for col in mocom.columns if 'test' in col], inplace=True)
comb.to_csv('param_test.txt', index=True, header=False, sep=' ')

print comb


# %% [markdown]
# ### Compute objective function using high/low data

# %%
def dparse(yr, mo, dy):
    # Date parser for working with the date format from PRMS files

    # Convert to integer first
    yr, mo, dy = [int(x) for x in [yr, mo, dy]]

    dt = datetime.datetime(yr, mo, dy)
    return dt

def nrmse(obsdata, simdata):
    # Compute the Normalized Root Mean Square error
    # NOTE: This routine assumes that both
    #       obsdata and simdata are pandas time series.
    os_diff = obsdata - simdata
    lt_mean = obsdata.mean(axis=0)

    sq_os_diff = os_diff**2
    sq_olt_diff = (obsdata - lt_mean)**2

    return np.sqrt(sq_os_diff.sum() / sq_olt_diff.sum())


# %%
hl_file = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/PRMS_master/06191500/HIGH6.LOW7_subdivide.06191500'
thecols = ['year', 'month', 'day', 'efc', 'obsQ', 'subclass']
df = pd.read_csv(hl_file, sep=r"\s+",
                 header=None, names=thecols,
                 parse_dates={'thedate': ['year', 'month', 'day']},
                 date_parser=dparse, index_col='thedate')
df.head()


# %%
# Get the start and end dates 
st = datetime.datetime(1985,10,1)
en = datetime.datetime(1990,9,30)

# %%
# Open run results
rundir='/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t3/06191500/runs/2015-03-05_1452/01332'
runfile = '%s/default.statvar' % rundir

# Load the data
data = prms.statvar(runfile).data
data.drop(['rec'], axis=1, inplace=True)


# %%
df_sub = df[st:en]
data_sub = data[st:en]

# %%
#combined = pd.merge(df_sub['efc'], data_sub['runoff'], how='left', left_on='thedate', right_on='thedate')
combined = pd.concat([df_sub['efc'], data_sub['runoff'], data_sub['basin_cfs']], axis=1)

print combined.head()
#print df_sub['efc'].head()
#print data_sub['runoff'].head()

# %%
combined[combined['efc'] == 7].iloc[:,1].head()
lows = combined[combined['efc'] == 7][['runoff','basin_cfs']]
highs = combined[combined['efc'] == 6][['runoff', 'basin_cfs']]
print lows.head()
print highs.head()

# %%
nrmse_low = nrmse(lows['runoff'], lows['basin_cfs'])
nrmse_high = nrmse(highs['runoff'], highs['basin_cfs'])
print " Low:", nrmse_low
print "High:", nrmse_high

# %% [markdown]
# <b>Configuration stuff</b>

# %%
configfile = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t4/basin.cfg'
basinConfigFile = 'basin.cfg'

config = ConfigParser.SafeConfigParser()
config.read(configfile)

myconfig = collections.OrderedDict()
#itemsects = {}


#print config.defaults()
for ss in config.sections():
    # Read all other sections
    for ii in config.items(ss):
        #itemsects[ii[0]] = ss
        myconfig[ii[0]] = {}
        myconfig[ii[0]]['sect'] = ss
        if len(ii[1].split(',')) == 1:
            if ii[1].isdigit():
                if float(ii[1]).is_integer():
                    myconfig[ii[0]]['value'] = int(ii[1])
                else:
                    myconfig[ii[0]]['value'] = float(ii[1])
            else:
                myconfig[ii[0]]['value'] = ii[1]
            
        else:
            # TODO: check for float, int, or string types in list
            myconfig[ii[0]]['value'] = ii[1].split(',')

for kk,vv in config.defaults().items():
    # Read the DEFAULT section
    myconfig[kk]['value'] = config.get('DEFAULT', kk)
    myconfig[kk]['sect'] = 'DEFAULT'
    
print myconfig
#print itemsects

# %%
