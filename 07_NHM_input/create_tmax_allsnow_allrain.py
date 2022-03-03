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
# %matplotlib inline

import datetime
import pandas as pd
import numpy as np
from prms_objfcn import dparse

import matplotlib.pyplot as plt

# %%
region = '10U'
base_daymet_dir = '/media/scratch/PRMS/datasets/daymet'
base_snodas_dir = '/media/scratch/PRMS/datasets/snodasPrecip'

tmaxfile = '%s/r%s/daymet_1980_2011_tmax.cbh' % (base_daymet_dir, region)
prcpfile = '%s/r%s/daymet_1980_2011_prcp.cbh' % (base_daymet_dir, region)

#tmaxfile = '/media/scratch/PRMS/datasets/daymet/r10U_daymet_2004-2010_tmax.cbh'
#prcpfile = '/media/scratch/PRMS/datasets/daymet/r10U_daymet_2004-2010_prcp.cbh'
#vpfile = '/media/scratch/PRMS/datasets/daymet/r10u_daymet_vp_2004-01-01_2010-12-31.csv'

event_type = 'snow' # one of snow, liquid
prcpmaskfile = '%s/r%s/%s_only_SNODASPRCP_Daily_r%s_byHRU_2004-01-01_2014-12-31.csv' % (base_snodas_dir, region, event_type, region)
#prcpmaskfile = '%s/%s_only_SNODASPRCP_Daily_r10U_byHRU_2004-01-01_2014-12-31.csv' % (base_snodas_dir, event_type)

# See Aiguo Dai, 2008
rain_max = 38.    # based on 90% probability of rain
rain_min = 34.7   # based on 50% 
snow_max = 34.25  # based on 50% probability of snow
snow_min = 28.2   # based on 95%

st = datetime.datetime(2004,1,1)
en = datetime.datetime(2010,12,31)


# %%
def read_cbh(filename):
    # Read a CBH file
    missing = [-99.0, -999.0]

    infile = open(filename, 'r')
    fheader = ''

    for ii in range(0,3):
        line = infile.readline()

        if line[0:4] in ['prcp', 'tmax', 'tmin']:
            # Change the number of HRUs included to one
            numhru = int(line[5:])
            fheader += line[0:5] + ' 1\n'
        else:
            fheader += line
    print fheader
    print 'numhru:', numhru

    df1 = pd.read_csv(infile, sep=' ', skipinitialspace=True, header=None, 
                      date_parser=dparse, parse_dates={'thedate': [0,1,2,3,4,5]}, index_col='thedate')

    infile.close()

    # Renumber/rename columns to reflect HRU number
    df1.rename(columns=lambda x: df1.columns.get_loc(x)+1, inplace=True)
    return df1



# %% [markdown]
# ## Read the daymet tmax and prcp cbh file

# %%
# Read the tmax CBH file and restrict to the given date range
df_tmax = read_cbh(tmaxfile)
df_tmax = df_tmax[st:en]

# Read the prcp CBH file and restrict to the given date range
df_prcp = read_cbh(prcpfile)
df_prcp = df_prcp[st:en]

# Create mask from precip values
df_prcp[df_prcp == 0.0] = np.nan
df_prcp[df_prcp > 0.0] = 1

# %%
# TRASH me
df_prcp.head()

# %% [markdown]
# ## Read SNODAS precip event mask file

# %%
# Read the SNODAS precip event mask
df_mask = pd.read_csv(prcpmaskfile, na_values=['MAXIMUM(none)', 255.0], header=0, skiprows=[0,2])
# df_mask.drop(0, inplace=True)

df_mask.rename(columns={df_mask.columns[0]:'thedate'}, inplace=True)
df_mask['thedate'] = pd.to_datetime(df_mask['thedate'])
df_mask.set_index('thedate', inplace=True)
df_mask = df_mask[st:en]

# Rename the columns to something a little more descriptive
#df2.rename(columns=lambda x: 'HRU_%s' % x, inplace=True)

# Fill in any missing days with NAN
df_mask = df_mask.resample('D', how='max')
df_mask[df_mask > 0.0] = 1

# %%
df_mask.head(20)

# %% [markdown]
# ## Mask tmax to precip events (rain, snow, or mixed)

# %%
# Use numpy (np) to multiply df1 and df3. 
# For some reason pandas multiply is horribly broken
df_tmax_m1 = pd.np.multiply(df_tmax,df_mask)
df_tmax_m1[df_tmax_m1 == 0.0] = np.nan

# Create secondary mask 
df_tmax_m2 = pd.np.multiply(df_tmax, df_prcp)

if event_type == 'snow':
    # Mask out temperature values outside the range we need for both masked sets
    df_tmax_m1[df_tmax_m1 > snow_max] = np.nan
    df_tmax_m1[df_tmax_m1 < snow_min] = np.nan
    
    df_tmax_m2[df_tmax_m2 > snow_max] = np.nan
    df_tmax_m2[df_tmax_m2 < snow_min] = np.nan
elif event_type == 'liquid':
    df_tmax_m1[df_tmax_m1 > rain_max] = np.nan
    df_tmax_m1[df_tmax_m1 < rain_min] = np.nan
    
    df_tmax_m2[df_tmax_m2 > rain_max] = np.nan
    df_tmax_m2[df_tmax_m2 < rain_min] = np.nan
    

# %%
#TRASH me
df_tmax_m2.head()

# %% [markdown]
# ## Compute monthly and mean monthly dataframes for both masked sets

# %%
# Resample to monthly mean
df_tmax_m1_mon = df_tmax_m1.resample('M', how='mean')

# compute mean monthly tmax for each hru
df_tmax_m1_mnmon = df_tmax_m1_mon.groupby(df_tmax_m1_mon.index.month).mean()

# Do the same for the secondary masked set
df_tmax_m2_mon = df_tmax_m2.resample('M', how='mean')
df_tmax_m2_mnmon = df_tmax_m2_mon.groupby(df_tmax_m2_mon.index.month).mean()

#plt.figure(figsize=(20,30))
#df6.iloc[:,0].plot(kind='hist')

print 'tmax_m1'
print 'Max:', df_tmax_m1_mnmon.max().max()
print 'Min:', df_tmax_m1_mnmon.min().min()

print 'tmax_m2'
print 'Max:', df_tmax_m2_mnmon.max().max()
print 'Min:', df_tmax_m2_mnmon.min().min()

# %%
# Where tmax_m1_mnmon is NaN replace with value from tmax_m2_mnmon if it exists
# otherwise replace with default value
# print df_tmax_m1_mnmon.head(12).iloc[:,10289]
# print df_tmax_m2_mnmon.head(12).iloc[:,10289]

# First replace any Nan entries with a valid entry from df_tmax_m2_mnmon
# Then replace any remaining NaN entries by padding them with the last valid value
df_tmax_final = df_tmax_m1_mnmon.fillna(value=df_tmax_m2_mnmon)
#df_tmax_final.fillna(method='pad', inplace=True)

if event_type == 'snow':
    df_tmax_final.fillna(value=32., inplace=True)
elif event_type == 'liquid':
    df_tmax_final.fillna(value=38., inplace=True)

df_tmax_final.index.name = 'HRU'
#df_tmax_final.head(12)
# print df_tmax_final.head(12).iloc[:,10289]

# %%
# Write the result to a file
df_tmax_final.T.to_csv('r%s_tmax_%s_mean_v3.csv' % (region, event_type))

# %%
snow_hold = df_tmax_final.copy()

# %%
diff = pd.np.subtract(df_tmax_final, snow_hold)
diff.T.to_csv('r%s_tmax_diff_mean_v3.csv' % region)

# %% [markdown]
# ## Update the input parameter file

# %%
# print df_tmax_final.values[0,:].flatten()
# print df_tmax_final.values[0,:]
print df_tmax_final.iloc[:,2280]

# %%
import prms_lib as prms

prms_region_dir = '/media/scratch/PRMS/regions/r%s' % region
prms_region_file = '%s/daymet.control.param' % prms_region_dir
prms_region_outfile = '%s/daymet.params.tmax_allsnow' % prms_region_dir

# Open input parameter file
params = prms.parameters(prms_region_file)

# Replace either tmax_allrain or tmax_allsnow values with the new ones
if event_type == 'snow':
    params.replace_values('tmax_allsnow', df_tmax_final.T.values)
elif event_type == 'liquid':
    params.replace_values('tmax_allrain_offset', df_tmax_final.T.values)

# Write out new input parameter file
#params.write_param_file(prms_region_outfile)
params.write_select_param_file('crap_new_param', ['tmax_allsnow', 'tmax_allrain_offset'])

# %% [markdown]
# ## Unused below here

# %%
# Get min and max for each HRU
themins = df5.min(axis=0)
themaxs = df5.max(axis=0)
themeds = df5.median(axis=0)
themns = df5.mean(axis=0)

mm = pd.concat({'themin': themins, 'themax': themaxs, 'themedian': themeds, 'themean': themns}, axis=1)
mm.head()
mm.to_csv('hru_ranges_rain_v2.csv')

# %%
# Read the daymet vapor pressure
df_vp1 = pd.read_csv(vpfile, na_values=['MEAN(Pa)'], header=True)
df_vp1.drop(0, inplace=True)

df_vp1.rename(columns={df_vp1.columns[0]:'thedate'}, inplace=True)
df_vp1['thedate'] = pd.to_datetime(df_vp1['thedate'])
df_vp1.set_index('thedate', inplace=True)
df_vp1 = df_vp1[st:en]
# Rename the columns to something a little more descriptive
#df2.rename(columns=lambda x: 'HRU_%s' % x, inplace=True)

# Fill in the missing days with NAN
df_vp1 = df_vp1.resample('D', how='mean')
df_vp = df_vp1
# df_vp = pd.np.multiply(df_vp1, df3)
# df_vp[df_vp == 0.0] = np.nan

# %%
# Compute saturation vapor pressure
df_es = df1.copy()

# Convert Fahrenheit to Celsius
df_es = (df_es - 32.) * (5./9.)
# print df_es.head()

df_es = 611.2 * np.exp(df_es * 17.67 / (df_es + 243.5))

# print df_es.head()

# %%
# Compute relative humidity
df_rh = pd.np.divide(df_vp, df_es)

# print df_rh.head()

# %%
df_rh_mon = df_rh.resample('M', how='max')
# print df_rh_mon.head()

df_rh_mnmon = df_rh_mon.groupby(df_rh_mon.index.month).max()
print df_rh_mnmon.iloc[:,100]

# %%
# Use df7 and compute the approximate RH given the mean monthly tmax

df_estRH = df7.copy()
df_estRH = (df_estRH - 32.) * (5./9.)
df_estRH = 9.5 * np.exp(-17.27 * df_estRH / (df_estRH + 238.3)) * (10.5 - df_estRH)
print df_estRH.iloc[:,100]

# %%
print df7.iloc[:,100]

# %%
temps_c = np.linspace(0,10,101)
temps_f = temps_c * 1.8 + 32
rh = []

for tt in temps_c:
    tcalc = 9.5 * np.exp(-17.27 * tt / (tt + 238.3)) * (10.5 - tt)
    rh.append(tcalc)

print rh

# %%
for tt, rr in zip(temps_c,rh):
    print tt, rr

# %%
# Get the index for the rh value closest to given value
theval = 30
ii = min(range(len(rh)), key=lambda i: abs(rh[i]-theval))

# Return the temperature associated with that RH
print temps_c[ii]

# %%
# Copy RH and convert to a percentage
tmp1 = df_rh_mnmon.copy() * 100.

def get_tf(x):
    return temps_f[min(range(len(rh)), key=lambda i: abs(rh[i] - x))]


def get_tc(x):
    return temps_c[min(range(len(rh)), key=lambda i: abs(rh[i] - x))]

tmp2 = tmp1.applymap(get_tf)
#tmp2 = temps_c[min(range(len(rh)), key=lambda i: abs(rh[i] - tmp1))]
print tmp2.head(12)


# %%

# %%
