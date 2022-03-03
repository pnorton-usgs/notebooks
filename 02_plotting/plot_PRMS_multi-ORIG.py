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
# Rough process
#     specify the runs/MOCOMid directory
#     specify the runid
#     modify the control file to include additional statvars
#        basin_tmin, basin_tmax, basin_ppt, basin_cfs, basin_potet, basin_actet,
#        basin_potsw, basin_gwflow_cfs, basin_sroff_cfs, basin_ssflow_cfs,
#        basin_storage, SR
#     run PRMS
#     process the statvar file
#     create daily, annual, and monthly values as needed


# %%
# Uncomment the following line to write to pdf without opening a window
#matplotlib.use('Agg')

# %matplotlib inline
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as dates
from matplotlib.dates import DayLocator, MonthLocator, YearLocator
from matplotlib.ticker import ScalarFormatter
import mpld3
from mpld3 import plugins, utils

import prms_lib as prms
import prms_objfcn as objfcn
import pandas as pd
import numpy as np
import math as mth
import datetime
#import lmoments as lmom

# %%
def mean_monthly(data):
    # Compute the mean monthly values for observations
    # NOTE: This assumes the data passed in is a pandas time series of
    #       monthly mean values

    # Copy the data so we don't inadvertently end up with a reference to the data
    mn_monthly = pd.DataFrame(data.copy())
    mn_monthly = mn_monthly.groupby(mn_monthly.index.month).mean().iloc[:,0]

    return mn_monthly

def flow_duration(ds):
    """Compute the flow duration for a given dataset"""
    # NOTE: This exists in prms_objfcn.py
    
    # See http://pubs.usgs.gov/sir/2008/5126/section3.html 
    # for the approach used to compute the flow duration
    
    # We only want valid values, sort the values in descending order
    rankedQ = sorted(ds[ds.notnull()], reverse=True)
    
    # Compute the exceedence probability for each Q
    prob = np.arange(len(rankedQ), dtype=np.float_) + 1.0
    prob = 100 * (prob / (len(rankedQ) + 1.0))

    # Return a dataframe of the flow duration / exceedence probability
    return pd.DataFrame({'exceedence': prob, 'Q': rankedQ}, columns=['exceedence', 'Q'])

    

def set_ylim(subplot, dataset, padding):
    # Set the y-axis limits
    
    # Dataset is assumed to be a single column pandas dataset
    maxy = dataset.max()
    miny = dataset.min()
    
    # Shut off automatic offsets for the y-axis
    subplot.get_yaxis().get_major_formatter().set_useOffset(False)

    # Set the limits on the y-axis
    yrange = maxy - miny
    ypad = abs(yrange * padding)
    subplot.set_ylim([miny - ypad, maxy + ypad])
    

def get_complete_wyears(ds):
    # Steps to remove "bad" years from dataset
    # "bad" is defined as any year where one or more days are missing

    ds['wyear'] = ds.index.year
    ds['month'] = ds.index.month
    ds['wyear'] = np.where(ds['month']>9, ds['wyear']+1, ds['wyear'])

    b = ds[ds['runoff'].isnull()]['wyear'].tolist()

    # Create set of unique 'bad' years
    badyears = {x for x in b}

    c = ds.loc[~ds['wyear'].isin(badyears)]
    c.drop(['wyear', 'month'], axis=1, inplace=True)
    return c


# %% [markdown]
# ### Plot observed versus simulated monthly and daily streamflow

# %%
sim_var = 'basin_cfs'
obs_var = 'runoff'

# Setup model run information
basedir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t4'
basinid = '4.39_03410500'
runid = '2015-04-29_1424'
#modelrunid = '02651'

workdir = '%s/%s/runs/%s' % (basedir, basinid, runid)
opt_file = '%s/optim_fixed.log' % workdir
opt_out = '%s/optim_fixed.log' % workdir

statvar_file = 'daymet.statvar'

nwis_file = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/PRMS_master/nwis_sites.tab'

# TODO: Read start and end times from the config file
st = datetime.datetime(1982,10,1)
en = datetime.datetime(2010,9,30)

# %%
nwis = pd.read_csv(nwis_file, sep='\t', dtype=np.str_)

#print nwis['station_nm'].loc[nwis['site_no'] == '01010070'].tolist()[0]

# %%
reload(prms)
reload(objfcn)

# Read in mocom file and use regex to specify variable whitespace between fields
#mocom_file = '%s/optim_fixed.log' % workdir
mocom = pd.read_csv(opt_file, sep=',')
maxgen = max(mocom['gennum'])    # Get the number of the last generation
print "Last generation = %d" % maxgen

# Get list of solutions from the last generation in the optimization log
modelrunids = mocom['soln_num'].loc[mocom['gennum'] == maxgen].tolist()
#modelrunid = '%05d' % modelrunids[-1]

# Get the 'best' pareto member based on Nash-Sutcliffe
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

modelrunid = '%05d' % df.sort(['NS','set']).iloc[-1].set
#------------------------------




site_no = basinid.split('_')[1]
#site_no = basinid
site_name = nwis['station_nm'].loc[nwis['site_no'] == site_no].tolist()[0]


pdf_filename = '%s/%s/pdf/%s_%s_%s_statvar.pdf' % (basedir, basinid, basinid, runid, modelrunid)

# Load the statvar file
sv = prms.statvar("%s/%s/%s" % (workdir, modelrunid, statvar_file))
statvar_data = sv.data

first_wyr = statvar_data.index.min().year+1
last_wyr = statvar_data.index.max().year
print 'Valid water year range: %d to %d' % (first_wyr, last_wyr)

plotvars = statvar_data.columns.tolist()
plotvars.remove(sim_var)
plotvars.remove(obs_var)
plotvars.remove('basin_ppt')
plotvars.remove('basin_tmax')
plotvars.remove('basin_tmin')
print 'Highest daily value on', statvar_data.idxmax(axis=0)['runoff']

title_size = 13    # Font size for the plot titles

# Get years of daily data without missing values
good_years = get_complete_wyears(statvar_data)

# Nash-Sutcliffe
print 'Nash-Sutcliffe:', -1 * objfcn.compute_objfcn('NS', 'daily', statvar_data, 'runoff', 'basin_cfs')


# Get the year with greatest annual flow
#ann_tmp = statvar_data.resample('A-SEP', how='mean')
ann_tmp = good_years.resample('A-SEP', how='mean')

plot_yr = {'High': ann_tmp.idxmax(axis=0)['runoff'].year,
           'Low': ann_tmp.idxmin(axis=0)['runoff'].year,
           'Median': ann_tmp.ix[(ann_tmp.runoff-ann_tmp['runoff'].median()).abs().argsort()[:1]].index.year}

outpdf = None    # outpdf is instantiated in the loop below

for kk, yy in plot_yr.iteritems():
    # Reset date range to maximum year
    st = datetime.datetime(yy-1,10,1)
    en = datetime.datetime(yy,9,30)

    plot_data = statvar_data[st:en]
    
    maxrows = 5
    maxcols = 1

    cc = maxrows * maxcols    # used to track plots on a page
    newYear = True
    
    for ii,vv in enumerate(plotvars):
        if cc == maxrows * maxcols:
            if outpdf is None:
                outpdf = PdfPages(pdf_filename)
                newYear = False
            elif newYear:
                newYear = False
            else:
                plt.suptitle('basin: %s (%s)\nrunid: %s (%s)\nYear: %d (%s flow)' % (basinid, site_name, runid, modelrunid, yy, kk), fontsize=title_size)
                outpdf.savefig()
                
            cc = 0
            
            fig, axes = plt.subplots(nrows=maxrows, ncols=maxcols, figsize=(17,11), sharex=True)
            ax = axes.flatten()
            
            # Plot obs/sim plot first on new page
            stuff = ax[cc].plot(plot_data.index.to_pydatetime(), plot_data[obs_var], color='grey', label=obs_var)
            stuff = ax[cc].plot(plot_data.index.to_pydatetime(), plot_data[sim_var], color='red', label=sim_var)
            
            ax[cc].xaxis.set_major_locator(YearLocator(base=1, month=1, day=1))
            #ax[cc].xaxis.set_minor_locator(MonthLocator(bymonth=[1,4,7,10], bymonthday=1))
            #ax[cc].get_xaxis().set_minor_formatter(dates.DateFormatter('%b'))
            ax[cc].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
            ax[cc].set_xlim([st, en])
            ax[cc].set_title('Simulated v. observed daily streamflow', fontsize=10)
            ax[cc].legend(loc='upper right', framealpha=0.5)

            # Create a secondary y-axis
            ax2 = ax[cc].twinx()
            ax2.set_xlim([st, en])

            # Set the secondary y-axis limit
            set_ylim(ax2, plot_data['basin_ppt'], 0.02)

            # Plot precipitation as a series of vertical lines
            ax2.vlines(plot_data.index.to_pydatetime(), [0], plot_data['basin_ppt'], color='blue', alpha=0.4)
            ax2.invert_yaxis()
            
            cc += 1    # increment cc after each page header plot
            
            
        # Shut off automatic offsets for the y-axis
        #ax[ii+1].get_yaxis().get_major_formatter().set_useOffset(False)

        # Set the limits on the y-axis
        set_ylim(ax[cc], plot_data[vv], 0.02)

        stuff = ax[cc].plot(plot_data.index.to_pydatetime(), plot_data[vv], color='green', label=vv)
        ax[cc].xaxis.set_major_locator(YearLocator(base=1, month=1, day=1))
        ax[cc].xaxis.set_minor_locator(MonthLocator(bymonth=[1,4,7,10], bymonthday=1))
        ax[cc].get_xaxis().set_minor_formatter(dates.DateFormatter('%b'))
        ax[cc].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))

        ax[cc].set_xlim([st, en])
        ax[cc].legend(loc='upper right', framealpha=0.5)
        #ax[ii+1].set_title(vv, fontsize=12)

        # Create a secondary y-axis
        ax2 = ax[cc].twinx()
        ax2.set_xlim([st, en])

        # Set secondary y-axis limit
        set_ylim(ax2, plot_data['basin_ppt'], 0.02)

        # Plot precipitation as a series of vertical lines
        ax2.vlines(plot_data.index.to_pydatetime(), [0], plot_data['basin_ppt'], color='blue', alpha = 0.4)
        ax2.invert_yaxis()
        
        cc += 1

    plt.suptitle('basin: %s (%s)\nrunid: %s (%s)\nYear: %d (%s flow)' % (basinid, site_name, runid, modelrunid, yy, kk), fontsize=title_size)
    #fig.gca().set_xlim([st, en])
    #fig.gca().xaxis.set_major_locator(MonthLocator(bymonth=[1,4,7,10], bymonthday=1))
    #fig.gca().get_xaxis().set_major_formatter(dates.DateFormatter('%Y\n%b'))
    #fig.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    
    # Remove any unused subplots
    while cc < maxrows * maxcols:
        fig.delaxes(ax[cc])
        cc += 1
        
    outpdf.savefig()

# Compute the exceedence curves for obs and sim for the entire period
st = datetime.datetime(first_wyr-1,10,1)
en = datetime.datetime(last_wyr,9,30)
plot_data = statvar_data[st:en]

fd_obs_period = flow_duration(plot_data.runoff)
fd_sim_period = flow_duration(plot_data.basin_cfs)

if min(fd_obs_period.Q) == 0.0:
    fd_obs_period.Q += 0.1
if min(fd_sim_period.Q) == 0.0:
    fd_sim_period.Q += 0.1
    
maxrows = 5
maxcols = 2
fig, axes = plt.subplots(nrows=maxrows, ncols=maxcols, figsize=(17,11), sharex=True) 
ax = axes.flatten()

cc = 0

for yy in range(first_wyr, last_wyr+1):
    # Plot the flow duration curves
    if cc == maxrows * maxcols:
        cc = 0
        plt.suptitle('basin: %s (%s)\nrunid: %s (%s)\nExceedence Curves' % (basinid, site_name, runid, modelrunid), fontsize=title_size)
        outpdf.savefig()
        fig, axes = plt.subplots(nrows=maxrows, ncols=maxcols, figsize=(17,11), sharex=True) 
        ax = axes.flatten()
        
    st = datetime.datetime(yy-1,10,1)
    en = datetime.datetime(yy,9,30)
    
    plot_data = statvar_data[st:en]
    
    fd_obs = flow_duration(plot_data.runoff)
    fd_sim = flow_duration(plot_data.basin_cfs)
    
    if min(fd_obs.Q) == 0.0:
        fd_obs.Q += 0.1
    if min(fd_sim.Q) == 0.0:
        fd_sim.Q += 0.1

    ax[cc].semilogy(fd_obs_period['exceedence'], fd_obs_period['Q'], 'k--', color='black', alpha=0.5)
    ax[cc].semilogy(fd_sim_period['exceedence'], fd_sim_period['Q'], 'k--', color='red', alpha=0.5)
    ax[cc].semilogy(fd_obs['exceedence'], fd_obs['Q'], color='grey', label=obs_var)
    ax[cc].semilogy(fd_sim['exceedence'], fd_sim['Q'], color='red', label=sim_var)
    ax[cc].xaxis.set_major_formatter(ScalarFormatter())
    ax[cc].yaxis.set_major_formatter(ScalarFormatter())
    ax[cc].get_xaxis().get_major_formatter().set_useOffset(False)
    ax[cc].set_xlim([-.7,100.5])
    ax[cc].set_title('WY%d' % yy, fontsize=12)
    ax[cc].legend(loc='upper right', framealpha=0.5)

    cc += 1
    
# Remove any unused subplots
while cc < maxrows * maxcols:
    fig.delaxes(ax[cc])
    cc += 1
    
plt.suptitle('basin: %s (%s)\nrunid: %s (%s)\nExceedence Curves' % (basinid, site_name, runid, modelrunid), fontsize=title_size)
outpdf.savefig()


# Close the pdf file and finish
outpdf.close()


# %%
from matplotlib.ticker import ScalarFormatter
    

s1 = sorted(plot_data['runoff'], reverse=True)

p1 = np.arange(len(s1), dtype=np.float_) + 1.0
p1 = 100 * (p1 / (len(s1) + 1.0))

df = {'exceedence': p1,
      'Q': s1}

newdata = pd.DataFrame(df, columns=['exceedence', 'Q'])
newdata.head()

#p1 = []
#for ii, kk in enumerate(s1):
#    p1.append(100*(float(ii+1)/float(len(s1)+1)))


#s2 = sorted(plot_data['basin_cfs'], reverse=True)

#p2 = []
#for ii, kk in enumerate(s2):
#    p2.append(100*(float(ii+1)/float(len(s2)+1)))

#fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(17,11), sharex=True)
#axes.semilogy(p1, s1, color='grey')
#axes.semilogy(p2, s2, color='red')
#axes.yaxis.set_major_formatter(ScalarFormatter())
#axes.get_xaxis().get_major_formatter().set_useOffset(False)
#axes.set_xlim([-.9,100.2])

# %%
#ann_tmp[ann_tmp['runoff'] == ann_tmp['runoff'].median()]
#print ann_tmp['runoff']
a = ann_tmp['runoff'].median()
print a

# Return the record that has the first closest runoff value to the median
b = ann_tmp.ix[(ann_tmp.runoff-a).abs().argsort()[:1]]
print b

x = []
x.append(['test',1])
x.append(['test2',2])

print x
print [y[1] for y in x]

# %% [markdown]
# ### Magnificent 7 work

# %%
dd = obs_subset.tolist()

# %%
# First four moments of the distribution of streamflow are:
#     mean, coefficient of variation, skewness, and kurtosis
lmom.samlmu(dd, 4)

# %%
dd_mn_mon = mean_monthly(obs_mon_mn)
dd_mn_mon

# %%
obs_subset.head()

# %%
#gg = obs_subset.groupby(pd.TimeGrouper('M'))
#gg = obs_subset.groupby(obs_subset.index.month).transform('mean')

# Compute the seasonally adjusted daily streamflow
# The following one-liner subtracts the longterm mean monthly value
# from each daily value in that month.
gg = obs_subset.groupby(obs_subset.index.month).transform(lambda x: x - x.mean())
gg.plot()

# %%
sdev = gg.std()
print sdev
ltmean = gg.mean()
print ltmean

# %%
print gg.head()

# %%
# Standardize the streamflow values using the mean and stdev
hh = (gg - ltmean) / sdev
print hh.head()

# %%
import scipy.stats as stats
import statsmodels.formula.api as smf

# %%
# Return the lag-1 autocorrelation
hh.autocorr()

# %%
dat_doy = obs_subset.index.dayofyear

tt = 2 * mth.pi * dat_doy
dat_cos = np.cos(tt)
dat_sin = np.sin(tt)
# hh is the data

outdata = {}
outdata['data'] = hh
outdata['dcos'] = dat_cos
outdata['dsin'] = dat_sin

testdf = pd.DataFrame(outdata)
testdf.head()

#print dat_cos
#print dat_sin
mod = smf.gls(formula='data ~ dcos + dsin', data=testdf)
res = mod.fit()
print res.summary()

# %%

# %%
res.params

# %%
