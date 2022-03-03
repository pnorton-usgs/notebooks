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

import prms_lib as prms
import pandas as pd
import numpy as np
import math as mth
import datetime
import lmoments as lmom

# %%
import jinja2
import json
import numpy as np
import matplotlib.pyplot as plt

import mpld3
from mpld3 import plugins, utils


class HighlightLines(plugins.PluginBase):
    """A plugin to highlight lines on hover"""

    JAVASCRIPT = """
    mpld3.register_plugin("linehighlight", LineHighlightPlugin);
    LineHighlightPlugin.prototype = Object.create(mpld3.Plugin.prototype);
    LineHighlightPlugin.prototype.constructor = LineHighlightPlugin;
    LineHighlightPlugin.prototype.requiredProps = ["line_ids"];
    LineHighlightPlugin.prototype.defaultProps = {alpha_bg:0.3, alpha_fg:1.0}
    function LineHighlightPlugin(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };

    LineHighlightPlugin.prototype.draw = function(){
      for(var i=0; i<this.props.line_ids.length; i++){
         var obj = mpld3.get_element(this.props.line_ids[i], this.fig),
             alpha_fg = this.props.alpha_fg;
             alpha_bg = this.props.alpha_bg;
         obj.elements()
             .on("mouseover", function(d, i){
                            d3.select(this).transition().duration(50)
                              .style("stroke-opacity", alpha_fg); })
             .on("mouseout", function(d, i){
                            d3.select(this).transition().duration(200)
                              .style("stroke-opacity", alpha_bg); });
      }
    };
    """

    def __init__(self, lines):
        self.lines = lines
        self.dict_ = {"type": "linehighlight",
                      "line_ids": [utils.get_id(line) for line in lines],
                      "alpha_bg": lines[0].get_alpha(),
                      "alpha_fg": 1.0}


#N_paths = 50
#N_steps = 100

#x = np.linspace(0, 10, 100)
#y = 0.1 * (np.random.random((N_paths, N_steps)) - 0.5)
#y = y.cumsum(1)

#fig, ax = plt.subplots(subplot_kw={'xticks': [], 'yticks': []})
#lines = ax.plot(x, y.T, color='blue', lw=4, alpha=0.1)
#plugins.connect(fig, HighlightLines(lines))
#mpld3.display()
#mpld3.show()

# %%
def mean_monthly(data):
    # Compute the mean monthly values for observations
    # NOTE: This assumes the data passed in is a pandas time series of
    #       monthly mean values

    # Copy the data so we don't inadvertently end up with a reference to the data
    mn_monthly = pd.DataFrame(data.copy())
    mn_monthly = mn_monthly.groupby(mn_monthly.index.month).mean().iloc[:,0]

    return mn_monthly


# %% [markdown]
# ### Plot observed versus simulated monthly and daily streamflow

# %%

# %%
model_id = '2015-03-06_1003'

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/tmp_data'
mocom_file = '%s/optim_fixed.log' % workdir

# Read in mocom file and use regex to specify variable white space between fields
mocom = pd.read_csv(mocom_file, sep=',')
maxgen = max(mocom['gennum'])    # Get the number of the last generation
print "Last generation = %d" % maxgen

# Get list of solutions from the last generation in the optimization log
runids = mocom['soln_num'].loc[mocom['gennum'] == maxgen].tolist()


workdir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t3/06191500/runs/%s' % model_id

#luca_file = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/tmp_data/test1_luca.statvar'
st = datetime.datetime(1985,10,1)
en = datetime.datetime(1990,9,30)

first = True

#fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(11,17), subplot_kw={'xticks': [], 'yticks': []})
fig, axes = plt.subplots(nrows=1, ncols=1, subplot_kw={'xticks': [], 'yticks': []})
#ax = axes.flatten()


# Plot each of the MOCOM solutions against observed and the LUCA final solution
for rr in runids:
    if first:
        # Only need to read the observations once
        obs_streamflow_file = '%s/%05d/streamflow.data' % (workdir, rr)
        
        obs_streamflow = prms.streamflow(obs_streamflow_file).data
        obs_subset = obs_streamflow[st:en].iloc[:,0]
        
        obs_mon_mn = obs_subset.resample('M', how='mean')
        obs_ann = obs_subset.resample('A-SEP', how='mean')
        
        # Plot the observed streamflow in blue
        monplots = axes.plot(obs_mon_mn.index.to_pydatetime(), obs_mon_mn, color='grey', label='runoff')
        
        first = False
        
    sim_streamflow_file = '%s/%05d/default.statvar' % (workdir, rr)
    
    sim_streamflow = prms.statvar(sim_streamflow_file).data
    sim_streamflow.drop(['rec'], axis=1, inplace=True)

    sim_subset = sim_streamflow[st:en].iloc[:,0]

    # Create monthly and annual subsets
    sim_mon_mn = sim_subset.resample('M', how='mean')
    sim_ann = sim_subset.resample('A-SEP', how='mean')
    
    # Plot the MOCOM pareto set simulated streamflow in red
    monplots = axes.plot(sim_mon_mn.index.to_pydatetime(), sim_mon_mn, color='red', alpha=0.1, label='MOCOM')

#ax[0].set_title('Monthly Streamflow', fontsize=12)
#plt.suptitle('Model_ID = %s' % model_id, fontsize=15)

# Replot observed
#ax[0].plot(obs_mon_mn.index.to_pydatetime(), obs_mon_mn, color='grey', label='runoff')

#ax[0].xaxis.set_major_locator(YearLocator(base=1, month=1, day=1))
#ax[0].xaxis.set_minor_locator(MonthLocator(bymonth=[1,4,7,10], bymonthday=1))
#ax[0].get_xaxis().set_minor_formatter(dates.DateFormatter('%b'))
#ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
#ax[0].set_xlim([st, en])

# Clean up the list of labels
#handles, labels = ax[0].get_legend_handles_labels()
#newLabels, newHandles = [], []

#for handle, label in zip(handles, labels):
#  if label not in newLabels:
#    newLabels.append(label)
#    newHandles.append(handle)
#ax[0].legend(newHandles, newLabels, loc='upper right', framealpha=0.5)
# ------------------------------------------------------------------------------------


plugins.connect(fig, HighlightLines(monplots))
mpld3.display()



# Plot the daily values for the observed, LUCA final, and last MOCOM run
#ax[1].plot(obs_subset.index.to_pydatetime(), obs_subset, color='grey', label='runoff')
 #ax[1].plot(luca_sub.index.to_pydatetime(), luca_sub, color='green', label='LUCA')
#ax[1].plot(sim_subset.index.to_pydatetime(), sim_subset, color='red', label='MOCOM')
#ax[1].xaxis.set_major_locator(YearLocator(base=1, month=1, day=1))
#ax[1].xaxis.set_minor_locator(MonthLocator(bymonth=[1,4,7,10], bymonthday=1))
#ax[1].get_xaxis().set_minor_formatter(dates.DateFormatter('%b'))
#ax[1].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
#ax[1].set_xlim([st, en])
#ax[1].legend(loc='upper right', framealpha=0.5)

#ax[1].set_title('Daily Streamflow', fontsize=12)

#outpdf.savefig()
#outpdf.close()


# %%
# Plot the annual values for the observed, LUCA final, and last MOCOM run
#obs_ann.plot(color='blue', figsize=(15,5))
#luca_ann.plot(color='green', figsize=(15,5))
#sim_ann.plot(color='red', figsize=(15,5))
fig, ax = plt.subplots(1)
ax.plot(obs_ann.index.to_pydatetime(), obs_ann, color='grey', label='runoff')
ax.plot(luca_ann.index.to_pydatetime(), luca_ann, color='green', label='LUCA')
ax.plot(sim_ann.index.to_pydatetime(), sim_ann, color='red', label='MOCOM')
ax.xaxis.set_major_locator(YearLocator(base=1, month=1, day=1))
#plt.xaxis.set_minor_locator(MonthLocator(bymonth=[1,4,7,10], bymonthday=1))
#plt.get_xaxis().set_minor_formatter(dates.DateFormatter('%b'))
ax.get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax.set_xlim([st, en])
ax.legend(loc='upper right', framealpha=0.5)

ax.set_title('Annual Streamflow', fontsize=12)

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
