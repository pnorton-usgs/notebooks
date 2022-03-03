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
import efc_lib as efc
import pandas as pd
import numpy as np
import math as mth
import datetime

# %%
sim_var = 'basin_cfs'
obs_var = 'runoff'

# Setup model run information
templatedir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/PRMS_master'
basedir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t3'
basinid = '06191500'
runid = '2015-03-30_1538'
modelrunid = '01549'

workdir = '%s/%s/runs/%s' % (basedir, basinid, runid)
modeldir = '%s/%s' % (templatedir, basinid)

st = datetime.datetime(1981,10,1)
en = datetime.datetime(2010,9,30)

# %%



# Read in mocom file and use regex to specify variable whitespace between fields
mocom_file = '%s/optim_fixed.log' % workdir
mocom = pd.read_csv(mocom_file, sep=',')
maxgen = max(mocom['gennum'])    # Get the number of the last generation
print "Last generation = %d" % maxgen

# Get list of solutions from the last generation in the optimization log
modelrunids = mocom['soln_num'].loc[mocom['gennum'] == maxgen].tolist()

# Load the statvar file
sv = prms.statvar("%s/%s/default.statvar" % (workdir, modelrunid))
statvar_data = sv.data
statvar_data.head()

#statvar_data.to_csv("%s/%s/statvar.csv" % (workdir, modelrunid))
#statvar_data.reset_index(inplace=True)
#statvar_data['thedate'] = pd.to_pydatetime(statvar_data['thedate'])
#statvar_data.set_index(inplace=True)

#plotvars = statvar_data.columns.tolist()
#plotvars.remove(sim_var)
#plotvars.remove(obs_var)



# %% [markdown]
# <b>Examples using efc_lib</b>

# %%
# Grab a timeseries of the runoff 
sf = statvar_data.loc[:,obs_var]

# The efc class requires a TimeSeries
a = efc.efc(sf)

# Make a copy of the sim_var and obs_var, add the H6L7 results to it
b = statvar_data.loc[:,[sim_var,obs_var]]
b['H6L7'] = a.H6L7
b['efc'] = a.efc
b['highlow'] = a.highlow

st = datetime.datetime(1985,10,1)
en = datetime.datetime(2000,9,30)

c = b[st:en]


# %%
# %matplotlib inline
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as dates
from matplotlib.dates import DayLocator, MonthLocator, YearLocator

def set_ylim(subplot, dataset, padding):
    # Set the y-axis limits
    
    # Dataset is assumed to be a single column pandas dataset
    maxy = dataset.max()
    miny = 0.
    #miny = dataset.min()
    
    # Shut off automatic offsets for the y-axis
    subplot.get_yaxis().get_major_formatter().set_useOffset(False)

    # Set the limits on the y-axis
    yrange = maxy - miny
    ypad = abs(yrange * padding)
    subplot.set_ylim([miny - ypad, maxy + ypad])
    
pdf_filename = 'H6L7_%s.pdf' % (basinid)
outpdf = PdfPages(pdf_filename)
    
    
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(17,11), sharex=True)
ax = axes.flatten()

mkrsize = 9.1

data = b[st:en]
cmap = ['#ff0033','#FF9933','#cc0066','#00cc66','#003366']
labels = ['Large flood', 'Small flood', 'High flow pulse', 'Low flow', 'Extreme low flow']

set_ylim(ax[0], data[obs_var], .02)
ax[0].plot(data.index.to_pydatetime(), data[obs_var], c='grey', lw=.5, alpha=0.5)

for xx in range(0,5):
    sdf = data[data['efc'] == xx+1]
    if sdf.shape[0] == 0:
        continue
    
    ax[0].scatter(sdf.index.to_pydatetime(), sdf[obs_var], c=cmap[xx], s=mkrsize, lw=0, alpha=0.7, label=labels[xx])
    
ax[0].xaxis.set_major_locator(YearLocator(base=1, month=1, day=1))
ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim([st, en])
ax[0].set_title('EFC', fontsize=10)
ax[0].legend(loc='upper left', framealpha=0.5)



# **********
# highlow
cmap = ['#00cc66','#ff9933','#9933ff']
labels = ['Low flow', 'Ascending flow', 'Descending flow']

set_ylim(ax[1], data[obs_var], .02)
ax[1].plot(data.index.to_pydatetime(), data[obs_var], c='grey', lw=.5, alpha=0.5)

for xx in range(0,3):
    sdf = data[data['highlow'] == xx+1]
    if sdf.shape[0] == 0:
        continue
    
    ax[1].scatter(sdf.index.to_pydatetime(), sdf[obs_var], c=cmap[xx], s=mkrsize, lw=0, alpha=0.7, label=labels[xx])
    
ax[1].xaxis.set_major_locator(YearLocator(base=1, month=1, day=1))
#ax[cc].xaxis.set_minor_locator(MonthLocator(bymonth=[1,4,7,10], bymonthday=1))
#ax[cc].get_xaxis().set_minor_formatter(dates.DateFormatter('%b'))
ax[1].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[1].set_xlim([st, en])
ax[1].set_title('HighLow', fontsize=10)
ax[1].legend(loc='upper left', framealpha=0.5)


# ***********
# HIGH6.LOW7
data = b[st:en]
cmap = ['#cc9966','#00cc66']
labels = ['High flow', 'Low flow']

set_ylim(ax[2], data[obs_var], .02)
ax[2].plot(data.index.to_pydatetime(), data[obs_var], c='grey', lw=.5, alpha=0.5)

for xx in range(0,2):
    sdf = data[data['H6L7'] == xx+6]
    if sdf.shape[0] == 0:
        continue
    
    ax[2].scatter(sdf.index.to_pydatetime(), sdf[obs_var], c=cmap[xx], s=mkrsize, lw=0, alpha=0.7, label=labels[xx])
    
ax[2].xaxis.set_major_locator(YearLocator(base=1, month=1, day=1))
#ax[cc].xaxis.set_minor_locator(MonthLocator(bymonth=[1,4,7,10], bymonthday=1))
#ax[cc].get_xaxis().set_minor_formatter(dates.DateFormatter('%b'))
ax[2].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[2].set_xlim([st, en])
ax[2].set_title('HIGH6/LOW7', fontsize=10)
ax[2].legend(loc='upper left', framealpha=0.5)

outpdf.savefig()
outpdf.close()

#plt.show()

# %%
def dparse_HL(yr, mo, dy):
    # Date parser for working with the date format from PRMS files

    # Convert to integer first
    yr, mo, dy = [int(x) for x in [yr, mo, dy]]

    dt = datetime.datetime(yr, mo, dy)
    return dt



st = datetime.datetime(1988,10,1)
en = datetime.datetime(1990,9,30)

hl_file = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/PRMS_master/06191500/HIGH6.LOW7_subdivide.06191500'

# Load the high/low data
thecols = ['year', 'month', 'day', 'H6L7', 'obsQ', 'efc']
hl_data = pd.read_csv(hl_file, sep=r"\s+", header=None, names=thecols,
                      parse_dates={'thedate': ['year', 'month', 'day']},
                      date_parser=dparse_HL, index_col='thedate')
hl_ss = hl_data[st:en]

# Grab a timeseries of the runoff from statvar file 
sf = statvar_data.loc[:,obs_var]
sf_ss = sf[st:en]

# Create a combined dataset from sim, obs, and hl
combined = pd.concat([hl_ss.loc[:,['efc', 'H6L7']], sf_ss], axis=1)
print combined.head()

# %%
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(17,11), sharex=True)
ax = axes.flatten()

data = combined

# ***********
# efc
cmap = ['#ff0033','#FF9933','#cc0066','#00cc66','#003366']
labels = ['Large flood', 'Small flood', 'High flow pulse', 'Low flow', 'Extreme low flow']

set_ylim(ax[0], data[obs_var], .02)
ax[0].plot(data.index.to_pydatetime(), data[obs_var], c='grey', lw=.5, alpha=0.5)

for xx in range(0,5):
    sdf = data[data['efc'] == xx+1]
    if sdf.shape[0] == 0:
        continue
    
    ax[0].scatter(sdf.index.to_pydatetime(), sdf[obs_var], c=cmap[xx], marker='s', s=9.1, lw=0, alpha=0.7, label=labels[xx])
    
ax[0].xaxis.set_major_locator(YearLocator(base=1, month=1, day=1))
ax[0].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[0].set_xlim([st, en])
ax[0].set_title('EFC', fontsize=10)
ax[0].legend(loc='upper left', framealpha=0.5)


# ***********
# HIGH6.LOW7
data = b[st:en]
cmap = ['#cc9966','#00cc66']
labels = ['High flow', 'Low flow']

set_ylim(ax[1], data[obs_var], .02)
ax[1].plot(data.index.to_pydatetime(), data[obs_var], c='grey', lw=.5, alpha=0.5)

for xx in range(0,2):
    sdf = data[data['H6L7'] == xx+6]
    if sdf.shape[0] == 0:
        continue
    
    ax[1].scatter(sdf.index.to_pydatetime(), sdf[obs_var], c=cmap[xx], marker='s', s=9.1, lw=0, alpha=0.7, label=labels[xx])
    
ax[1].xaxis.set_major_locator(YearLocator(base=1, month=1, day=1))
#ax[cc].xaxis.set_minor_locator(MonthLocator(bymonth=[1,4,7,10], bymonthday=1))
#ax[cc].get_xaxis().set_minor_formatter(dates.DateFormatter('%b'))
ax[1].get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax[1].set_xlim([st, en])
ax[1].set_title('HIGH6/LOW7', fontsize=10)
ax[1].legend(loc='upper left', framealpha=0.5)


plt.show()

# %%
pd.set_option('mode.chained_assignment', None)

pdf_filename = 'H6L7_sim_v_obs_%s.pdf' % (basinid)
outpdf = PdfPages(pdf_filename)
    
    
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(17,11), sharex=True)
#ax = axes.flatten()

mkrsize = 9.1

st = datetime.datetime(1984,10,1)
en = datetime.datetime(1990,9,30)

# ***********
# HIGH6.LOW7 
data = b[st:en]


cmap = ['#cc9966','#00cc66']
labels = ['High flow', 'Low flow']

set_ylim(ax, data[obs_var], .05)
ax.plot(data.index.to_pydatetime(), data[obs_var], c='grey', lw=.5, alpha=0.5)


# This sets simulated streamflow to nan whenever it is not a high flow (which is 6)
sdf = data.copy()
sdf.loc[:, sim_var][sdf['H6L7'] == 7] = np.nan
ax.plot(sdf.index.to_pydatetime(), sdf[sim_var], c='red', lw=.5, alpha=1.0)

# This sets simulated streamflow to nan whenever it is not a low flow (which is 7)
sdf = data.copy()
sdf.loc[:, sim_var][sdf['H6L7'] == 6] = np.nan
ax.plot(sdf.index.to_pydatetime(), sdf[sim_var], c='blue', lw=.5, alpha=1.0)

for xx in range(0,2):
    sdf = data[data['H6L7'] == xx+6]
    if sdf.shape[0] == 0:
        continue
    
    ax.scatter(sdf.index.to_pydatetime(), sdf[obs_var], c=cmap[xx], marker='s', s=mkrsize, lw=0, alpha=0.7, label=labels[xx])

ax.xaxis.set_major_locator(YearLocator(base=1, month=1, day=1))
ax.xaxis.set_minor_locator(MonthLocator(bymonth=[1,4,7,10], bymonthday=1))
ax.get_xaxis().set_minor_formatter(dates.DateFormatter('%b'))
ax.get_xaxis().set_major_formatter(dates.DateFormatter('\n%Y'))
ax.set_xlim([st, en])
ax.set_title('HIGH6.LOW7 sim vs. obs', fontsize=10)
ax.legend(loc='upper left', framealpha=0.5)

outpdf.savefig()
outpdf.close()

#plt.show()

# %%
