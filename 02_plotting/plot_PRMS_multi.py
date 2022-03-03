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
import pandas as pd
import numpy as np
import math as mth
import datetime
import lmoments as lmom

import bokeh
from bokeh.plotting import output_notebook
#from bokeh.models import LinearAxis, PrintfTickFormatter
#from bokeh.plotting import *
from bokeh import plotting
output_notebook()
#from bokeh.plotting import line, show

print "bokeh version:", bokeh.__version__


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
from bokeh import plotting
colors = [
      "#1f77b4",
      "#ff7f0e", 
      "#2ca02c", "#98df8a",
      "#d62728", "#ff9896",
      "#9467bd", "#c5b0d5",
      "#8c564b", "#c49c94",
      "#e377c2", "#f7b6d2",
      "#7f7f7f", "#ffbb78",
      "#bcbd22", "#dbdb8d",
      "#17becf", "#9edae5"
    ]
def make_color_cycle():
    ic = iter(colors)
    while True:
        try:
            yield ic.next()
        except StopIteration as e:
            ic = iter(colors)


def multiplot(df, colsets, **kwargs):
    x_range = [None]

    x=df.index
    columns = dict([[c, df[c]] for c in df.columns])
    columns['x'] = df.index
    datasource = plotting.ColumnDataSource(columns)
    def default_plot(c__, title=False):
        defaults = dict(f=plotting.line)
        defaults.update(kwargs)
        if type(c__) == type({}):
            defaults.update(c__)
            colname = defaults.pop('col')
        else:
            colname = c__
        plotf = defaults.pop('f')
        l = plotf('x', colname, source=datasource,
                x_range=x_range[0], color=color_cycle.next(), 
                legend=colname, title=title or colname, **defaults)
        x_range[0] = l.x_range
        return l
    rows = []

    color_cycle = make_color_cycle()
    for col_or_colset in colsets:
        if type(col_or_colset) == type([]):
            plotting.hold(True)
            l = default_plot(col_or_colset[0])
            for c in col_or_colset[1:]:
                default_plot(c)
            plotting.hold(False)
        else:
            l = default_plot(col_or_colset)

        rows.append([l])
        plotting.figure()
    plotting.gridplot(rows)
    plotting.show()


# %%

pdf_filename = '%s/%s/pdf/%s_%s_%s_statvar.pdf' % (basedir, basinid, basinid, runid, modelrunid)

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
statvar_data.drop(['rec'], axis=1, inplace=True)
#statvar_data.reset_index(inplace=True)
#statvar_data['thedate'] = pd.to_pydatetime(statvar_data['thedate'])
#statvar_data.set_index(inplace=True)

plotvars = statvar_data.columns.tolist()
plotvars.remove(sim_var)
plotvars.remove(obs_var)

# Plot the observed streamflow in blue
#multiplot(statvar_data, [obs_var, sim_var], plot_width=900, plot_height=200, tools="pan,box_zoom,wheel_zoom,select,resize,crosshair")

plot_config = dict(plot_width=900, plot_height=200, h_symmetry=False, 
                   min_border=2, min_border_left=60, tools="pan,box_zoom,wheel_zoom,box_select,tap,crosshair,reset")

columns = dict(([[c, statvar_data[c]] for c in statvar_data.columns]))
columns['x'] = statvar_data.index
datasource = plotting.ColumnDataSource(columns)

rows = []

sym_size=6.
xr = None

cfig = plotting.figure(x_range=xr, x_axis_type="datetime", title='', **plot_config)
cfig.below[0].formatter.formats = dict(years=['%Y'],
                                       months=['%b %Y'],
                                       days=['%d %b %Y'])
for ii, pv in enumerate([sim_var, obs_var]):
    colors = ['red', 'green']
    cfig.scatter('x', pv, source=datasource, legend=pv, color=colors[ii], size=sym_size, alpha=.6, nonselection_alpha = 0.1)
    
    if xr is None:
        xr = cfig.x_range
        
rows.append([cfig])

for pv in plotvars:
    cfig = plotting.figure(x_range=xr, x_axis_type="datetime", title='', **plot_config)
    cfig.below[0].formatter.formats = dict(years=['%Y'],
                                           months=['%b %Y'],
                                           days=['%d %b %Y'])
    
    #yformatter = PrintfTickFormatter(format="%' 0.2f")
    #yaxis = LinearAxis(formatter=yformatter)
    #cfig.yaxis()[0].formatter = LinearAxis(formatter=yformatter)
        
    #yaxis()[0].formatter=PrintfTickFormatter()
    #yaxis=LinearAxis()
    #cfig.left[0].formatter.formats = yaxis
    #cfig.add_layout(LinearAxis(formatter=PrintfTickFormatter("%' 0.2f")), 'left')
    
    cfig.scatter('x', pv, source=datasource, legend=pv, color='green', size=sym_size, nonselection_alpha = 0.1, 
                 nonselection_fill_color = "pink")
    

    
    
    if xr is None:
        xr = cfig.x_range
        
    rows.append([cfig])
    #plotting.figure(x_range=xr)
    #plotting.figure(title=pv)
    
grd = plotting.gridplot(rows)
plotting.show(grd)



# %%
from bokeh.plotting import *
from bokeh.models import BoxSelectTool

N = 100

x = np.linspace(0, 4*np.pi, N)
y = np.sin(x)

#output_file("scatter_selection.html", title="scatter_selection.py example")

TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select"

p1 = figure(title="selection on mouseup", tools=TOOLS)

p1.circle(x, y, color="red", size=6)
select_tool = p1.select(dict(type=BoxSelectTool))
select_tool.select_every_mousemove = False

p2 = figure(title="selection on mousemove", tools=TOOLS)
p2.circle(x, y, marker="square", color="green", size=6)

select_tool = p2.select(dict(type=BoxSelectTool))
select_tool.select_every_mousemove = True

p3 = figure(title="default highlight", tools=TOOLS)
p3.circle(x,y, color="#FF00FF", size=6)

p4 = figure(title="custom highlight", tools=TOOLS)
p4.square(x,y, color="blue", size=6,
    nonselection_fill_color="#FFFF00", nonselection_fill_alpha=1)

show(VBox(p1,p2,p3,p4))  # open a browser

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
