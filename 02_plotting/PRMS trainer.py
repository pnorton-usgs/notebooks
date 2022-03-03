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
from IPython.html import widgets
from IPython.display import display

import collections
import ConfigParser
import datetime
import os
import re
import shutil
import subprocess
import numpy as np
import prms_lib as prms
import prms_cfg
reload(prms)
reload(prms_cfg)

# %%
os.chdir("/Users/pnorton/Projects/National_Hydrology_Model/notebooks")


# %%
def to_prms_datetime(date):
    """Takes a datetime object and converts it to a string of form
       YYYY,MM,DD,HH,mm,ss"""
    return date.strftime('%Y,%m,%d,%H,%M,%S')

def to_datetime(date_str):
    """Takes a date string of the form 'YYYY-MM-DD HH:mm:ss' (and variations thereof)
       and converts it to a datetime"""
    return datetime.datetime(*[int(x) for x in re.split('-| |:', date_str)])

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

def redistribute_mean(old_vals, new_mean):
    # Redistribute mean value to set of multiple initial values
    # see Hay and Umemoto, 2006 (p. 11)

    ZC = 10.    # Constant to avoid zero values
    new_vals = []

    old_mean = sum(old_vals) / float(len(old_vals))

    for vv in old_vals:
        new_vals.append((((new_mean + ZC) * (vv + ZC)) / (old_mean + ZC)) - ZC)

    return new_vals




def update_input_param(in_param_file, param_file, params):
    # For each parameter, redistribute the mean to the original parameters
    # and update the input parameter file
    #print "Using input parameter file:", prms_input_file
    pobj = prms.parameters(in_param_file)

    for kk, vv in params.iteritems():
        #if vv['newval'] != vv['cmean']:
        orig_vals = pobj.get_var(kk)['values']

        if len(orig_vals) > 1:
            # an array
            new_vals = redistribute_mean(orig_vals, vv['newval'])

            #print 'Size of parameter array, %s: %d' % (pp, len(orig_vals))
            pobj.replace_values(kk, new_vals)
        else:
            # scalars
            #print pobj.get_var(pp)
            pobj.replace_values(kk, vv['newval'])

    # Write the new values back to the input parameter file
    pobj.write_param_file(param_file)


# %%
basinid = '4.39_03410500'
runid = '2015-05-05_1251'
modelrunid = '00920'
configfile = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t4/%s/basin.cfg' % basinid

# config = ConfigParser.SafeConfigParser()
# config.read(configfile)
# cfg = prms_cfg.cfg_get_vars(config)


cfg = prms_cfg.cfg(configfile)

# %%
base_calib_dir = cfg.get_value('base_calib_dir')
working_dir = '%s/%s/runs/%s/%s' % (base_calib_dir, basinid, runid, modelrunid)

start_dir = os.getcwd() 

os.chdir('%s/%s' % (base_calib_dir, basinid))

# Read the param_range_file
param_limits_file = cfg.get_value('param_range_file')
infile = open(param_limits_file, 'r')

param_limits = {}
for row in infile:
    flds = row.strip().split(' ')
    param_limits[flds[0]] = {}
    param_limits[flds[0]]['minval'] = float(flds[2])
    param_limits[flds[0]]['maxval'] = float(flds[1])
infile.close()

#for kk,vv in param_limits.iteritems():
#    print '% 15s: %0.5f %0.5f' % (kk, vv['minval'], vv['maxval'])

# Change into the modelrunid directory
os.chdir(working_dir)

# Open the control file and get the name of the statvar file
ctl = prms.control(cfg.get_value('prms_control_file'))
statvar_file = os.path.basename(ctl.get_var('stat_var_file')['values'][0])
print statvar_file

# Copy the original statvar file to orig_*
cmd_opts = " %s orig_%s" % (statvar_file, statvar_file)
subprocess.call(cfg.get_value('cmd_cp') + cmd_opts, shell=True)
print cfg.get_value('cmd_cp') + cmd_opts

# Make copy of existing *.param file
param_file = cfg.get_value('prms_input_file')
wk_param_file = 'wk_%s' % param_file

# Make copy of starting statvar file
cmd_opts = " %s %s" % (param_file, wk_param_file)
subprocess.call(cfg.get_value('cmd_cp') + cmd_opts, shell=True)



# Compute current mean value for parameters from .param file
pobj = prms.parameters('wk_%s' % param_file)

#print "\nCurrent mean values"
#print '-'*40
for kk, vv in param_limits.iteritems():
    curr_vals = pobj.get_var(kk)['values']
    curr_mn = sum(curr_vals) / float(len(curr_vals))
    param_limits[kk]['cmean'] = curr_mn
    param_limits[kk]['newval'] = curr_mn

#for kk,vv in param_limits.iteritems():
#    print '% 15s: %0.5f %0.5f %0.5f' % (kk, vv['minval'], vv['maxval'], vv['cmean'])

os.chdir(start_dir)

# %%
# %matplotlib inline
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as dates
from matplotlib.dates import DayLocator, MonthLocator, YearLocator
from matplotlib.ticker import ScalarFormatter

from IPython.html import widgets
from IPython.display import display

#def testplot(a=-2, b=2.):
#    plt.plot([a,b],[b,a])

def select_param(sel):
    #stuff.options = param_limits[sel]
    print param_limits[sel]
    
    if limits.max < param_limits[sel]['minval']:    
        limits.max = param_limits[sel]['maxval']
        limits.min = param_limits[sel]['minval']
    else:
        limits.min = param_limits[sel]['minval']
        limits.max = param_limits[sel]['maxval']
        
    
    limits.value = param_limits[sel]['cmean']
    #print '***',sel
    #print '---',param_widget.value
    
def showresult(res):
    print res

def update(a):
    print param_limits[param_widget], a
    param_limits[param_widget]['newval'] = limits
    print param_limits[param_widget]
    
param_widget = widgets.Select(options=param_limits.keys())
init = param_widget.value

limits = widgets.FloatSlider(min=param_limits[init]['minval'], 
                             max=param_limits[init]['maxval'], 
                             value=param_limits[init]['cmean'])
b = widgets.Button(description = 'Set')
b.on_click(update)

j = widgets.interactive(showresult, res=limits)
i = widgets.interactive(select_param, sel=param_widget)
#b = widgets.interactive(update, param=param_widget, val=limits)
display(i)
display(j)
display(b)

print limits.keys
# Next add button to apply changes
# run PRMS
# plot results



#a_slider = widgets.IntSlider(description='First', orientation='horizontal', min=-5, max=5, step=1, value=0)
#b_slider = widgets.FloatSlider(description='Bob', orientation='horizontal', min=-5, max=5, step=0.3, value=0)
#widgets.interact_manual(testplot,a=a_slider,b=b_slider)


# %%
def on_value_change(name, value):
    print name, value
    
def createSliders(widgetlist):
    # Use this is you want multiple widgets per 'line'
    #containerList = [widgets.Box(children=[widgets.FloatSlider(description=kk, min=vv['minval'],
    #                                                           max=vv['maxval'], value=vv['cmean'], step=0.001)])
    #                for kk,vv in widgetlist.iteritems()]
    containerList = [widgets.FloatSlider(description=kk, min=vv['minval'],
                                                         max=vv['maxval'], value=vv['cmean'], step=0.001)
                    for kk,vv in widgetlist.iteritems()]
    
    #for ii in containerList:
    #    ii.on_trait_change(on_value_change, 'value')
    
    
    bigContainer = widgets.Box(children=containerList)
    
    display(bigContainer)
    return bigContainer

def update_params(a):
    for c in bigContainer.children:
        if c.value != param_limits[c.description]['newval']:
            print 'Updated %s' % c.description
            
            param_limits[c.description]['newval'] = c.value
            
def reset_params(a):
    for c in bigContainer.children:
        param_limits[c.description]['newval'] = param_limits[c.description]['cmean']
        c.value = param_limits[c.description]['cmean']
        
bigContainer = createSliders(param_limits)

b = widgets.Button(description = 'Update')
b.on_click(update_params)
display(b)

rst = widgets.Button(description = 'Reset')
rst.on_click(reset_params)
display(rst)

        

# %%
os.chdir(working_dir)
#print param_limits
update_input_param(param_file, wk_param_file, param_limits)

# Convert the config file start/end date strings to datetime objects
st_date_model = to_datetime(cfg.get_value('start_date_model'))
st_date_calib = to_datetime(cfg.get_value('start_date'))
en_date = to_datetime(cfg.get_value('end_date'))

# Run PRMS
cmd_opts = " -C%s -set param_file %s -set start_time %s -set end_time %s" % (cfg.get_value('prms_control_file'), wk_param_file, to_prms_datetime(st_date_model), to_prms_datetime(en_date))
print cfg.get_value('cmd_prms') + cmd_opts


subprocess.call(cfg.get_value('cmd_prms') + cmd_opts, shell=True)


# %%
def bigplot():
    #statvar_file = 'daymet.statvar'
    sv_orig = prms.statvar("%s/orig_%s" % (working_dir, statvar_file))
    statvar_orig_data = sv_orig.data
    
    sv = prms.statvar("%s/%s" % (working_dir, statvar_file))
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
    
    good_years = get_complete_wyears(statvar_data)

    # Nash-Sutcliffe
    #print 'Nash-Sutcliffe:', -1 * objfcn.compute_objfcn('NS', 'daily', statvar_data, 'runoff', 'basin_cfs')


    # Get the year with greatest annual flow
    #ann_tmp = statvar_data.resample('A-SEP', how='mean')
    ann_tmp = good_years.resample('A-SEP', how='mean')

    #plot_yr = {'High': ann_tmp.idxmax(axis=0)['runoff'].year,
    #           'Low': ann_tmp.idxmin(axis=0)['runoff'].year,
    #           'Median': ann_tmp.ix[(ann_tmp.runoff-ann_tmp['runoff'].median()).abs().argsort()[:1]].index.year}
    
    plot_yr = {'High': ann_tmp.idxmax(axis=0)['runoff'].year}
    #plot_yr = {'Median': ann_tmp.ix[(ann_tmp.runoff-ann_tmp['runoff'].median()).abs().argsort()[:1]].index.year}

    outpdf = None    # outpdf is instantiated in the loop below

    for kk, yy in plot_yr.iteritems():
        # Reset date range to maximum year
        st = datetime.datetime(yy-1,10,1)
        en = datetime.datetime(yy,9,30)

        plot_data = statvar_data[st:en]
        orig_plot_data = statvar_orig_data[st:en]
        
        maxrows = 5
        maxcols = 1

        cc = maxrows * maxcols    # used to track plots on a page
        newYear = True

        for ii,vv in enumerate(plotvars):
            if cc == maxrows * maxcols:
                if outpdf is None:
                    #outpdf = PdfPages(pdf_filename)
                    newYear = False
                elif newYear:
                    newYear = False
                else:
                    plt.suptitle('basin: %s (%s)\nrunid: %s (%s)\nYear: %d (%s flow)' % (basinid, site_name, runid, modelrunid, yy, kk), fontsize=title_size)
                    #outpdf.savefig()

                cc = 0

                fig, axes = plt.subplots(nrows=maxrows, ncols=maxcols, figsize=(17,11), sharex=True)
                ax = axes.flatten()

                # Plot obs/sim plot first on new page
                stuff = ax[cc].plot(plot_data.index.to_pydatetime(), plot_data[obs_var], color='grey', label=obs_var)
                stuff = ax[cc].plot(plot_data.index.to_pydatetime(), plot_data[sim_var], color='red', label=sim_var)
                stuff = ax[cc].plot(orig_plot_data.index.to_pydatetime(), orig_plot_data[sim_var], ':', color='red', label='original')

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
                try:
                    ma = min(plot_data[vv].min(), orig_plot_data[vv].min())
                    mi = max(plot_data[vv].max(), orig_plot_data[vv].max())
                    mm = [ma, mi]
                    set_ylim(ax[cc], mm, 0.02)
                except:
                    set_ylim(ax[cc], plot_data[vv], 0.02)


            stuff = ax[cc].plot(plot_data.index.to_pydatetime(), plot_data[vv], color='green', label=vv)
            stuff = ax[cc].plot(orig_plot_data.index.to_pydatetime(), orig_plot_data[vv], ':', color='green')
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

        # Remove any unused subplots
        while cc < maxrows * maxcols:
            fig.delaxes(ax[cc])
            cc += 1


# %%
sim_var = cfg.get_value('sim_var')
obs_var = cfg.get_value('obs_var')
site_name = ''
title_size = 13

bigplot()

# %%
