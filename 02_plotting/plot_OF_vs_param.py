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
# Uncomment the following line to write to pdf without opening a window
#matplotlib.use('Agg')

# %matplotlib inline
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm

import prms_lib as prms
import pandas as pd
import numpy as np
import math as mth
import datetime
import re

# %%
#workdir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/tmp_data'
#mocom_file = '%s/mocom_t2.txt' % workdir

# Read in mocom file and use regex to specify variable whitespace between fields
#mocom = pd.read_csv(mocom_file, sep=r"\s*", engine='python')
#mocom['objfcn'] = mocom['test0']*0.3 + mocom['test1']*0.7
#mocom.head()

# %%
# Setup model run information
templatedir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/PRMS_master'
basedir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t3'
basinid = '06191500'
runid = '2015-03-30_1538'

workdir = '%s/%s/runs/%s' % (basedir, basinid, runid)
modeldir = '%s/%s' % (templatedir, basinid)

# Plot last 25% of the generations? If false return all
top_qtr = True

# %% [markdown]
# ### Read in parameter file used for MOCOM calibration

# %%
#workdir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/tmp_data'
limits_file = '%s/param_limits.txt' % workdir
have_limits = False

try:
    limits = pd.read_csv(limits_file, header=None, names=['parameter','maxval','minval'], sep=r"\s*", engine='python')

    params = limits['parameter'].tolist()
    have_limits = True
    print limits
except IOError:
    print "ERROR: %s does not exist.\nLimits will not be plotted." % limits_file
    

# %% [markdown]
# ### Read input parameter file to get the default parameter values

# %%
# Read the initial input parameters values and compute the mean for each one
initparams = prms.parameters('%s/default.params' % modeldir)

# Build list of initial/default values for the calibration parameters
initvals = []
for vv in params:
    initvals.append(np.mean(initparams.get_var(vv)['values']))

# %% [markdown]
# ### Plot objective-function(s) versus parameters

# %%
# Plot all the generation results
if top_qtr:
    pdf_filename = '%s/%s/pdf/%s_%s_OF_v_params_lastQTR.pdf' % (basedir, basinid, basinid, runid)
else:
    pdf_filename = '%s/%s/pdf/%s_%s_OF_v_params.pdf' % (basedir, basinid, basinid, runid)

mocom_file = '%s/optim_fixed.log' % workdir

# Read in mocom file and use regex to specify variable whitespace between fields
mocom = pd.read_csv(mocom_file, sep=',')
objfcns = [col for col in mocom.columns if 'test' in col]
#of_names = {'test0':'NRMSE(lowflows)',
#            'test1':'NRMSE(highflows)'}

of_names = {'test0':'NRMSE(monthly)',
            'test1':'NRMSE(annual)',
            'test2':'NRMSE(mean monthly)'}

# --------------------------------------------------------------
# Create an array of colors based on the generation number field
maxgen = max(mocom['gennum'])
mingen = min(mocom['gennum'])

if top_qtr:
    # Just work with the last 25% of the generations
    mocom = mocom[mocom['gennum'] > int(maxgen * 0.75)]

# Normalize generation numbers to 0 to 1
nvals = (mocom['gennum'] - mingen) / (maxgen - mingen)

# Generate the colors array
colors = cm.rainbow(nvals)
# --------------------------------------------------------------

# Setup output to a pdf file
outpdf = PdfPages(pdf_filename)
numrows = int(round(len(params) / 4. + 0.5))
print numrows 

for oo in objfcns:
    miny = min(mocom[oo])
    maxy = max(mocom[oo])
    
    fig, axes = plt.subplots(nrows=numrows, ncols=4, figsize=(20,5*numrows))
    ax = axes.flatten()

    for ii,pp in enumerate(params):
        # Range of x-axis for parameter
        if have_limits:
            maxx = limits.iloc[ii,1]
            minx = limits.iloc[ii,2]
        else:
            maxx = max(mocom[pp])
            minx = min(mocom[pp])
        
        # Shut off automatic offsets for the y-axis
        ax[ii].get_yaxis().get_major_formatter().set_useOffset(False)

        # Set the limits on the x-axis
        xxrange = maxx - minx
        xpad = abs(xxrange * 0.1)
        ax[ii].set_xlim([minx - xpad, maxx + xpad])
        
        # Set the limits on the y-axis
        yrange = maxy - miny
        ypad = abs(yrange * 0.1)
        ax[ii].set_ylim([miny - ypad, maxy + ypad])

        ax[ii].scatter(mocom[pp], mocom[oo], color=colors, alpha=0.5)
        
        # Plot parameter range limits
        if have_limits:
            ax[ii].plot([maxx, maxx], [miny - ypad, maxy + ypad], 'k--', color='red')
            ax[ii].plot([minx, minx], [miny - ypad, maxy + ypad], 'k--', color='red')
        
        ax[ii].plot([initvals[ii], initvals[ii]], [miny, maxy], markeredgecolor='black', markerfacecolor='yellow', color='grey', marker='D')

        # Re-plot the final pareto set 
        ax[ii].scatter(mocom[pp].loc[mocom['gennum'] == maxgen], mocom[oo].loc[mocom['gennum'] == maxgen], color='black', marker='x')
        ax[ii].set_title(pp, fontsize=12)

    plt.suptitle('Basin: %s (%s)\nObjective Function vs. Parameters\n%s' % (basinid, runid, of_names[oo]), fontsize=14)
    #plt.subplots_adjust(top=0.75)
    #plt.subplots_adjust(hspace = 0.3)
    outpdf.savefig()
    
    
outpdf.close()

# %%

# %%
