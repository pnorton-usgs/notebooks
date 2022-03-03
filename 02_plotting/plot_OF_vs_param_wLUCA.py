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
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
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

# %% [markdown]
# ### Read in parameter file used for MOCOM calibration

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/tmp_data'
limits_file = '%s/param_limits.txt' % workdir

limits = pd.read_csv(limits_file, header=None, names=['parameter','maxval','minval'], sep=r"\s*", engine='python')

params = limits['parameter'].tolist()
print limits

# %% [markdown]
# ### Process LUCA summary file for comparison

# %%
# Deal with LUCA summary file (ugh)

luca_file = '%s/test1_summary3.txt' % workdir
luca_out = '%s/test1_summary3_fixed.txt' % workdir

infile = open(luca_file, 'r')
rawdata = infile.read().splitlines()
infile.close()
it = iter(rawdata)

outheader = ['parameter', 'round', 'mean']
outdata = []

outfile = open(luca_out, 'w')
outfile.write('parameter,round,mean\n')

for line in it:
    if line[0:29] == '>>> Objective Function Values':
        while next(it)[0:5] != 'Round':
            pass
        objfcn = re.findall(r"[\w.-]+", next(it))
        for ii,xx in enumerate(objfcn):
            outfile.write('objfcn,%d,%s\n' % (ii+1, xx))
        
    if line[0:14] == '>>> Parameter:':
        var = line.split(' ')[2]
        
        while True: 
            x = next(it)
            if x[0:4] == 'Mean':
                break
        
        words = re.findall(r"[\w.-]+", x) 
        #print words
        
        for xx in range(2,8):
            outfile.write('%s,%d,%s\n' % (var,xx-1,words[xx]))
        #outdata.append([var, words[2], words[3], words[4], words[5], words[6], words[7]])
outfile.close()

# %% [markdown]
# ### Read input parameters from last LUCA round that were used for MOCOM calibration

# %%
modeldir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/PRMS_master/06191500'

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
pdf_filename = 'OF_params_scatter.pdf'

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/tmp_data'
mocom_file = '%s/optim_fixed.log' % workdir

# Read in mocom file and use regex to specify variable whitespace between fields
mocom = pd.read_csv(mocom_file, sep=',')
mocom['objfcn'] = mocom['test0']*0.3 + mocom['test1']*0.7

# Read in 'fixed' version of LUCA summary file for comparison to MOCOM
luca_file = '%s/test1_summary3_fixed.txt' % workdir
luca = pd.read_csv(luca_file)
luca.head()
lucaPivot = luca.pivot(index='round', columns='parameter', values='mean')

# get minimum/maximum objective functions values for mocom and luca
if min(lucaPivot['objfcn']) < min(mocom['objfcn']):
    miny = min(lucaPivot['objfcn'])
else:
    miny = min(mocom['objfcn'])
    
if max(lucaPivot['objfcn']) > max(mocom['objfcn']):
    maxy = max(lucaPivot['objfcn'])
else:
    maxy = max(mocom['objfcn'])

# --------------------------------------------------------------
# Create an array of colors based on the generation number field
maxgen = max(mocom['gennum'])
mingen = min(mocom['gennum'])

# Normalize generation numbers to 0 to 1
nvals = (mocom['gennum'] - mingen) / (maxgen - mingen)

# Generate the colors array
colors = cm.rainbow(nvals)
# --------------------------------------------------------------

# Setup output to a pdf file
outpdf = PdfPages(pdf_filename)

numrows = int(round(len(params) / 4.))
fig, axes = plt.subplots(nrows=numrows, ncols=4, figsize=(20,15))
ax = axes.flatten()

for ii,pp in enumerate(params):
    # Shut off automatic offsets for the y-axis
    ax[ii].get_yaxis().get_major_formatter().set_useOffset(False)
    
    # Set the limits on the x-axis
    xlim_min = limits.iloc[ii,2] - abs(limits.iloc[ii,2] * .05)
    if abs(xlim_min) < 0.1:
        xlim_min -= (limits.iloc[ii,1] - limits.iloc[ii,2]) / 10.0
    xlim_max = limits.iloc[ii,1] + abs(limits.iloc[ii,1] * 0.05)
    ax[ii].set_xlim([xlim_min, xlim_max])
    
    # Set the limits on the y-axis
    ax[ii].set_ylim([miny - abs(miny*.02), maxy + abs(maxy * .02)])

    maxX = limits.iloc[ii,1]
    minX = limits.iloc[ii,2]
    
    ax[ii].scatter(mocom[pp], mocom['objfcn'], color=colors, alpha=0.5)
    ax[ii].plot([maxX, maxX], [miny, maxy], 'k--', color='red')
    ax[ii].plot([minX, minX], [miny, maxy], 'k--', color='red')
    try:
        ax[ii].scatter(lucaPivot[pp], lucaPivot['objfcn'], color='green', marker=u's')
    except:
        # LUCA may not have the parameter, but we don't care; move on.
        pass
    
    ax[ii].plot([initvals[ii], initvals[ii]], [miny, maxy], markeredgecolor='black', markerfacecolor='yellow', color='grey', marker='D')

    # Re-plot the final pareto set 
    ax[ii].plot(mocom[pp].loc[mocom['gennum'] == maxgen], mocom['objfcn'].loc[mocom['gennum'] == maxgen], color='black', marker='x')
    ax[ii].set_title(pp, fontsize=12)

plt.suptitle('Objective Function vs. Parameters', fontsize=18)
outpdf.savefig()
outpdf.close()

# %%

# %%
