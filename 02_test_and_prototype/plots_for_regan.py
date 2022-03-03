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
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/tmp'
pdf_filename = '%s/06191500_old_v_prms402_b.pdf' % workdir

st = datetime.datetime(1981,10,1)
en = datetime.datetime(1985,9,30)

# %%
orig = prms.statvar('%s/daymet.statvar_orig' % workdir).data
test = prms.statvar('%s/daymet.statvar_test4' % workdir).data

orig_mon = orig[st:en].resample('M', how='mean')
test_mon = test[st:en].resample('M', how='mean')

# %%
print orig.columns
#print test.head()

# %%
outpdf = PdfPages(pdf_filename)

fig, axes = plt.subplots(nrows=19, ncols=1, figsize=(17,100))
ax = axes.flatten()

for ii, cc in enumerate(orig.columns):
    #print cc
    #ax[ii].scatter(orig.loc[:,cc], test.loc[:,cc])
    ax[ii].plot(orig_mon.loc[:,cc].index.to_pydatetime(), orig_mon.loc[:,cc], linewidth=1.25, color='black', alpha=.8, label='orig')
    ax[ii].plot(test_mon.loc[:,cc].index.to_pydatetime(), test_mon.loc[:,cc], linewidth=0.25, color='red', alpha=1., label='test')
    
    ax[ii].set_title(cc, fontsize=12)
    
outpdf.savefig()
outpdf.close()

# %%
