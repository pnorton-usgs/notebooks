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

import prms_lib as prms
import pandas as pd
import numpy as np
import math as mth
import datetime
import lmoments as lmom
import statsmodels.formula.api as smf
reload(prms)

# %%
workdir = 'test_data'
obs_streamflow_file = '%s/X.data' % workdir
sim_streamflow_file = '%s/X.statvar' % workdir

st = datetime.datetime(1982,1,1)
en = datetime.datetime(2009,12,31)

# %%
obs_streamflow = prms.streamflow(obs_streamflow_file).data
#obs_subset = obs_streamflow.resample('D').iloc[:,0]
obs_subset = obs_streamflow[st:en].iloc[:,0]
obs_subset *= 0.02832
obs_subset.head()

# %%
luca_sim = prms.statvar(sim_streamflow_file).data
luca_sim.drop(['rec'], axis=1, inplace=True)
luca_subset = luca_sim[st:en].iloc[:,2]

# %%
luca_subset.head()

# %%
# Compute the first four L-moments
#lmom.samlmu(obs_subset.dropna().tolist(), 4)

# NOTE: The mag7.f code from Lauren doesn't ignore missing data so the l-moments are skewed
#       Is this what she is wanting?
#       Will need to divide the 2nd moment by the first to get the Coeff of variation
lmom.samlmu(obs_subset.tolist(), 4)

# %%
# Compute seasonally adjusted streamflow
obs_tmp = obs_subset.groupby(obs_subset.index.month).transform(lambda x: x - x.mean())

# %%
obs_sdev = obs_tmp.std()
print 'obs_sdev:', obs_sdev
obs_ltmean = obs_tmp.mean()
print 'obs_ltmean:', obs_ltmean

obs_adj = (obs_tmp - obs_ltmean) / obs_sdev
print obs_adj.head()

# %%
# NOTE: missing data is included in the computation which skews the result
obs_adj.autocorr()

# %%
dat_doy = obs_subset.index.dayofyear

tt = 2 * mth.pi * dat_doy
dat_cos = np.cos(tt)
dat_sin = np.sin(tt)
# hh is the data

outdata = {}
outdata['data'] = obs_adj
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
res.params

# %%
amplitude = mth.sqrt(res.params[1]^2 + res.params[2]^2)

# %%
res.params[1]

# %%
