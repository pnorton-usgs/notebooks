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
import numpy as np
import os
import pandas as pd
import prms_lib as prms

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet'
thevar = 'prcp'   # one of tmin, tmax, prcp

# Create the output hdf5 file
hdf = pd.HDFStore('{0}/r10U_{1}.h5'.format(workdir, thevar))

# %%
# reload(prms)
# Open the CBH file and read it in
filename = '{0}/r10U_{1}.cbh.gz'.format(workdir, thevar)
% time df1 = pd.read_csv(filename, sep=' ', skipinitialspace=True, skiprows=3, engine='c', header=None)
# % time df1 = prms.read_cbh(filename)

# %%
os.path.splitext(filename)

# %%
# Add the data to the hdf5 file

# % time hdf.put('prcp', df1, format='table', data_columns=True)
# % time hdf.append(thevar, df1)
% time hdf.put(thevar, df1)

# %%

# %%
hdf[thevar].head()

# %%
# Close the hdf5 file
hdf.close()

# %%
h5file = '{0}/r10U_{1}.h5'.format(workdir, thevar)
% time hdf2 = pd.read_hdf(h5file)

# %%
print hdf2.info()

# %%
reload(prms)
import split_by_HRU as sbh
reload(sbh)

# %%
sbh.full_cbh_subset('/media/scratch/PRMS/regions/rGCPO/daymtTmin1980_2014.txt', 'crap', 'GCPO', 'tmin', 2, hdf=True)

# %%
df1.

# %%
hdf

# %%
