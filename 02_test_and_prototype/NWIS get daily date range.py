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

# %% [markdown]
# ####Reads reference gages (probably gages-II) file and get the available date range for daily streamflow from NWIS for each gage in the file.

# %%
import pandas as pd
import numpy as np
import urllib2
import re
from time import strftime

# %% [markdown]
# Read the reference gages file

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing'
ref_file = '%s/refPARKER_gages' % workdir

colnames = ['region', 'seq', 'siteno', 'area', 'ns1', 'ns2', 'ref']
types = {'region': np.dtype(str), 'seq': np.dtype(int), 'siteno': np.dtype(str),
         'area': np.dtype(float), 'ns1': np.dtype(float), 'ns2': np.dtype(float),
         'ref': np.dtype(str)}

refdata = pd.read_csv(ref_file, sep=r"\s+", header=None, names=colnames, dtype=types)

# Create a list of the reference gage sites
thesites = ','.join(refdata['siteno'].tolist())

# %% [markdown]
# Download the daily streamflow date range for each streamgage

# %%
# Setup the NWIS services URL stuff
base_url = 'http://waterservices.usgs.gov/nwis'
stn_url = '%s/site/?format=rdb&sites=%s&seriesCatalogOutput=true&outputDataTypeCd=dv&parameterCd=00060' % (base_url, thesites)

streamgagesPage = urllib2.urlopen(stn_url)

t1 = re.compile('^#.*$\n?', re.MULTILINE)   # remove comment lines
t2 = re.compile('^5s.*$\n?', re.MULTILINE)  # remove field length lines

# Get the list of streamgages within the specified region
streamgagesFromREST = streamgagesPage.read()

# Strip the comment lines and field length lines from the result
streamgagesFromREST = t1.sub('', streamgagesFromREST, 0)
streamgagesFromREST = t2.sub('', streamgagesFromREST, 0)


# %% [markdown]
# Output date range for just the mean daily value for each streamgage

# %%
aa = streamgagesFromREST.splitlines()

outfile = open('crap.tab', 'w')

for line in aa:
    if line.split('\t')[13] in ['parm_cd', '00060'] and line.split('\t')[14] in ['stat_cd', '00003']:
        outfile.write(line + '\n')

outfile.close()

# %% [markdown]
# Read the file created above, join with the original reference data, and write out to new file

# %%
siteinfo = pd.read_csv('crap.tab', sep='\t', usecols=['site_no','begin_date','end_date'],
                      dtype={'site_no': np.dtype(str)})
comb = pd.merge(refdata, siteinfo, how='left', left_on='siteno', right_on='site_no')
comb.head()

# %%
comb.to_csv('%s/ref_with_dates.csv' % workdir, index=False)

# %%
C = float(10**-10)
OM = 0.01138749304761900
NM = 0.00100783

numer = ((NM + C) * (0.010048891 + C)) / (OM+C) - C
print numer

# %%
