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
import pandas as pd
import datetime

def dparse(yr, mo, dy, hr, minute, sec):
    # Date parser for working with the date format from *.statvar files
    
    # Convert to integer first
    yr, mo, dy, hr, minute, sec = [int(x) for x in [yr, mo, dy, hr, minute, sec]]
    
    dt = datetime.datetime(yr, mo, dy, hr, minute, sec)
    return dt


# %%
# Setup model run information
templatedir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/PRMS_master'
basedir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t3'
basinid = '06191500'
runid = '2015-03-30_1538'

workdir = '%s/%s/runs/%s' % (basedir, basinid, runid)
modeldir = '%s/%s' % (templatedir, basinid)

# %%
filename = 'X.statvar'

infile = open(filename, 'r')

# The first line gives the number of variables that follow
numvars = int(infile.readline())
print "Number of variables: %d" % (numvars)

# The next numvar rows contain a variable name followed by a number which
# indicates the number of columns used by that variable.
# The relative order of the variable in the list indicates the column
# the variable data is found in.
varcol = {}
varlist = []

for rr in range(0,numvars):
    row = infile.readline()
    fields = row.rstrip().split(' ')

    varcol[fields[0]] = [rr, int(fields[1]), [] ]
    varlist.append(fields[0])

print varcol





# %%
# Read the remainder of the file in
# Because we pass the read_csv routine an open file handle instead of a filename
# it will start reading where the iterator was at after reading the variable names.

# The first 7 columns are [record year month day hour minute seconds]
thecols = ['rec', 'year', 'month', 'day', 'hour', 'min', 'sec']

# Add the remaining columns to the list
for xx in varlist:
    thecols.append(xx)

# Use pandas to read the data in from the remainder of the file
# We use a custom date parser to convert the date information to a datetime
thedata = pd.read_csv(infile, sep=r"\s+", header=None, names=thecols, 
                      parse_dates={'thedate': ['year', 'month', 'day', 'hour', 'min', 'sec']}, 
                      date_parser=dparse, index_col='thedate')


# %%

# %%
