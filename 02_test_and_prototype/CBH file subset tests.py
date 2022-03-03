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

# %%

# %%
workdir = "/Users/pnorton/USGS/Projects/National_Hydrology_Model/regions/r10U/input/cbh"
filename = '%s/daymet_1980_2010_tmin.cbh' % workdir
missing = [-99.0, -999.0]

infile = open(filename, 'r')
fheader = ''

for ii in range(0,3):
    line = infile.readline()
    
    if line[0:4] in ['prcp', 'tmax', 'tmin']:
        # Change the number of HRUs included to one
        numhru = int(line[5:])
        fheader += line[0:5] + ' 1\n'
    else:
        fheader += line
print fheader
print 'numhru:', numhru


# %%
# Read in the CBH data for the HRU we want to extract
hruindex = 1    # one-based hru index

df1 = pd.read_csv(infile, sep=' ', skipinitialspace=True,
                  #usecols=[0, 1, 2, 3, 4, 5, hruindex+5],
                  header=None)
# df1 = pd.read_csv(infile, sep=r"\s*", engine='python',
#                   skiprows=3, usecols=[0, 1, 2, 3, 4, 5, hruindex+6],
#                   header=None)
infile.close()

df1.head(10)

# %%
df1.loc[:,[0,1,2,8]]

# %%
# Write the subsetted CBH data out
outfile = open('crap.cbh', 'w')
outfile.write(fheader)

df1.to_csv(outfile, sep=' ', float_format='%0.4f', header=False, index=False)
outfile.close()

# %%

# %%
workdir = "/Users/pnorton/Projects/National_Hydrology_Model/tmp"
filename = '%s/daymet_1980_2011_prcp.cbh' % workdir
missing = [-99.0, -999.0]

infile = open(filename, 'r')
fheader = ''

for ii in range(0,3):
    line = infile.readline()
    
    if line[0:6] in ['precip', 'tmax', 'tmin']:
        # Change the number of HRUs included to one
        numhru = int(line[7:])
        fheader += line[0:5] + ' 1\n'
    else:
        fheader += line
print fheader
print 'numhru:', numhru

# %%
df1 = pd.read_csv(infile, sep=' ', skipinitialspace=True,
                  #usecols=[0, 1, 2, 3, 4, 5, hruindex+5],
                  header=None)
# df1 = pd.read_csv(infile, sep=r"\s*", engine='python',
#                   skiprows=3, usecols=[0, 1, 2, 3, 4, 5, hruindex+6],
#                   header=None)
infile.close()

df1.head(10)

# %%
# Check for precip values less than 0.001
df2 = df1[df1.iloc[:,6:] < 0.001]
df2.sum().sum()

# %%

# %%
