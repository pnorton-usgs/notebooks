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
from pyPRMS import Streamflow as sf

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20180927_ARPO'
filename = '{}/sf_data'.format(workdir)

sf_data = sf.Streamflow(filename, verbose=True)

# %%
df = sf_data.data
df.head()

# %%
aa = df.iloc[:,1976]
aa

# %%
aa.info()

# %%
infile = open(filename, 'r')
rawdata = infile.read().splitlines()
infile.close()

it = iter(rawdata)
hdr_count = 0

# Skip to the data section
for line in it:
    hdr_count += 1
    if line[0:10] == '##########':
        hdr_count += 1  # plus one for good measure
        break
        
print('Data starts at line: {}'.format(hdr_count))

# %%
cline = 0
pcols = []
plines = []
pdata = []

for line in it:
    flds = line.split()
    cline += 1
    
    for idx, ff in enumerate(flds):
        try:
            aa = float(ff)
        except ValueError:
            pcols.append(idx)
            plines.append(cline)
            pdata.append(ff)
#             print('Error in line {}, column {}, data = {}'.format(cline, idx, ff))

print('Bad columns: {}'.format(set(pcols)))
print('Bad lines: {}'.format(set(plines)))
print('Bad data: {}'.format(set(pdata)))

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/paramdb_v2'
filename = '{}/poi_gage_id.csv'.format(workdir)

infile = open(filename, 'r')
rawdata = infile.read().splitlines()
infile.close()

it = iter(rawdata)
it.next()

gageids = {}

for line in it:
    flds = line.split(',')
    
    if flds[1] in gageids:
        print '{} ({}): Duplicate gageid {} ({})'.format(flds[1], gageids[flds[1]], flds[1], flds[0])
    else:
        gageids[flds[1]] = flds[0]
        

# %%
