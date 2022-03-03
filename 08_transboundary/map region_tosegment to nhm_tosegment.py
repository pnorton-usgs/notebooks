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
#     display_name: Python [conda env:bandit]
#     language: python
#     name: conda-env-bandit-py
# ---

# %% [markdown]
# ## Create a csv file which maps regional tosegment and region ID to the nhm tosegment

# %%
def data_it(filename):
    """Get iterator to a parameter db file.

    :returns: iterator
    """

    # Read the data
    fhdl = open(filename)
    rawdata = fhdl.read().splitlines()
    fhdl.close()
    return iter(rawdata)


# %%
REGIONS = ['r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09',
           'r10L', 'r10U', 'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18']

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb'
local_param = 'tosegment'
nhm_param = 'tosegment_nhm'

# %%
# Read the tosegment parameters 
tmp_data = []

for rr in REGIONS:
    # Read parameter information
    cdir1 = '{}/{}/{}'.format(workdir, local_param, rr)
    cdir2 = '{}/{}/{}'.format(workdir, nhm_param, rr)

    it1 = data_it('{}/{}.csv'.format(cdir1, local_param))
    next(it1)    # Skip the header row

    it2 = data_it('{}/{}.csv'.format(cdir2, nhm_param))
    next(it2)    # Skip the header row

    # Read the parameter values
    for rec1, rec2 in zip(it1, it2):
        
        idx1, val1 = rec1.split(',')
        idx2, val2 = rec2.split(',')

        tmp_data.append([rr[1:], val1, val2])


# %%
# Write the csv file
outdir = '/Users/pnorton/Projects/National_Hydrology_Model/GIS'
outhdl = open(f'{outdir}/tosegment_region_to_nhm.csv', 'w')

outhdl.write('region,tosegment_region,tosegment_nhm\n')
for xx in tmp_data:
    outhdl.write(f'\'{xx[0]}\',{xx[1]},{xx[2]}\n')
    
outhdl.close()

# %%
