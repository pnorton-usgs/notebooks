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

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
workdir_tb = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/tbparamdb'


# %%
def read_csv(filename):
    fhdl = open(filename)
    rawdata = fhdl.read().splitlines()
    fhdl.close()
    it = iter(rawdata)

    data = []

    # Skip header
    next(it)

    # Read the parameter values
    for rec in it:
        try:
            data.append(int(rec.split(',')[0]))
        except ValueError:
            data.append(int(float(rec.split(',')[0])))
    return data

def read_csv2(filename):
    fhdl = open(filename)
    rawdata = fhdl.read().splitlines()
    fhdl.close()
    it = iter(rawdata)

    data = []

    # Skip header
    next(it)

    # Read the parameter values
    for rec in it:
        try:
            data.append(int(rec.split(',')[1]))
        except ValueError:
            data.append(int(float(rec.split(',')[1])))
    return data

def check_skips(var):
    prev = var[0] - 1

    for idx, val in enumerate(var):
        if prev != val - 1:
            print(f'{idx}: prev = {prev}, cval = {val-1}')
        prev += 1


# %%
# 2020-02-26 PAN: check seg and hru ids for skips in v11e geofabric (nhru_v11 and nsegment_v11 ids)
# NOTE: This doesn't check for monotonic ordering of the native file; only for skips in ID numbering
workdir = "/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/2019-11-25_work"
seg_file = f'{workdir}/SEG_mapping.csv'
hru_file = f'{workdir}/HRU_mapping.csv'

segs = read_csv(seg_file)
hrus = read_csv(hru_file)

segs.sort()
hrus.sort()

print('segments')
check_skips(segs)

print('HRUs')
check_skips(hrus)

# %%
# 2020-02-28 PAN: Check for ordering of nhm_seg and nhm_id parameters
# NOTE: This DOES check for skips in the original file order of IDs
workdir = "/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/tbparamdb"
seg_file = f'{workdir}/nhm_seg.csv'
hru_file = f'{workdir}/nhm_id.csv'

segs = read_csv2(seg_file)
hrus = read_csv2(hru_file)

print('segments')
check_skips(segs)

print('HRUs')
check_skips(hrus)

# %%
# 2020-02-28 PAN: Check for ordering of nhm_seg and nhm_id parameters
# NOTE: This DOES check for skips in the original file order of IDs
workdir = "/Users/pnorton/tmp/tmp_paramdb"
seg_file = f'{workdir}/nhm_seg.csv'
hru_file = f'{workdir}/nhm_id.csv'

segs = read_csv2(seg_file)
hrus = read_csv2(hru_file)

print('segments')
check_skips(segs)

print('HRUs')
check_skips(hrus)

# %% [markdown]
# ### Other checks

# %%
transp_beg = read_csv('{}/op_flow_thres.csv'.format(workdir_tb))

prev = transp_beg[0] - 1

for idx, val in enumerate(transp_beg):
    if prev != val - 1:
        print('{}: prev = {}, cval = {}'.format(idx, prev, val-1))
    prev += 1
    
print('first: {}'.format(transp_beg[0]))
print('last: {}'.format(transp_beg[-1]))

# %%
nhm_id = read_csv2('{}/nhm_id.csv'.format(workdir_tb))

prev = nhm_id[0] - 1

for idx, val in enumerate(nhm_id):
    if prev != val - 1:
        print('{}: prev = {}, cval = {}'.format(idx, prev, val))
        prev = val
    else:
        prev += 1
    
print('first: {}'.format(nhm_id[0]))
print('last: {}'.format(nhm_id[-1]))

# %%
nhm_seg = read_csv2('{}/nhm_seg.csv'.format(workdir_tb))

prev = nhm_seg[0] - 1

for idx, val in enumerate(nhm_seg):
    if prev != val - 1:
        print('{}: prev = {}, cval = {}'.format(idx, prev, val))
        prev = val
    else:
        prev += 1
    
print('first: {}'.format(nhm_seg[0]))
print('last: {}'.format(nhm_seg[-1]))

# %%
param = read_csv('{}/nhm_to_GFv11_HRU.csv'.format(workdir_tb))

print(set(nhm_id) ^ set(param))
print(nhm_id == param)

# %%

# %%
# Segment parameter
param = read_csv('{}/nhm_to_GFv11_SEG.csv'.format(workdir_tb))

print(set(nhm_seg) ^ set(param))
print(nhm_seg == param)

# %%

# %%

# %%
nhm_id = read_csv2('{}/nhm_id.csv'.format(workdir_tb))
type_nca = read_csv2('{}/../type_nca.csv'.format(workdir_tb))

for idx, nhmid in enumerate(nhm_id):
    if type_nca[idx] == 1:
        print(nhmid)

# %%
coastal_hru = read_csv2('{}/../coastal_hru.csv'.format(workdir_tb))

for idx, nhmid in enumerate(nhm_id):
    if coastal_hru[idx] == 1:
        print(nhmid)

# %%
# Non-route HRUs that are being removed from the NHM
# hru_repl_v1.hru_segment.csv
fhdl = open('{}/../hru_repl_v1.hru_segment.csv'.format(workdir_tb))
rawdata = fhdl.read().splitlines()
fhdl.close()
it = iter(rawdata)

# data = []

# Skip header
next(it)

# Read the parameter values
for rec in it:
    nhmid, seg = rec.split(',')
    
    if int(seg) == 0:
        print(nhmid)
        

# %%
