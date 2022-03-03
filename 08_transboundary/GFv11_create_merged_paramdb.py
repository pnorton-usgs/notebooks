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

# %%
import glob
import os

# %%
workdir_tb_ids = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/2019-11-25_work'
workdir_tb_parameters = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/tbparamdb'

workdir_v1_parameters = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_CONUS'


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
            data.append(int(rec.split(',')[1]))
        except ValueError:
            data.append(int(float(rec.split(',')[1])))
    return data

def read_to_dict(filename):
    fhdl = open(filename)
    rawdata = fhdl.read().splitlines()
    fhdl.close()
    it = iter(rawdata)

    data = {}

    # Skip header
    next(it)

    # Read the parameter values
    for rec in it:
        rec_ar = rec.split(',')
        
        try:
            data[int(rec_ar[0])] = int(rec_ar[1])
        except ValueError:
            data[int(rec_ar[0])] = float(rec_ar[1])
    return data    


# %%
version_ids = read_to_dict('{}/GFv11_seg_version.csv'.format(workdir_tb_ids))
seg_id_nhm = read_to_dict('{}/GFv11_seg_id_nhm.csv'.format(workdir_tb_ids))



# %%

# %%
# Get a list of just the GFv11 segment ids
seg_id_list = sorted(list(seg_id_nhm.keys()))

if len(seg_id_list) != seg_id_list[-1]:
    print('Size of list does not match id numbering')

# %%
tb_files = []

file_it = glob.glob('{}/*.csv'.format(workdir_tb_parameters))
for kk in file_it:
    tb_files.append(kk)
tb_files.sort()

# %%
for ff in tb_files:
    filename = os.path.basename(ff)
    
    cname = os.path.basename(os.path.splitext(ff)[0])
    print(cname)
    
    v1_paramfile = '{}/{}'.format(workdir_v1_parameters, filename)
    
    for ii in seg_id_list:
        if version_ids[ii] == 1:
            print('Read version 1 parameters')
        else:
            print('Read version 1.1 parameters')

    # %%
    filename = os.path.basename(ff)
    
    cname = os.path.basename(os.path.splitext(ff)[0])
    print(cname)
    
    # Read the parameter data
    tmp_data = []
    
    # ===================================================================
    # Read the NHM v1 parameter first
    it = _data_it('{}/{}'.format(workdir_v1, filename))
    next(it)
        
    for rec in it:
        idx, val = rec.split(',')
        tmp_data.append(val)
    
    # ===================================================================
    # Now add the transboundary parameter information
    if cname == 'hru_segment_nhm':
        tb_type_nca = read_csv('{}/../type_nca.csv'.format(workdir_tb))
        
    try:
        it = _data_it(ff)
        next(it)  # Skip the header row
    except IOError:
        print('Skipping parameter: {}. File does not exist.'.format(ff))
        continue
            
    # Read the parameter values
    if cname == 'hru_segment_nhm':
        cidx = 0
        
        for rec in it:
            idx, val = rec.split(',')
            
            if tb_type_nca[cidx] == 1:
                print('  - reset idx:{} to zero'.format(cidx))
                tmp_data.append(0)
            else:
                tmp_data.append(val)
            cidx += 1
    else:
        for rec in it:
            idx, val = rec.split(',')
            tmp_data.append(val)

    # ===================================================================
    # Write the merged parameter to the output
    outfile = open('{}/{}'.format(workdir_out, filename), 'w')
    outfile.write('$id,{}\n'.format(cname))
    
    for idx,val in enumerate(tmp_data):
        outfile.write('{},{}\n'.format(idx+1, val))
        
    outfile.close()

# %%
