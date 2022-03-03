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
import glob
import os


# %%
def _data_it(filename):
    """Get iterator to a parameter db file.

    :returns: iterator
    """

    # Read the data
    fhdl = open(filename)
    rawdata = fhdl.read().splitlines()
    fhdl.close()
    return iter(rawdata)


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


# %%
workdir_tb = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/tbparamdb'
workdir_v1 = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_CONUS'
workdir_out = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/paramdb_merged'

# %%
tb_files = []

file_it = glob.glob('{}/*.csv'.format(workdir_tb))
for kk in file_it:
    tb_files.append(kk)
tb_files.sort()

# %%
for ff in tb_files:
    filename = os.path.basename(ff)
    if filename in ['nhm_to_GFv11_HRU.csv', 'nhm_to_GFv11_SEG.csv', 'obsout_segment.csv', 'seg_id_nhm.csv']:
        continue
    
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
#     if cname == 'hru_segment_nhm':
#         tb_type_nca = read_csv('{}/../type_nca.csv'.format(workdir_tb))
        
    try:
        it = _data_it(ff)
        next(it)  # Skip the header row
    except IOError:
        print('Skipping parameter: {}. File does not exist.'.format(ff))
        continue
            
    # Read the parameter values
#     if cname == 'hru_segment_nhm':
#         cidx = 0
        
#         for rec in it:
#             idx, val = rec.split(',')
            
#             if tb_type_nca[cidx] == 1:
#                 print('  - reset idx:{} to zero'.format(cidx))
#                 tmp_data.append(0)
#             else:
#                 tmp_data.append(val)
#             cidx += 1
#     else:
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
len(tb_type_nca)

# %%
idx

# %%
