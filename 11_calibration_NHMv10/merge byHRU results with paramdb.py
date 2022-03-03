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
#     display_name: Python [conda env:bandit_38]
#     language: python
#     name: conda-env-bandit_38-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import glob
import numpy as np
import os
import pandas as pd
import sys

from collections import Counter

from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ControlFile import ControlFile
from pyPRMS.ValidParams import ValidParams

# %%
base_dir = '/Users/pnorton/Projects/National_Hydrology_Model'

byHRU_dir = f'{base_dir}/calibrations/NHMv10/byHRU_test1/ALL'
paramdb_dir = f'{base_dir}/calibrations/NHMv10/byHRU_test1/paramdb_v10_daymet_CONUS'


# %%
pdb_orig = ParamDb(paramdb_dir, verbose=True)

# %% [markdown]
# Example of reading a single-HRU parameter file from a calibration run

# %%
hru_param_file = f'{byHRU_dir}/ALL.SCE_HRU18325'
pfile_hru = ParameterFile(hru_param_file, verify=False)

# %%
# List parameters that were declared more than once in the parameter file
print(pfile_hru.updated_parameters)

# %%
# Show total number of updated/duplicated parameters
len(pfile_hru.updated_parameters)

# %%
# Check consistency
pfile_hru.parameters.check()

# %%
pfile_hru.parameters['snowpack_init'].modified

# %%
pfile_hru.parameters['soil2gw_max'].modified

# %%
filelist = glob.glob(f'{byHRU_dir}/ALL.SCE_HRU*')

print(f'Processing {len(filelist)} HRUs')
for cfile in sorted(filelist):
    chru = int(os.path.basename(cfile).split('HRU')[1])
    
    # Zero-based chru number for indexing
    chru0 = chru - 1
    
    sys.stdout.write(f'\rUpdating HRU {chru}: {cfile}    ')

    # Setting verify to False speeds up processing because the master parameter list will not
    # be loaded for each HRU parameter file.
    pfile_hru = ParameterFile(cfile, verify=False)
    
    for cparam in pfile_hru.updated_parameters:
        # Arrays are 0-based in python
        if pfile_hru.parameters[cparam].modified:
            if cparam == 'snarea_curve':
                pass
                # Special handling for snarea_curve
#                 snow_index = pfile_orig.parameters['hru_deplcrv'].data[chru0]
#                 pfile_orig.parameters[cparam].data.reshape((-1, 11))[snow_index-1, :] = pfile_hru.parameters[cparam].data
            else:
                pdb_orig.parameters.update_element(cparam, chru, pfile_hru.parameters[cparam].data[0])
#                 pdb_orig.parameters[cparam].update_element(chru0, pfile_hru.parameters[cparam].data[0])
#                 pdb_orig.parameters[cparam].data[chru0] = pfile_hru.parameters[cparam].data[0]

# %%
for cparam in pfile_hru.updated_parameters:
    print(f'{cparam}: {pdb_orig.parameters[cparam].modified}')

# %%
pdb_orig.parameters['jh_coef'].data[18324]

# %%
pfile_hru.parameters['snowpack_init'].modified

# %%
pdb_orig.dimensions.get('nhru').size

# %%
aa = Counter()

# %%
aa

# %%
aa.update(['joe', 'bob'])

# %%
aa

# %%
aa.update(['joe', 'steve'])

# %%
aa

# %%
pdb_orig.degenerate_parameters()

# %%
