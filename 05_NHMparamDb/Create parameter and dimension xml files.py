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

from pyPRMS.ParameterSet import ParameterSet


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
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_CONUS'

# %%
# Get the list of parameter files to read
file_list = []

file_it = glob.glob('{}/*.csv'.format(workdir))
for kk in file_it:
    file_list.append(kk)
    
file_list.sort()

# %%
new_ps = ParameterSet(verbose=True)

# =======================================================
# Set global nhru dimension based on nhm_id parameter
it = _data_it('{}/nhm_id.csv'.format(workdir))
next(it)

tmp_data = []

# Read the parameter values
for rec in it:
    idx, val = rec.split(',')
    tmp_data.append(val)
    
new_ps.dimensions.add('nhru', size=len(tmp_data))
print('nhm_id: {}'.format(len(tmp_data)))

# =======================================================
# Set global nsegment dimension based on nhm_seg parameter
it = _data_it('{}/nhm_seg.csv'.format(workdir))
next(it)

tmp_data = []

# Read the parameter values
for rec in it:
    idx, val = rec.split(',')
    tmp_data.append(val)
    
new_ps.dimensions.add('nsegment', size=len(tmp_data))
print('nhm_seg: {}'.format(len(tmp_data)))

# =======================================================
# Set global npoigages dimension based on poi_gage_id parameter
it = _data_it('{}/poi_gage_id.csv'.format(workdir))
next(it)

tmp_data = []

# Read the parameter values
for rec in it:
    idx, val = rec.split(',')
    tmp_data.append(val)
    
new_ps.dimensions.add('npoigages', size=len(tmp_data))
print('poi_gage_id: {}'.format(len(tmp_data)))

# =============================================================
# Set global ndeplval dimension based on snarea_curve parameter
it = _data_it('{}/snarea_curve.csv'.format(workdir))
next(it)

tmp_data = []

# Read the parameter values
for rec in it:
    idx, val = rec.split(',')
    tmp_data.append(val)
    
new_ps.dimensions.add('ndeplval', size=len(tmp_data))
print('ndeplval: {}'.format(len(tmp_data)))

global_dims = new_ps.dimensions

# =======================================================
# Load the parameters
for ff in file_list:
    cname = os.path.basename(os.path.splitext(ff)[0])
    print(cname)
    
    # Read the parameter data
    tmp_data = []

    # Read parameter information
    try:
        it = _data_it(ff)
        next(it)  # Skip the header row
    except IOError:
        print('Skipping parameter: {}. File does not exist.'.format(ff))
        continue

    if new_ps.master_parameters is not None:
        new_ps.parameters.add(cname, info=new_ps.master_parameters[cname])
    else:
        new_ps.parameters.add(cname)

    # Add the dimension sizes
    master_info = new_ps.master_parameters[cname]
    for dd in master_info.dimensions.values():
        if dd.name in ['nhru', 'nsegment']:
            new_ps.parameters[cname].dimensions.add(dd.name, size=global_dims[dd.name].size)
        elif dd.name == 'npoigages':
            new_ps.parameters[cname].dimensions.add(dd.name, size=global_dims[dd.name].size)
        elif dd.name == 'ndeplval':
            new_ps.parameters[cname].dimensions.add(dd.name, size=global_dims[dd.name].size)
        else:
            new_ps.parameters[cname].dimensions.add(dd.name)
            
    # Read the parameter values
    for rec in it:
        idx, val = rec.split(',')
        tmp_data.append(val)

    new_ps.parameters[cname].data = tmp_data

# %%
print(new_ps.dimensions)

# %%
# Loop through all the parameter dimensions; add any missing global dimensions; 
# check that the sizes of dimensions don't change compared to the global dimension
for pp in new_ps.parameters.values():
    for dd in pp.dimensions.values():
        if new_ps.dimensions.exists(dd.name):
            if new_ps.dimensions[dd.name].size != dd.size:
                print("{}: Sizes don't match ({} != {})".format(dd.name, new_ps.dimensions[dd.name].size, dd.size))
        else:
            new_ps.dimensions.add(dd.name, dd.size)

print(new_ps.dimensions)

# %%
# Check the parameters for validity
new_ps.parameters.check()

# %%
# Write the parameteter and dimension XML files for the parameter database
new_ps.write_parameters_xml(workdir)
new_ps.write_dimensions_xml(workdir)

# %%
