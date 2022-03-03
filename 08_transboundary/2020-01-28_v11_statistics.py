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
import numpy as np
import os

import io
import pkgutil
import xml.dom.minidom as minidom
import xml.etree.ElementTree as xmlET


# %%
def str_to_float(data):
    # Convert provide list of data to float
    try:
        return [float(vv) for vv in data]
    except ValueError as ve:
        print(ve)

def str_to_int(data):
    # Convert list of data to integer
    try:
        return [int(vv) for vv in data]
    except ValueError as ve:
        # Perhaps it's a float, try converting to float and then integer
        try:
            tmp = [float(vv) for vv in data]
            return [int(vv) for vv in tmp]
        except ValueError as ve:
            print(ve)

def str_to_str(data):
    # nop for list of strings
    # 2019-05-22 PAN: For python 3 force string type to byte
    #                 otherwise they are treated as unicode
    #                 which breaks the write_netcdf() routine.
    # 2019-06-26 PAN: Removed the encode because it broken writing the ASCII
    #                 parameter files. Instead the conversion to ascii is
    #                 handled in the write_netcdf routine of ParameterSet
    # data = [dd.encode() for dd in data]
    return data


def get_dimensions(search_name):
#     search_name = 'radmax'
    dimensions = []

    for elem in xml_root.findall('parameter'):
        # print(elem.attrib.get('name'))
        name = elem.attrib.get('name')
        dtype = elem.find('type').text

        if search_name == name:
            for cdim in elem.findall('.dimensions/dimension'):
                dimensions.append(cdim.attrib.get('name'))

    #         desc = elem.find('desc').text
    #         units = elem.find('units').text
    #         dimname = dim[1:]

            break
    return dimensions, dtype


def read_parameter(filename):
    # Open and read the parameter file
    fhdl = open(filename, 'r')

    rawdata = fhdl.read().splitlines()
    fhdl.close()
    it = iter(rawdata)

    paramdata = []

    # Skip header
    next(it)

    # Read rows
    for rec in it:
        fld = rec.split(',')
        paramdata.append(fld[1])
    
    return paramdata


# %%
workdir_tb_ids = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/2019-11-25_work'
workdir_tb_parameters = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/tbparamdb'

workdir_v1_parameters = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/v1_paramdb'


datatype_conv = {'I': str_to_int, 'F': str_to_float, 'D': str_to_float, 'S': str_to_str}

# %% [markdown]
# ## Read the parameter information from parameters.xml

# %%
xml_fh = io.StringIO(pkgutil.get_data('pyPRMS', 'xml/parameters.xml').decode('utf-8'))
xml_tree = xmlET.parse(xml_fh)
xml_root = xml_tree.getroot()

# %% [markdown]
# ## Output basic statistics for each TB or v1 file

# %%
# Use the transboundary parameter file names to drive the processing

tb_files = []

file_it = glob.glob('{}/*.csv'.format(workdir_tb_parameters))
for kk in file_it:
    tb_files.append(kk)
tb_files.sort()

output_type = 'v1'  # one of 'v1', 'TB'

print('Parameter,tb_min,tb_max,tb_mean,tb_median,v1_min,v1_max,v1_mean,v1_median')
for ff in tb_files:
    filename = os.path.basename(ff)
    
    cname = os.path.basename(os.path.splitext(ff)[0])
    if cname in ['nhm_to_GFv11_HRU', 'nhm_to_GFv11_SEG', 'seg_id_nhm', 'obsout_segment',
                 'seg_depth', 'seg_length', 'seg_slope', 'mann_n', 'poi_gage_id']:
        continue

    tb_paramfile = '{}/{}'.format(workdir_tb_parameters, filename)
    v1_paramfile = '{}/{}'.format(workdir_v1_parameters, filename)

    # Get the dimensions for the parameter
    dimensions, dtype = get_dimensions(cname)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Read the TB parameter file
    tb_paramdata = read_parameter(tb_paramfile)
    tb_paramdata = datatype_conv[dtype](tb_paramdata)

    tb_data_np = np.array(tb_paramdata)
    tb_data_min = np.min(tb_data_np)
    tb_data_max = np.max(tb_data_np)
    tb_data_mean = np.mean(tb_data_np)
    tb_data_median = np.median(tb_data_np)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Read the v1 parameter file
    v1_paramdata = read_parameter(v1_paramfile)
    v1_paramdata = datatype_conv[dtype](v1_paramdata)

    v1_data_np = np.array(v1_paramdata)
    v1_data_min = np.min(v1_data_np)
    v1_data_max = np.max(v1_data_np)
    v1_data_mean = np.mean(v1_data_np)
    v1_data_median = np.median(v1_data_np)
    
    print(f'{cname},{tb_data_min},{tb_data_max},{tb_data_mean},{tb_data_median},{v1_data_min},{v1_data_max},{v1_data_mean},{v1_data_median}')
    
#     if data_min == data_max:
#         print(f'{cname} | {data_min} |')
#     else:
#         print(f'{cname} | {data_min} | {data_max} | {np.mean(data_np)} | {np.median(data_np)}')

# %%

# %%

# %%

# %%
