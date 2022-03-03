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

from collections import OrderedDict


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

def get_param_default(search_name):
    dimensions = []

    for elem in xml_root.findall('parameter'):
        # print(elem.attrib.get('name'))
        name = elem.attrib.get('name')
        dtype = elem.find('type').text

        if search_name == name:
            def_value = elem.find('default').text
            
            if dtype == 'I':
                def_value = int(def_value)
            elif dtype in ['F', 'D']:
                def_value = float(def_value)
                
#             for cdim in elem.findall('.dimensions/dimension'):
#                 dimensions.append(cdim.attrib.get('name'))

    #         desc = elem.find('desc').text
    #         units = elem.find('units').text
    #         dimname = dim[1:]

            break
    return def_value    
    
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
# workdir_v1_parameters = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_CONUS'


datatype_conv = {'I': str_to_int, 'F': str_to_float, 'D': str_to_float, 'S': str_to_str}

# %% [markdown]
# ## Read the parameter information from parameters.xml

# %%
xml_fh = io.StringIO(pkgutil.get_data('pyPRMS', 'xml/parameters.xml').decode('utf-8'))
xml_tree = xmlET.parse(xml_fh)
xml_root = xml_tree.getroot()

# %% [markdown]
# ## Get mapping of local nhm_id to local index for TB and v1

# %%
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TB map of local nhm_id to local index
fhdl = open('{}/nhm_id.csv'.format(workdir_tb_parameters))

rawdata = fhdl.read().splitlines()
fhdl.close()
it = iter(rawdata)

# Dictionary with keys=nhm_id, values=idx0
tb_map = {}

# Skip header
next(it)

# Read rows
# For now the idx is left as 0-based
for idx, rec in enumerate(it):
    fld = rec.split(',')
    tb_map[int(fld[1])] = idx
    

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1 paramdb map of nhm_id to local index
fhdl = open('{}/nhm_id.csv'.format(workdir_v1_parameters))

rawdata = fhdl.read().splitlines()
fhdl.close()
it = iter(rawdata)

v1_map = {}

# Skip header
next(it)

# Read rows
# For now the idx is left as 0-based
for idx, rec in enumerate(it):
    fld = rec.split(',')
    v1_map[int(fld[1])] = idx

# %%

    

# %%
# HRU_mapping fields
# Field names: GFv11_id,Version,hru_id_nat,nhm_id,hru_segment_GFv11
fhdl = open('{}/HRU_mapping.csv'.format(workdir_tb_ids))

rawdata = fhdl.read().splitlines()
fhdl.close()
it = iter(rawdata)

hru_mapping = {}

# Skip header
next(it)

# Read rows
# For now the idx is left as 0-based
for rec in it:
    fld = rec.split(',')
    gfv11_id = int(fld[0])
    version = float(fld[1])
    v1_nhm_id = int(fld[2])
    tb_nhm_id = int(fld[3])
    
    if version < 1.1:
        v1_idx = v1_map[v1_nhm_id]
        tb_idx = -1
    else:
        try:
            if gfv11_id == 18868:
                # This is marked v1.1 only because the HRU was resized?
                print('ERROR: v:{} tb_nhm_id={}, v1_nhm_id={}, gfv11_id={} <<forced add entry>>'.format(version, tb_nhm_id, v1_nhm_id, gfv11_id))
                v1_idx = v1_map[v1_nhm_id]
                tb_idx = -1
            else:
                v1_idx = -1
                tb_idx = tb_map[tb_nhm_id]
        except KeyError:
            print('ERROR: v:{} tb_nhm_id={}, v1_nhm_id={}, gfv11_id={}'.format(version, tb_nhm_id, v1_nhm_id, gfv11_id))
            tb_idx = -2

#     if gfv11_id in v1_map:
#         v1_idx = v1_map[gfv11_id]
#     else:
#         v1_idx = -1
        
#     if gfv11_id in tb_map:
#         v11_idx = tb_map[gfv11_id]
#     else:
#         v11_idx = -1
       
    # List names are: version, v1_nhm_id, tb_nhm_id, hru_segment_GFv11, v1_idx, v11_idx
    hru_mapping[gfv11_id] = [version, v1_nhm_id, tb_nhm_id, int(fld[4]), v1_idx, tb_idx]
    

# %% [markdown]
# ## Get mapping of local nhm_seg to local index for TB and v1

# %%
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TB map of nhm_seg to local index
fhdl = open('{}/nhm_seg.csv'.format(workdir_tb_parameters))

rawdata = fhdl.read().splitlines()
fhdl.close()
it = iter(rawdata)

# Dictionary with keys=nhm_seg, values=idx0
tb_segid_map = {}

# Skip header
next(it)

# Read rows
# For now the idx is left as 0-based
for idx, rec in enumerate(it):
    fld = rec.split(',')
    tb_segid_map[int(fld[1])] = idx
    
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version 1 paramdb map of nhm_id to local index
fhdl = open('{}/nhm_seg.csv'.format(workdir_v1_parameters))

rawdata = fhdl.read().splitlines()
fhdl.close()
it = iter(rawdata)

v1_segid_map = {}

# Skip header
next(it)

# Read rows
# For now the idx is left as 0-based
for idx, rec in enumerate(it):
    fld = rec.split(',')
    v1_segid_map[int(fld[1])] = idx

# %%
# SEGMENT mapping fields
# Fields in SEG_mapping.csv: GFv11_id, Version, seg_id_nhm, tosegment_GFv11
fhdl = open('{}/SEG_mapping.csv'.format(workdir_tb_ids))

rawdata = fhdl.read().splitlines()
fhdl.close()
it = iter(rawdata)

seg_mapping = {}
seg_old_to_new = {}

# Skip header
next(it)

# Read rows
# For now the idx is left as 0-based
for rec in it:
    fld = rec.split(',')
    gfv11_id = int(fld[0])
    old_seg_id_nhm = int(fld[2])
    version = float(fld[1])
    
    if version < 1.1:
        try:
            v1_idx = v1_segid_map[old_seg_id_nhm]
            tb_idx = -1
        except KeyError:
            print(f'old_seg_id_nhm: {old_seg_id_nhm}, is missing in nhm_seg.csv for CONUS (v. 1)')            
    else:
        try:
            v1_idx = -1
            tb_idx = tb_segid_map[old_seg_id_nhm]
        except KeyError:
            print(f'old_seg_id_nhm: {old_seg_id_nhm}, is missing in nhm_seg.csv for the transboundary')
           
    # List names are: version, old_seg_id_nhm, tosegment_GFv11, v1_idx, tb_idx
    seg_mapping[gfv11_id] = [version, old_seg_id_nhm, int(fld[3]), v1_idx, tb_idx]
    seg_old_to_new[old_seg_id_nhm] = gfv11_id

# %%

# %% [markdown]
# ## Merge the nhru-dimensioned parameters

# %%
# Use the transboundary parameter file names to drive the processing

tb_files = []

file_it = glob.glob('{}/*.csv'.format(workdir_tb_parameters))
for kk in file_it:
    tb_files.append(kk)
tb_files.sort()

nhru = max(list(hru_mapping.keys()))
nsegment = max(list(seg_mapping.keys()))
    
for ff in tb_files:
    filename = os.path.basename(ff)
    cname = os.path.basename(os.path.splitext(ff)[0])
    
    if cname in ['nhm_to_GFv11_HRU', 'nhm_to_GFv11_SEG', 'seg_id_nhm', 'obsout_segment',
                 'poi_gage_id', 'poi_gage_segment', 'poi_type']:
        continue
        
    print(cname)

    tb_paramfile = '{}/{}'.format(workdir_tb_parameters, filename)
    v1_paramfile = '{}/{}'.format(workdir_v1_parameters, filename)

    # Get the dimensions for the parameter
    dimensions, dtype = get_dimensions(cname)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Read the TB parameter file
    tb_paramdata = read_parameter(tb_paramfile)
    tb_paramdata = datatype_conv[dtype](tb_paramdata)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Read the v1 parameter file
    v1_paramdata = read_parameter(v1_paramfile)
    v1_paramdata = datatype_conv[dtype](v1_paramdata)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Now create the merged paramdb file
    if len(dimensions) == 1:
        # Open a file to write to
        crap = open('/Users/pnorton/tmp/tmp_paramdb/{}.csv'.format(cname), 'w')
        crap.write('$id,{}\n'.format(cname))

        if 'nhru' in dimensions:
            # Loop through the GFv11_id's in order and write a new v11 parameter file
            
            if cname == 'nhm_id':
                # The nhm_id is replaced with GFv11_id
                print('    - replaced with GFv11_id')
                for xx in range(0, nhru):
                    idx1 = xx + 1
                    crap.write('{},{}\n'.format(idx1, idx1))
            elif cname in ['hru_segment', 'hru_segment_nhm']:
                # hru_segment_nhm is replaced with hru_segment_GFv11
                print('    - replaced with hru_segment_GFv11')
                for xx in range(0, nhru):
                    idx1 = xx + 1
                    crap.write('{},{}\n'.format(idx1, hru_mapping[idx1][3]))
            elif cname == 'hru_deplcrv':
                # 2020-02-28 NOTE: This is just a placeholder to get some snow depletion
                #                  curves in the v1.1 database. They will get replaced
                print(f'    - using default value = 1 for CONUS')
#                 for xx in range(0, nhru):
#                     idx1 = xx + 1
#                     crap.write('{},{}\n'.format(idx1, 1))
                
                for xx in range(0, nhru):
                    idx1 = xx + 1
                    curr_rec = hru_mapping[idx1]
                    if curr_rec[0] < 1.1:
                        crap.write(f'{idx1},1\n')
                    else:
                        crap.write(f'{idx1},{tb_paramdata[curr_rec[5]]}\n')
#                         crap.write('{},{}\n'.format(idx1, tb_paramdata[curr_rec[5]]))
            elif cname == 'potet_sublim':
                # 2020-02-27 This has the wrong default in TB
                def_value = get_param_default(cname)
                print(f'    - using default value = {def_value}')
                for xx in range(0, nhru):
                    idx1 = xx + 1
                    crap.write('{},{}\n'.format(idx1, def_value))                
            elif cname in ['dprst_frac_init', 'gwstor_init', 'snowpack_init',
                           'ssstor_init_frac', 'soil_moist_init_frac', 'soil_rechr_init_frac']:
                # Use default value for init parameters
                def_value = get_param_default(cname)
                print(f'    - using default value = {def_value}')
                for xx in range(0, nhru):
                    idx1 = xx + 1
                    crap.write('{},{}\n'.format(idx1, def_value))
            else:
                # Names for hru_mapping entries are: version, v1_nhm_id, v11_nhm_id, hru_segment_GFv11, v1_idx, v11_idx
                for xx in range(0, nhru):
                    idx1 = xx + 1
                    curr_rec = hru_mapping[idx1]
                    if curr_rec[0] < 1.1:
                        crap.write('{},{}\n'.format(idx1, v1_paramdata[curr_rec[4]]))
                    else:
                        crap.write('{},{}\n'.format(idx1, tb_paramdata[curr_rec[5]]))
        elif 'one' in dimensions:
            # Just use the v1 parameter value
            crap.write('1,{}\n'.format(v1_paramdata[0]))
        elif 'ndeplval' in dimensions:
            # Use the depletions curves from TB params only
            for ii, xx in enumerate(tb_paramdata):
                crap.write(f'{ii+1},{xx}\n')
        elif 'npoigages' in dimensions:
            # POI processing is handled below
            continue
        elif 'nsegment' in dimensions:
            if cname == 'nhm_seg':
                # The nhm_seg is replaced with GFv11_id from seg_mapping
                print('    - replace with GFv11_id')
                for xx in range(0, nsegment):
                    idx1 = xx + 1
                    crap.write('{},{}\n'.format(idx1, idx1))
            elif cname in ['tosegment', 'tosegment_nhm']:
                # tosegment_nhm is replaced with tosegment_GFv11
                print('    - replaced with tosegment_GFv11')
                for xx in range(0, nsegment):
                    idx1 = xx + 1
                    crap.write('{},{}\n'.format(idx1, seg_mapping[idx1][2]))
            elif cname in ['segment_flow_init']:
                # Use default value for init parameters
                def_value = get_param_default(cname)
                print(f'    - using default value = {def_value}')
                for xx in range(0, nsegment):
                    idx1 = xx + 1
                    crap.write('{},{}\n'.format(idx1, def_value))
            else:
                for xx in range(0, nsegment):
                    # [version, old_seg_id_nhm, tosegment_GFv11, v1_idx, v11_idx]
                    idx1 = xx + 1
                    curr_rec = seg_mapping[idx1]
                    try:
                        if curr_rec[0] < 1.1:
                            crap.write('{},{}\n'.format(idx1, v1_paramdata[curr_rec[3]]))
                        else:
                            crap.write('{},{}\n'.format(idx1, tb_paramdata[curr_rec[4]]))
                    except KeyError:
                        print('ERROR: seg: version={}, old_seg_id={}, tosegment={}, v1_idx={}, tb_idx={}'.format(curr_rec[0], curr_rec[1], curr_rec[2], curr_rec[3], curr_rec[4]))
        crap.close()
        
    elif len(dimensions) == 2:
        # Convert tb and v1 parameter data to Fortran-ordered arrays
        
        # Open a file to write to
        crap = open('/Users/pnorton/tmp/tmp_paramdb/{}.csv'.format(cname), 'w')
        crap.write('$id,{}\n'.format(cname))
        
        if 'nmonths' in dimensions:
            tb_np = np.array(tb_paramdata).reshape((-1, 12,), order='F')
            v1_np = np.array(v1_paramdata).reshape((-1, 12,), order='F')
            
            if 'nhru' in dimensions:
                # Build a merged output array
                out_np = np.empty((nhru, 12), order='F')

                for xx in range(0, nhru):
                    idx1 = xx + 1
                    curr_rec = hru_mapping[idx1]

                    if curr_rec[0] < 1.1:
                        out_np[xx, :] = v1_np[curr_rec[4], :]
                    else:
                        out_np[xx, :] = tb_np[curr_rec[5], :]

                for ii, dd in enumerate(out_np.ravel(order='F').tolist()):
                    crap.write('{},{}\n'.format(ii+1, dd))
            elif 'nsegment' in dimensions:
                # Build a merged output array
                out_np = np.empty((nsegment, 12), order='F')

                for xx in range(0, nsegment):
                    idx1 = xx + 1
                    curr_rec = seg_mapping[idx1]

                    if curr_rec[0] < 1.1:
                        out_np[xx, :] = v1_np[curr_rec[3], :]
                    else:
                        out_np[xx, :] = tb_np[curr_rec[4], :]

                for ii, dd in enumerate(out_np.ravel(order='F').tolist()):
                    crap.write('{},{}\n'.format(ii+1, dd))
        else:
            print('ERROR: Unable to handle dimension, {}'.format(dimensions[1]))
            
        crap.close()

# %%
# Check for POIs common between TB and v1 

tb_files = []

file_it = glob.glob('{}/*.csv'.format(workdir_tb_parameters))
for kk in file_it:
    tb_files.append(kk)
tb_files.sort()

for ff in tb_files:
    filename = os.path.basename(ff)
    
    cname = os.path.basename(os.path.splitext(ff)[0])
    if cname != 'poi_gage_id':
        continue

    # Get the dimensions for the parameter
    dimensions, dtype = get_dimensions(cname)
    
    tb_paramfile = '{}/{}'.format(workdir_tb_parameters, filename)
    v1_paramfile = '{}/{}'.format(workdir_v1_parameters, filename)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Read the TB parameter file
    tb_paramdata = read_parameter(tb_paramfile)
    tb_paramdata = datatype_conv[dtype](tb_paramdata)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Read the v1 parameter file
    v1_paramdata = read_parameter(v1_paramfile)
    v1_paramdata = datatype_conv[dtype](v1_paramdata)

    common = set(v1_paramdata) & set(tb_paramdata)

print(f'Number of POIs common between TB and v1: {len(common)}')

# %%

# %%
# Process the POIs
poi_tmp = {}
repl_count = 0

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in the version 1 POIs
print('Reading v1 POIs')
poi_type = read_parameter(f'{workdir_v1_parameters}/poi_type.csv')
poi_gage_id = read_parameter(f'{workdir_v1_parameters}/poi_gage_id.csv')
poi_gage_segment = read_parameter(f'{workdir_v1_parameters}/poi_gage_segment.csv')

for ii, tt, ss in zip(poi_gage_id, poi_type, poi_gage_segment):
    if ii in poi_tmp:
#         print(f'  {ii} already exists ({poi_tmp[ii][0]}, {poi_tmp[ii][1]}), updating with new values ({ss}, {tt})')
        repl_count += 1
        
    poi_tmp[ii] = [ss, tt]

print(f'  Number of POIs: {len(poi_tmp)}')
print(f'  Number of updated POIs: {repl_count}')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in the transboundary POIs
print('Reading transboundary POIs')
poi_type = read_parameter(f'{workdir_tb_parameters}/poi_type.csv')
poi_gage_id = read_parameter(f'{workdir_tb_parameters}/poi_gage_id.csv')
poi_gage_segment = read_parameter(f'{workdir_tb_parameters}/poi_gage_segment.csv')

for ii, tt, ss in zip(poi_gage_id, poi_type, poi_gage_segment):
    if ii in poi_tmp:
#         print(f'  {ii} already exists ({poi_tmp[ii][0]}, {poi_tmp[ii][1]}), updating with new values ({ss}, {tt})')
        repl_count += 1
        
    poi_tmp[ii] = [ss, tt]

print(f'  Number of POIs: {len(poi_tmp)}')
print(f'  Number of updated POIs: {repl_count}')

# Writing the POI parameters
print('Writing the POI parameter files')
poi_type_hdl = open('/Users/pnorton/tmp/tmp_paramdb/poi_type.csv', 'w')
poi_type_hdl.write('$id,poi_type\n')

poi_gage_id_hdl = open('/Users/pnorton/tmp/tmp_paramdb/poi_gage_id.csv', 'w')
poi_gage_id_hdl.write('$id,poi_gage_id\n')

poi_gage_segment_hdl = open('/Users/pnorton/tmp/tmp_paramdb/poi_gage_segment.csv', 'w')
poi_gage_segment_hdl.write('%id,poi_gage_segment\n')

idx1 = 0
for kk, vv in poi_tmp.items():
    idx1 += 1
    poi_type_hdl.write(f'{idx1},{vv[1]}\n')
    poi_gage_id_hdl.write(f'{idx1},{kk}\n')
    poi_gage_segment_hdl.write(f'{idx1},{vv[0]}\n')
    

poi_type_hdl.close()
poi_gage_id_hdl.close()
poi_gage_segment_hdl.close()

# %%
