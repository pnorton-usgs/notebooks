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
from collections import OrderedDict
import numpy as np
import pandas as pd

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS'
gf_poi_filename = f'{workdir}/poi_agency.csv'

paramdb_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS'


# %%
def read_csv(filename):
    # Read the 2nd column of a paramdb file into a list
    fhdl = open(filename)
    rawdata = fhdl.read().splitlines()
    fhdl.close()
    it = iter(rawdata)

    data = []

    # Skip header
    next(it)

    # Read the parameter values
    for rec in it:
        data.append(rec.split(',')[1])
    return data


# %%
def read_csv_pairs(filename1, filename2):
    # Given two paramdb files this will return an ordered dictionary of the 2nd columns of each file
    # where the first file is the key, second file is the value
    fhdl = open(filename1)
    rawdata = fhdl.read().splitlines()
    fhdl.close()
    it = iter(rawdata)
    
    fhdl2 = open(filename2)
    rawdata2 = fhdl2.read().splitlines()
    fhdl2.close()
    it2 = iter(rawdata2)
    
    data = OrderedDict()
    next(it)
    next(it2)
    
    for lh,rh in zip(it, it2):
        data[lh.split(',')[1]] = int(rh.split(',')[1])
    return data


# %%
# Get the POI station IDs from the paramdb
pdb_poi_id = read_csv(f'{paramdb_dir}/poi_gage_id.csv')

print(f'poi_ids: {len(pdb_poi_id)}')
print(f'unique: {len(set(pdb_poi_id))}')

# %%
# GNIS_Name	Type_Gage	Type_Ref	Gage_Source	poi_segment_v1_1
col_names = ['GNIS_Name', 'Type_Gage', 'Type_Ref', 'Gage_Source', 'poi_segment_v1_1']
col_types = [np.str_, np.str_, np.str_, np.str_, np.int]

cols = dict(zip(col_names, col_types))

df_poi = pd.read_csv(gf_poi_filename, sep='\t', dtype=cols, index_col=1)

# %%
df_poi.loc[:, ('poi_segment_v1_1')].to_dict()

# %%
# df_poi.loc[:, ('Gage_Source')].tolist()
df_poi.index.tolist()

# %%
fhdl = open(gf_poi_filename, 'r')

all_gages = []
ec_gages = []

for row in fhdl:
    flds = row.strip().split(',')
    gage_id = flds[1]
    gage_src = flds[3]
    
    if gage_id in pdb_poi_id:
        if gage_src != 'EC':
            if gage_src == '0':
                print(f'{gage_id} in PARAMDB but has no gage_src')
            try:
                gage_id_int = int(gage_id)
#                 print(f'{gage_id} in PARAMDB ({gage_src})')
                
                if len(gage_id) > 15:
                    print(f'{gage_id} is not a USGS gage')
            except ValueError:
                print(f'{gage_id} incorrectly sourced to {gage_src}')

                if gage_id != 'Type_Gage':
                    ec_gages.append(gage_id)
        elif gage_src == 'EC':
            if len(gage_id) > 7:
                print(f'{gage_id} incorrectly sourced to {gage_src}')
            else:
                ec_gages.append(gage_id)
        all_gages.append(gage_id)
    else:
        print(f'{gage_id} not included in paramdb ({gage_src})')

# %%
missing = set(pdb_poi_id) - set(all_gages)
len(missing)

# %%
list(missing)

# %%
set(all_gages) - set(pdb_poi_id)

# %%
# Get the POI station IDs from the paramdb
pdb_pois = read_csv_pairs(f'{paramdb_dir}/poi_gage_id.csv', f'{paramdb_dir}/poi_gage_segment.csv')

# %%

# %%
# Build ordered dictionary of geospatial fabric POIs
# key->gage_id, value->gage_seg

fhdl = open(gf_poi_filename, 'r')
rawdata = fhdl.read().splitlines()
fhdl.close()
it = iter(rawdata)
next(it)

gf_pois = OrderedDict()

for row in it:
    flds = row.strip().split(',')
    gage_id = flds[1]
    gage_seg = int(flds[4])
    
    gf_pois[gage_id] = gage_seg


# %%

# %%
# Check POIs for incorrect gage segment; compares paramdb gage segment to GF gage segment
mis_cnt = 0
mis_v11_cnt = 0

for nhm_poi, nhm_seg in pdb_pois.items():
    if nhm_poi not in gf_pois:
#         print(f'NHM POI: {nhm_poi} not in GFv11')
        mis_v11_cnt += 1
    else:
        if nhm_seg != gf_pois[nhm_poi]:
            print(f'{nhm_poi} v11:{nhm_seg}, GF:{gf_pois[nhm_poi]}, v10:{v1_pdb_pois.get(nhm_poi)}, tb:{tb_pdb_pois.get(nhm_poi)}')
            mis_cnt += 1

# %%
print(f'Number of incorrect POIs segments: {mis_cnt}')
print(f'Number of POIs missing from GF: {mis_v11_cnt}')

# %%

# %%
v1_pdb_src = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/v1_paramdb'
v1_pdb_pois = read_csv_pairs(f'{v1_pdb_src}/poi_gage_id.csv', f'{v1_pdb_src}/poi_gage_segment.csv')

# %%
tb_pdb_src = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/tbparamdb'
tb_pdb_pois = read_csv_pairs(f'{tb_pdb_src}/poi_gage_id.csv', f'{tb_pdb_src}/poi_gage_segment.csv')

# %%
# Lookup segment information for a single POI
poi = '01011000'

print(f'       GF v11 segment: {gf_pois[poi]}')
print(f'  Paramdb v11 segment: {pdb_pois[poi]}')
print(f'Transboundary segment: {tb_pdb_pois[poi]}')
print(f'  Paramdb v10 segment: {v1_pdb_pois[poi]}')


# %%
import datetime

# %%
ed = datetime.datetime(2020,7,20)
st = datetime.datetime(1911,1,26)

# %%
(ed - st) + datetime.timedelta(days=1)

# %%
aa  = 'Downloading observations for streamgage:                '
len(aa)

# %%

# %%

# %%

# %%
