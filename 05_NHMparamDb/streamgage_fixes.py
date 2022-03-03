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
from future.utils import iteritems, itervalues

from itertools import izip
from collections import OrderedDict
import pandas as pd
import numpy as np

import pyPRMS.NhmParamDb as nhm

REGIONS = ['r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09',
           'r10L', 'r10U', 'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18']

srcdir = '/Users/pnorton/Projects/National_Hydrology_Model/tmp'
nhm_srcdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb'

# %%
gage_chg_file = '{}/NewGageSeg.txt'.format(srcdir)

col_names = ['region', 'site_no', 'segment']
col_types = [np.str_, np.str_, np.int_]
cols = dict(zip(col_names, col_types))

df_gage_chg = pd.read_csv(gage_chg_file, sep=' ', dtype=cols)
df_gage_chg.head()

# %%
# df_gage_chg[df_gage_chg.loc[:,'region'] == 'r01']['site_no'].tolist()
df_gage_chg[df_gage_chg.loc[:,'region'] == 'r01']

# %%

# %%
for rr in REGIONS:
# for rr in ['r01']:
    # Get subset of segment changes for current region
    df2 = df_gage_chg[df_gage_chg.loc[:,'region'] == rr]
    gage_changes = dict(zip(df2.site_no,df2.segment))
    
    gages_dict = OrderedDict()
    
    # First read the poi_gage_id file
    fhdl = open('{}/poi_gage_id/{}/poi_gage_id.csv'.format(nhm_srcdir, rr), 'r')
    rawdata_gages = fhdl.read().splitlines()
    fhdl.close()
    it_gages = iter(rawdata_gages)
    next(it_gages)
    
    # Next read the poi_gage_segment file
    fhdl = open('{}/poi_gage_segment/{}/poi_gage_segment.csv'.format(nhm_srcdir, rr), 'r')
    rawdata_segs = fhdl.read().splitlines()
    fhdl.close()
    it_segs = iter(rawdata_segs)
    next(it_segs)   
    
    for rec_gages, rec_segs in izip(it_gages, it_segs):
        idx_gages, val_gages = rec_gages.split(',')
        idx_segs, val_segs = rec_segs.split(',')
        
        gages_dict[val_gages] = val_segs
        
    # Now update sites that have new segments
    for site, newseg in iteritems(gage_changes):
        try:
            if gages_dict[site] != newseg:
                print('Changing segment for site {} in region {} from {} to {}'.format(site, rr, gages_dict[site], newseg))
                gages_dict[site] = newseg
        except KeyError:
            print('WARNING: Streamgage {} does not exist in nhmparamDb'.format(site))

    # Lastly write the changed files for the nhmParamDb
#     outhdl = open('{}/20181026_nhm_updates/{}_poi_gage_segment.csv'.format(srcdir, rr), 'w')
    outhdl = open('{}/poi_gage_segment/{}/poi_gage_segment.csv'.format(nhm_srcdir, rr), 'w')
    outhdl.write('$id,poi_gage_segment\n')
    
    idx = 1
    for val in itervalues(gages_dict):
        outhdl.write('{},{}\n'.format(idx, val))
        idx += 1
        
    outhdl.close()
                  

# %%
# Update poi_type for each streamgage

for rr in REGIONS:
# for rr in ['r01']:
    gages_dict = OrderedDict()
    
    # First read the poi_gage_id file
    fhdl = open('{}/poi_gage_id/{}/poi_gage_id.csv'.format(nhm_srcdir, rr), 'r')
    rawdata_gages = fhdl.read().splitlines()
    fhdl.close()
    it_gages = iter(rawdata_gages)
    next(it_gages)
    
    # Next read the poi_gage_segment file
    fhdl = open('{}/poi_type/{}/poi_type.csv'.format(nhm_srcdir, rr), 'r')
    rawdata_segs = fhdl.read().splitlines()
    fhdl.close()
    it_segs = iter(rawdata_segs)
    next(it_segs)   
    
    for rec_gages, rec_segs in izip(it_gages, it_segs):
        idx_gages, val_gages = rec_gages.split(',')
        idx_segs, val_segs = rec_segs.split(',')
        
        gages_dict[val_gages] = val_segs
        
    # Now update sites that have new poi_type
    fhdl = open('{}/Gage_poiType/{}_gage_poiType'.format(srcdir, rr), 'r')
    rawdata = fhdl.read().splitlines()
    fhdl.close()
    it_type = iter(rawdata)
    next(it_type)
    
    poitype_dict = OrderedDict()
    for rec_chgs in it_type:
        site, newtype = rec_chgs.split(',')
        poitype_dict[site] = newtype
        
    for site, newtype in iteritems(poitype_dict):
        try:
            if gages_dict[site] != newtype:
                print('Changing poi_type for site {} in region {} from {} to {}'.format(site, rr, gages_dict[site], newtype))
                gages_dict[site] = newtype
        except KeyError:
            print('WARNING: Streamgage {} does not exist in nhmparamDb'.format(site))
            
    # Lastly write the changed files for the nhmParamDb
#     outhdl = open('{}/20181026_nhm_updates/{}_poi_type.csv'.format(srcdir, rr), 'w')
    outhdl = open('{}/poi_type/{}/poi_type.csv'.format(nhm_srcdir, rr), 'w')
    outhdl.write('$id,poi_type\n')
    
    idx = 1
    for val in itervalues(gages_dict):
        outhdl.write('{},{}\n'.format(idx, val))
        idx += 1
        
    outhdl.close()
#     print(gages_dict)

# %%
