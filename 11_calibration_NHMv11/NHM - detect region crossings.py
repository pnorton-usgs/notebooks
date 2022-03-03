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
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import pandas as pd
import numpy as np

from collections import OrderedDict

# %%
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/FY2021_projects/extraction_requests'

filename = f'{work_dir}/segs_with_huc2.csv'

# %%
# col_names = ['nsegment_v1_1','tosegment_v1_1','Version','HUC2','NAME']
# col_types = [int, int, float, str, str]
# cols = dict(zip(col_names, col_types))

# df = pd.read_csv(filename, sep=',', usecols=[0, 1, 3], dtype=cols)

# %% [markdown]
# ## Read segs_with_huc2.csv and generate mappings of 1) segment to huc2, 2) segment to tosegment

# %%
fhdl = open(filename, 'r')
rawdata = fhdl.read().splitlines()
fhdl.close()

seg_huc = OrderedDict()
seg_toseg = OrderedDict()

it = iter(rawdata)
next(it)

for row in it:
    fields = row.split(',')
    
    seg_huc[int(fields[0])] = fields[3]
    seg_toseg[int(fields[0])] = int(fields[1])


# %% [markdown]
# ## Write seg_outlets.csv
# This file contains all outlets in the NHM with HUC2 designation (HUC2 == '00' indicates segment leaving the NHM)

# %%
outfile = open(f'{work_dir}/seg_outlets.csv', 'w')
outfile.write('seg_id,toseg,from_huc,to_huc\n')

for seg, toseg in seg_toseg.items():
    if toseg == 0:
        print(f'{seg}, HUC: {seg_huc[seg]} is an outlet')
        outfile.write(f'{seg},{toseg},{seg_huc[seg]},00\n')
    elif seg_huc[seg] != seg_huc[toseg]:
        print(f'{seg}, HUC: {seg_huc[seg]} ===> {toseg}, HUC: {seg_huc[toseg]}')
        outfile.write(f'{seg},{toseg},{seg_huc[seg]},{seg_huc[toseg]}\n')
        
outfile.close()

# %% [markdown]
# ## Create mapping of region segment outlets and region cutoff segments

# %%
region_seg_outlets = {}
region_outlet_cutoffs = {}

for seg, toseg in seg_toseg.items():
    chuc = seg_huc[seg]

    if toseg == 0:
        region_seg_outlets.setdefault(chuc, []).append(seg)
        print(f'{seg}, HUC: {chuc} is an outlet')
    elif chuc != seg_huc[toseg]:
        region_seg_outlets.setdefault(chuc, []).append(seg)
        print(f'{seg}, HUC: {chuc} ===> {toseg}, HUC: {seg_huc[toseg]}')
        
        region_outlet_cutoffs.setdefault(seg_huc[toseg], []).append(seg)


# %%
region_outlet_cutoffs

# %% [markdown]
# ## Create mapping of HRU-to-HUC2 from non-routed_HRUs_by_HUC2.csv

# %%
hru_file = f'{work_dir}/non-routed_HRUs_by_HUC2.csv'

fhdl = open(hru_file, 'r')
rawdata = fhdl.read().splitlines()
fhdl.close()

it = iter(rawdata)
next(it)

nonrouted_hrus = {}
for row in it:
    fields = row.split(',')
    nonrouted_hrus.setdefault(fields[1], []).append(int(fields[0]))
    
# nonrouted_hrus

# %% [markdown]
# ## Write bandit config files for each NHM region

# %%
from Bandit import bandit_cfg as bc

bandit_work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/FY2021_projects/extraction_requests/bandit_configs'

base_bandit_config = f'{bandit_work_dir}/bandit_v11_gridmet.cfg'

for kk, vv in region_seg_outlets.items():
    # Read the default configuration file
    config = bc.Cfg(base_bandit_config)
    
    # Update the outlets in the basin.cfg file and write into the headwater directory
    config.update_value('outlets', vv)

    if kk in nonrouted_hrus:
        config.update_value('hru_noroute', nonrouted_hrus[kk])
    # if noroute_hrus_by_loc is not None and kk in noroute_hrus_by_loc:
    #     config.update_value('hru_noroute', noroute_hrus_by_loc[kk])

    if kk in region_outlet_cutoffs:
        config.update_value('cutoffs', region_outlet_cutoffs[kk])
        
    config.write(f'{bandit_work_dir}/bandit_v11_gm_r{kk}.cfg')

# %%
