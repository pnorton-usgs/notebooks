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
from collections import OrderedDict

from Bandit.bandit_multi_locations import read_file

# %%
xwalk_file = '/Users/pnorton/Projects/National_Hydrology_Model/GIS/nhm_comid_crosswalk_v1.0.csv'

req_file = '/Users/pnorton/Projects/National_Hydrology_Model/FY2022_projects/extraction_requests/2021-10-19_vlah/priority1_COMIDs.csv'
req2_file = '/Users/pnorton/Projects/National_Hydrology_Model/FY2022_projects/extraction_requests/2021-10-19_vlah/priority2_COMIDs.csv'

# %%
# def read_file(filename):
#     """Read a csv"""

#     values_by_id = OrderedDict()

#     src_file = open(filename, 'r')
#     src_file.readline()

#     # Read in the non-routed HRUs by location
#     for line in src_file:
#         cols = line.strip().replace(' ', '').split(',')
        
#         # Assume first column is a number
#         cols = [int(xx) for xx in cols]

#         if cols[0] in values_by_id:
#             print(f'{cols[0]} already exists; additional seg_id: {cols[1]}')
#         else:
#             values_by_id[cols[0]] = cols[1]

#     return values_by_id

# %%
# Read the COMID to nhm_seg crosswalk file
comids = read_file(xwalk_file)

# %%
fhdl = open(req_file, 'r')
fhdl.readline()

reqs = []
for row in fhdl:
    reqs.append(int(row.strip()))
    
for rr in reqs:
    try:
        print(f'{rr}: {comids[rr]}')
    except KeyError:
        print(f'{rr}: NOT FOUND')    

# %%

# %%
fhdl = open(req2_file, 'r')
fhdl.readline()

reqs2 = []
for row in fhdl:
    reqs2.append(int(row.strip()))
    
for rr in reqs2:
    try:
        print(f'{rr}: {comids[rr]}')
    except KeyError:
        print(f'{rr}: NOT FOUND')    

# %%
