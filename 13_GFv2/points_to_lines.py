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
#     display_name: Python [conda env:gis_38]
#     language: python
#     name: conda-env-gis_38-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import pandas as pd
import numpy as np

from collections import OrderedDict

# %%
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/FY2021_projects/modeling_fabric'

st_filename = f'{work_dir}/seg_startpoints.csv'
en_filename = f'{work_dir}/seg_endpoints_v2.csv'

# %%
# cols = ['fid', 'OBJECTID', 'nsegment_v1_1', 'tosegment_v1_1', 'seg_id_nhm', 
#         'Version', 'Shape_Length', 'seg_cum_area_seg_cum_area', 'vertex_pos',
#         'vertex_index', 'vertex_part', 'vertex_part_index', 'distance', 'angle',
#         'xcoord', 'ycoord', 'zcoord', 'mvalue']

st_df = pd.read_csv(st_filename, sep=',', usecols=[2, 14, 15])
st_df.rename(columns={'x': 'st_lon', 'y': 'st_lat'}, inplace=True)
st_df.set_index('nsegment_v1_1', inplace=True)

# en_df = pd.read_csv(en_filename, sep=',', usecols=[2, 3, 14, 15])
en_df = pd.read_csv(en_filename, sep=',', usecols=[0, 1, 2, 3])
en_df.rename(columns={'x': 'en_lon', 'y': 'en_lat'}, inplace=True)
en_df.set_index('nsegment_v1_1', inplace=True)

# %%
en_df.info()

# %%
st_df.info()

# %%
df = st_df
df = pd.merge(df, en_df, how='left', left_index=True, right_index=True)
# poi_info_df = pd.merge(poi_info_df, poi_df, how='left', left_index=True, right_index=True)

# %%
df.info()

# %%
df.head()

# %%
# Create maps for coordinates and tosegs
# st_seg_to_coords = OrderedDict()
# st_seg_toseg = OrderDict()

# en_seg_to_coords = OrderedDict()

# for xx in df.index.tolist():
#     seg_to_coords[xx] = [df.loc[xx, 'xcoord'], df.loc[xx, 'ycoord']]
#     seg_toseg[xx] = df.loc[xx, 'tosegment_v1_1']

# %%
from math import sin, cos, sqrt, atan2, radians

R = 6373.0   # approx radius of earth in km

seg_pts = OrderedDict()
pt_pairs = OrderedDict()
dist_to_seg = OrderedDict()

for xx in df.index.tolist():
    lat_d = df.loc[xx, 'en_lat']
    lon_d = df.loc[xx, 'en_lon']
    lat1 = radians(df.loc[xx, 'en_lat'])
    lon1 = radians(df.loc[xx, 'en_lon'])
    toseg = df.loc[xx, 'tosegment_v1_1']
    
    if toseg == 0:
        continue
    
    lat2 = radians(df.loc[toseg, 'st_lat'])
    lon2 = radians(df.loc[toseg, 'st_lon'])
    
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    a = sin(dlat / 2.0)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2.0)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    
    dist_to_seg[xx] = R * c
    pt_pairs[xx] = toseg
    seg_pts[xx] = [lon_d, lat_d]

# %%
dist_to_seg

# %%
# For each segment write the distance between the endpoint of the segment to the starting point of the next segment
fhdl = open(f'{work_dir}/seg_dist.csv', 'w')

fhdl.write('seg_id,toseg_dist\n')

for xx, yy in dist_to_seg.items():
    fhdl.write(f'{xx},{yy}\n')
    
fhdl.close()

# %%
# Write pairs of points for each seg
fhdl = open(f'{work_dir}/seg_pairs.csv', 'w')

fhdl.write('seg_id,seg_id2\n')

for xx, pp in pt_pairs.items():
    fhdl.write(f'{xx},{pp}\n')
#     fhdl.write(f'{xx},{pp[0]},{pp[1]}\n')
#     fhdl.write(f'{xx},{pp[2]},{pp[3]}\n')

fhdl.close()


# %%
# Write the lat/lon for each segment
fhdl = open(f'{work_dir}/seg_pts.csv', 'w')

fhdl.write('seg_id,x,y\n')

for xx, pp in seg_pts.items():
    fhdl.write(f'{xx},{pp[0]},{pp[1]}\n')

fhdl.close()


# %%

# %%
import ogr

# %%
multiline = ogr.Geometry(ogr.wkbLineString)

# %%
for ss, cc in seg_to_coords.items():
    if seg_toseg[ss] == 0:
        continue
    
    pt1 = cc
    pt2 = seg_to_coords[seg_toseg[ss]]
    
    print(pt1, pt2)

# %%
