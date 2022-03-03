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
import pandas as pd
import geopandas as gd
from shapely.geometry import Point
from shapely import wkt

from pyPRMS.ParamDb import ParamDb

# %%
# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/FY2021_projects/modeling_fabric/GFv2.0'
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/FY2021_projects/modeling_fabric/GFv2.0/20210812_POIs'
gfv2_file = f'{work_dir}/GFv2_gages_cleaned.csv'
# gfv2_file = f'{work_dir}/20210713_gages.csv'

paramdb_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS'

# %%
# For 20210713_gages.csv
# col_names = ["objectid", "COMID", "Type_HUC12", "Type_Gages", "Type_TE", "Type_NID",
#              "Type_WBIn", "Type_WBOut", "Type_Conf", "longitude", 'latitude']
# col_types = [int, int, str, str, str, str, str, str, str, float, float]

# For 20210812 GFv2_gages_cleaned.csv
col_names = ['idx', 'Gage_no', 'COMID.x', 'REACHCODE', 'REACH_meas', 'drain_area', 
             'contrib_dr', 'begin_date', 'end_date', 'count_nu', 'CLASS', 'rank_id', 
             'geom', 'Nav_COMID', 'TotDASqKM_orig', 'TotDASqKM_Nav', 'COMID.y', 
             'ref_ID', 'rpu', 'TotDASqKM_Ref', 'seg_id', 'split', 'DAR_Nav', 'DAR_Ref']

# NOTE: Must use Int64 dtype for nullable-integer fields
col_types = [int, str, pd.Int64Dtype(), str, float, float, float, str,
             str, pd.Int64Dtype(), str, pd.Int64Dtype(), str, pd.Int64Dtype(), float, float,
             pd.Int64Dtype(), pd.Int64Dtype(), str, float, pd.Int64Dtype(), str, float, float]

cols = dict(zip(col_names, col_types))

# Handling date fields: https://stackoverflow.com/questions/21269399/datetime-dtypes-in-pandas-read-csv
date_flds = ['begin_date', 'end_date']

df = pd.read_csv(gfv2_file, sep=',', dtype=cols, parse_dates=date_flds)
# gdf = gd.read_file(f'{work_dir}/20210713_gages.csv')

# %%
df.info()

# %%
df.head()

# %%
# For a discussion on converting a dataframe with WKT coordinates to geopandas see:
#    https://geopandas.readthedocs.io/en/latest/gallery/create_geopandas_from_pandas.html
df['geom'] = gd.GeoSeries.from_wkt(df['geom'])

gdf = gd.GeoDataFrame(df, crs='epsg:4326', geometry=df.geom)

# Create a geodataframe from pandas having lat/lon fields
# gdf = gd.GeoDataFrame(df.drop(['longitude', 'latitude'], axis=1), 
#                       crs='epsg:4326', geometry=[Point(xy) for xy in zip(df.longitude, df.latitude)])

# %%
gdf.info()

# %%

# %%
gdf.plot()

# %%
gdf.head()

# %%

# %%

# %%
# %%time
pdb = ParamDb(paramdb_dir)

# %%
poi_gages = pdb.parameters['poi_gage_id'].data.tolist()

# %%

# %%
# v2_gages = df.Type_Gages.tolist()

v2_gages = df.Gage_no.tolist()

# %%
missing_from_v2 = set(poi_gages).difference(set(v2_gages))

# %%
len(missing_from_v2)

# %%
missing_from_v2

# %%
hydat_gages = []
nwis_gages = []

for xx in missing_from_v2:
    if len(xx) ==7:
        hydat_gages.append(xx)
    else:
        nwis_gages.append(xx)

# %%
len(hydat_gages)

# %%
len(nwis_gages)

# %%
nwis_gages.sort()

# %%
v2_gages.index('06755960')

# %%
nwis_gages

# %%
# Write the missing gages to file
ohdl = open(f'{work_dir}/tmp_v11_missing_v2.csv', 'w')
ohdl.write('gage_no,miss_flag\n')

for xx in nwis_gages:
    ohdl.write(f'{xx},1\n')
    
ohdl.close()

# Create matching csvt file
ohdl = open(f'{work_dir}/tmp_v11_missing_v2.csvt', 'w')
ohdl.write('String,Integer')
ohdl.close()


# %%
