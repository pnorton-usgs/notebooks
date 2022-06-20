# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python [conda env:gis_38]
#     language: python
#     name: conda-env-gis_38-py
# ---

# %%
from osgeo import ogr
import os

# %%
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20220214_gm_pipestem/GIS'
filename = f'{work_dir}/HRU_subset.shp'

# %%
driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(filename, 0)
layer = dataSource.GetLayer()

for feature in layer:
    geom = feature.GetGeometryRef()
    print(geom.Centroid().ExportToWkt())

# %%
print(geom.Centroid())

# %%

# %%
