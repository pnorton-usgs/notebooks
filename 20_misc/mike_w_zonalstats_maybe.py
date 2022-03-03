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

# %%
import rasterstats
from rasterstats import zonal_stats
from affine import Affine
import time
import rasterio
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal, ogr
from osgeo.gdalconst import *
from osgeo import gdal
# from gdal import *
import numpy as np
import sys
import json
import csv
import pandas as pd
import os
import math
import scipy
from scipy import stats
import os.path
from os import path
import pathlib
from pathlib import Path
driver = ogr.GetDriverByName("FileGDB")

# %%

# %%
hru = gpd.read_file('G:/WBEEP/SB_Layers/TGF.gdb', layer = 'nhru')

# %%
asp = rasterio.open("G:/WBEEP/OpenSource/aspect.tif")

# %%
asp_a = asp.read()

# %%
asp_a = np.ma.masked_equal(asp_a, -9999)

# %%
in_asp_sin = np.sin(asp_a * (math.pi / 180.0))

# %%
in_asp_sin.max()

# %%
in_asp_cos = np.cos(asp_a * (math.pi / 180.0))

# %%
in_asp_cos.max()

# %%
with rasterio.open(r'G:/WBEEP/OpenSource/Untitled_Folder/nhrug.tif') as src:
    transformb = src.meta['transform']
    print(type(transformb), src.meta)
    in_asp_cos = src.read(1)

# %%
hrudata = gpd.GeoDataFrame.from_file(r'G:/WBEEP/OpenSource/Untitled_Folder/nhru.shp')

# %%
cos_stats = zonal_stats(hrudata, in_asp_cos, transform=transformb.to_gdal(),prefix='cos_', 
                    all_touched=False, geojson_out=True, nodata=src.meta['nodata'], stats="sum")

# %%
stats_gdf = gpd.GeoDataFrame.from_features(cos_stats)

# %%
stats_gdf.head()

# %%
cos_stats

# %%
with rasterio.open(r'G:/WBEEP/OpenSource/Untitled_Folder/nhrug.tif') as src:
    transformb = src.meta['transform']
    print(type(transformb), src.meta)
    in_asp_sin = src.read(1)

# %%
sin_stats = zonal_stats(hru, in_asp_sin, transform=transformb.to_gdal(),  stats="sum")

# %%
sin_stats
