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
# import fiona

import geopandas
import networkx as nx
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib as mpl  

# %%
# pnc.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)
# geopandas.read_file(filename, layer=layer_name)

# %%
hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
hru_layer_name = 'nhruv11_sim30'
hru_shape_key='nhru_v11'

# %%
# c = fiona.open(hru_geodatabase, 'r')
# len(list(c))

# help(fiona.open)

# %%
# %%time
xx = geopandas.read_file(hru_geodatabase, layer=hru_layer_name)

# %%
xx.crs.name

# %%
# %%time
xx.crs = 'EPSG:5070'

# %%
xx.crs.name

# %%
xx.info()

# %%

# %%
import cartopy
import cartopy.feature as cfeature
import cartopy.crs as ccrs

import cartopy.io.shapereader as shpreader

# %%
resol = '110m'
shpfilename = shpreader.natural_earth(resolution=resol,
                                      category='physical',
                                      name='lakes')
reader = shpreader.Reader(shpfilename)

greatlakes = [lake for lake in reader.records() if lake.attributes['name_alt'] == 'Great Lakes']
# gl_feature = cartopy.feature.ShapelyFeature([greatlakes.geometry], crs_proj,
#                                             edgecolor=water_style['edgecolor'],
#                                             facecolor=water_style['facecolor'])

# %%
for lake in reader.records():
    if lake.attributes['name_alt'] == 'Great Lakes':
        print(lake.attributes['name_alt'])

# %%
greatlakes

# %%

# %%

# %%

# %%
