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
#     display_name: Python [conda env:gis]
#     language: python
#     name: conda-env-gis-py
# ---

# %%
import logging

import fiona




# %%
shpfile = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'

logging.basicConfig(level=logging.DEBUG)

with fiona.Env(CPL_DEBUG=True):
    fiona.open(shpfile)

# %%
shapes = fiona.open(shpfile)

# %%
shapes.crs

# %%
shapes.crs_wkt

# %%
shapes.bounds

# %%
shapes.meta

# %%
aa = shapes.crs_wkt
type(aa)

# %%
for xx in aa.split(','):
    print(xx)
#     for yy in xx.split(','):
#         print(yy)
        
#         for zz in yy.split('['):
#             print(zz)

# %%
