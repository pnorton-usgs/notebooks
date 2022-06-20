# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python [conda env:holoview]
#     language: python
#     name: conda-env-holoview-py
# ---

# %%
import geopandas

# %%
hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
hru_layer_name = 'nhruv11_sim30'
hru_shape_key='nhru_v11'

# %%
gdf = geopandas.read_file(hru_geodatabase, layer=hru_layer_name)

# %%
gdf.crs.name

# %%
gdf.crs

# %%
gdf.crs.coordinate_operation.method_code

# %%
gdf.geometry.total_bounds

# %%
gdf.crs = 'EPSG:5070'

# %%

# %%
import pyproj

from shapely.geometry import Point
from shapely.ops import transform

wgs84_pt = Point(-72.2495, 43.886)

wgs84 = pyproj.CRS('EPSG:4326')
utm = pyproj.CRS('EPSG:32618')

project = pyproj.Transformer.from_crs(wgs84, utm, always_xy=True).transform
utm_point = transform(project, wgs84_pt)
print(utm_point)

# %%

# %%
aea_pt = Point(-2356125.3294, 209654.8967)

aea = pyproj.CRS('EPSG:5070')

project = pyproj.Transformer.from_crs(aea, wgs84, always_xy=True).transform
wgs_point = transform(project, aea_pt)
print(wgs_point)

# %% [markdown]
# ## Transform to lon/lat with longitude between 0 and 360

# %%
import pyproj

from shapely.geometry import Point
from shapely.ops import transform

aea_pt = Point(-2356125.3294, 209654.8967)
aea = pyproj.CRS('EPSG:5070')

# Can't use EPSG:4326 for this
# if lon_wrap=0 then longitude is -180 to 180
# if lon_wrap=180 then longitude is 0 to 360
latlon = pyproj.CRS.from_user_input("+proj=longlat +lon_wrap=180")

xformer = pyproj.Transformer.from_crs(aea, latlon, always_xy=True).transform
wgs_point = transform(xformer, aea_pt)

print(wgs_point)

# %%

# %%
# %%time
gdf2 = gdf.to_crs(latlon)

# %%
gdf2.crs.name

# %%
gdf2.geometry.total_bounds

# %%
