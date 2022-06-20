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
# import hydrodata as hd
# from hydrodata import NWIS, plot
import pygeohydro as gh
# from pygeohydro import NWIS, plot
from pynhd import NLDI, NHDPlusHR, WaterData

# %% [markdown]
# # Search NLDI by navigating above USGS-06469400
# - grab all comids and gages on main flowline
# - use station comids to filter NHM POIs

# %%
reach_id = '948100389'
lim_distance = 1

nldi_flw_main = NLDI().navigate_byid(fsource="comid", fid=reach_id, 
                                     navigation="upstreamMain", source="flowlines", distance=lim_distance)

nldi_st_all = NLDI().navigate_byid(fsource="comid", fid=reach_id, 
                                   navigation="upstreamMain", source="nwissite", distance=lim_distance)

# %%
nldi_st_all.head()

# %%

# %%
poi_id = '06469400'
gage_id = f'USGS-{poi_id}'
lim_distance = 150

# fsource: comid, huc12pp, nwissite, wade, WQP
# fid: why does USGS- prepended to the gage ID?
# navigation: upstreamMain, upstreamTributaries
# source: return the data from another source after navigating the features using fsource
# distance: (km) limit search for navigation up to distance

nldi_flw_main = NLDI().navigate_byid(fsource="nwissite", fid=gage_id, 
                                     navigation="upstreamMain", source="flowlines", distance=lim_distance)

nldi_st_all = NLDI().navigate_byid(fsource="nwissite", fid=gage_id, 
                                   navigation="upstreamMain", source="nwissite", distance=lim_distance)

nldi_POI_all = NLDI().navigate_byid(fsource="nwissite", fid=gage_id, 
                                    navigation="upstreamMain", source="gfv11_pois", distance=lim_distance)

# Get the tributaries
nldi_flw_tribs = NLDI().navigate_byid(fsource="nwissite", fid=gage_id, 
                                      navigation="upstreamTributaries", source="flowlines", distance=lim_distance)

# Get the basin
# ?? why does this use the raw station ID instead of the USGS-* id?
basin = NLDI().get_basins(poi_id)

# Get the HUC12 pour points
pp = NLDI().navigate_byid(fsource="nwissite", fid=gage_id, navigation="upstreamTributaries", 
                          source="huc12pp",distance=lim_distance)

# %%
nldi_flw_main.head()

# %%
nldi_st_all.head()

# %%
nldi_POI_all.head(10)

# %%
# Plot the flowlines
ax = nldi_flw_main.plot(lw=1, color="darkblue", zorder=2, label="Mainstem")

nldi_st_all.plot(ax = ax, label="USGS stations", marker="o", markersize=3, zorder=4, color="r")
nldi_st_all.plot(ax = ax, label="gv11_pois", marker="+", markersize=40, zorder=3, color="black")

ax.legend(bbox_to_anchor=(0.6, -0.15))
ax.set_aspect("equal")
ax.figure.set_dpi(100)

# %% [markdown]
# ## Get comid for the gage-of-interest

# %%
# Get comid for gage of interest USGS-06469400
# st = nldi_st_all.loc[nldi_st_all['identifier'] == gage_id]
# st_comid = st['comid'].values[0]

st_comid = nldi_st_all.loc[nldi_st_all['identifier'] == gage_id]['comid'].values[0]
print(st_comid)

# %% [markdown]
# ## Find the NHM identifier in nldi_POI_all that matches the gage-of-interest comid

# %%
# find identifier in nldi_POI_all that has comid = gage of interest
poi = nldi_POI_all.loc[nldi_POI_all['comid'] == st_comid]

# NHM stream segment is the identifier in poi
# NOTE: for some reason the identifier is stored as a string instead of an integer
nhm_seg_id = [int(xx) for xx in poi['identifier'].values]
print(nhm_seg_id)

# %%
poi

# %%
ax = basin.plot(facecolor="none", edgecolor="k", figsize=(8, 8))
nldi_st_all.plot(ax=ax, label="USGS stations", marker="*", markersize=300, zorder=4, color="b")

# st_d20.plot(ax=ax, label="USGS stations up to 20 km", marker="v", markersize=100, zorder=5, color="darkorange",)

pp.plot(ax=ax, label="HUC12 pour points", marker="o", markersize=50, color="k", zorder=3)

nldi_flw_main.plot(ax=ax, lw=3, color="r", zorder=2, label="Main")
nldi_flw_tribs.plot(ax=ax, lw=1, zorder=1, label="Tributaries")
ax.legend(loc="best")
ax.set_aspect("auto")
ax.figure.set_dpi(100)
# ax.figure.savefig("_static/nhdplus_navigation.png", bbox_inches="tight", facecolor="w")

# %%

# %%
from pyPRMS.ParamDb import ParamDb
from Bandit.bandit_helpers import parse_gages, set_date, subset_stream_network

# %%
paramdb_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS'


# %%
print('Read parameter database')
pdb = ParamDb(paramdb_dir)
nhm_params = pdb.parameters
nhm_global_dimensions = pdb.dimensions

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get tosegment_nhm
nhm_seg = nhm_params.get('nhm_seg').tolist()

# Build the stream network
print('Building stream network')
dag_ds = pdb.parameters.stream_network(tosegment='tosegment_nhm', seg_id='nhm_seg')

# %%
print('Pulling subset of stream network')
uscutoff_seg = []
dsmost_seg = nhm_seg_id
dag_ds_subset = subset_stream_network(dag_ds, uscutoff_seg, dsmost_seg)

# %%
# Create list of toseg ids for the model subset
ss_nhm_seg = list(set(xx[0] for xx in dag_ds_subset.edges))
print(f'Number of segments in subset: {len(ss_nhm_seg)}')



# %%
ss_nhm_seg

# %%
nldi_flw_tribs

# %%

# %%
st_comid = nldi_st_all.loc[nldi_st_all['identifier'] == gage_id]['comid'].values[0]
print(st_comid)

# %%

# %%

# %%
poi = nldi_POI_all.loc[nldi_POI_all['comid'] == st_comid]

nhm_seg_id = [int(xx) for xx in poi['identifier'].values]
print(nhm_seg_id)

# %%
poi

# %%

# %%
# NLDI().navigate_byid(fsource="gfv11_pois", fid=gage_id, 
#                      navigation="upstreamMain", source="gfv11_pois", distance=lim_distance)

# %%
# Get dictionary of all valid feature sources for NLDI
NLDI().valid_fsources

# %%
