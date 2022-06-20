# ---
# jupyter:
#   jupytext:
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

# %%
# import math
# import resource
# import sys
import numpy as np
import pandas as pd
import netCDF4 as cf
import time

from numpy import cos, sin
from scipy.spatial import cKDTree


# %%
class Kdtree_fast(object):
    def __init__(self, ncfile, latvarname, lonvarname):
        self.ncfile = ncfile

        # Need to check for time-varying spatial coordinates (e.g. WRF output)
        if self.ncfile.variables[latvarname].ndim == 1:
            # 1D latitude array
            lats_1d = self.ncfile.variables[lonvarname]
        elif self.ncfile.variables[latvarname].ndim > 2:
            self.latvar = self.ncfile.variables[latvarname][0, ...]
        else:
            self.latvar = self.ncfile.variables[latvarname]

        if self.ncfile.variables[latvarname].ndim == 1:
            # 1D longitude array
            lons_1d = self.ncfile.variables[latvarname]
        elif self.ncfile.variables[lonvarname].ndim > 2:
            self.lonvar = self.ncfile.variables[lonvarname][0, ...]
        else:
            self.lonvar = self.ncfile.variables[lonvarname]

        if self.ncfile.variables[latvarname].ndim == 1:
            # For 1D lat/lon create 2D lat and 2D lon arrays
            self.latvar, self.lonvar = np.meshgrid(lons_1d, lats_1d, indexing='ij')

        # Read latitude and longitude from file into numpy arrays
        self.latvals = np.radians(self.latvar)
        self.lonvals = np.radians(self.lonvar)

        self.shape = self.latvals.shape

        clat, clon = cos(self.latvals), cos(self.lonvals)
        slat, slon = sin(self.latvals), sin(self.lonvals)

        triples = zip(np.ravel(clat*clon), np.ravel(clat*slon), np.ravel(slat))
        self.kdt = cKDTree(list(triples))
        del clat
        del clon
        del slat
        del slon

    def query(self, lat0, lon0):
        lat0_rad = np.radians(lat0)
        lon0_rad = np.radians(lon0)

        clat0, clon0 = cos(lat0_rad), cos(lon0_rad)
        slat0, slon0 = sin(lat0_rad), sin(lon0_rad)
        dist_sq_min, minindex_1d = self.kdt.query([clat0*clon0, clat0*slon0, slat0])
        # print "minindex_1d:", minindex_1d
        # print "dist_sq_min:", dist_sq_min
        # print "self.shape:", self.shape
        iy_min, ix_min = np.unravel_index(minindex_1d, self.shape)

        return iy_min, ix_min


# %%
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet_v3'
filename = f'{work_dir}/daymet_v3_landmask.nc'

# %%
# %%time
fhdl = cf.Dataset(filename, 'r')

lats = fhdl.variables['lat'][:]
lons = fhdl.variables['lon'][:]

# %%
# %%time
ns = Kdtree_fast(fhdl, 'lat', 'lon')

# %%

# %%
# compute neighborhood window
lat_winsize = 2
lon_winsize = 2
weight = np.zeros(lat_winsize * lon_winsize)
vals = np.zeros(lat_winsize * lon_winsize)
dist = np.zeros(lat_winsize * lon_winsize)
p = 2   # Arbitrary power parameter
print("lat_winsize =", lat_winsize)
print("lon_winsize =", lon_winsize)

# %%
test_lat = 47.284431
test_lon = -99.35701

# %%
lat_idx, lon_idx = ns.query(test_lat, test_lon)

print(lat_idx, lon_idx)

# %%
# Find the bounding box of cells around the given point
if test_lat > lats[lat_idx, lon_idx]:
    latstart_idx = lat_idx
    latend_idx = lat_idx + 1
else:
    latstart_idx = lat_idx - 1
    latend_idx = lat_idx

if test_lon > lons[lat_idx, lon_idx]:
    lonstart_idx = lon_idx
    lonend_idx = lon_idx + 1
else:
    lonstart_idx = lon_idx - 1
    lonend_idx = lon_idx
    
print(latstart_idx, latend_idx)
print(lonstart_idx, lonend_idx)

# %%
# Compute the weights for this point
tlons = lons[latstart_idx:latend_idx+1, lonstart_idx:lonend_idx+1]
tlats = lats[latstart_idx:latend_idx+1, lonstart_idx:lonend_idx+1]

# Compute the distance between the given point and each cell
dist = (tlons - test_lon)**2 + (tlats - test_lat)**2
dist = np.sqrt(dist)
dist_total = np.sum(dist**-p)

# Compute the weights for each cell
weights = dist**-p / dist_total

# %%
print(weights)

# %%

# %%
# Interpolate the data
interp_data = np.sum(thevals[:, latstart_idx:latend_idx+1, lonstart_idx:lonend_idx+1]*weights, axis=(1, 2))


# %%

# %%

# %%

# %%
