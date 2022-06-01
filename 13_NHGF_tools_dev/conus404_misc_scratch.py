# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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

# %%
import dask
import datetime
import fsspec
import numpy as np
import pandas as pd
import time
import xarray as xr
import geoviews as gv
import hvplot.xarray

import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf

from typing import Dict, Optional, Union

from dask.distributed import Client

# %%
work_dir = '/Volumes/USGS_NHM2/datasets/wrf-conus404/originals/WY2021'
filename = f'{work_dir}/wrfxtrm_d01_2020-10-01_00:00:00'

# work_dir = '/Volumes/USGS_NHM2/datasets/wrf-conus404/originals'
# filename = f'{work_dir}/wrfout_d01_2020-09-30_00:00:00'

# %%
client = Client()

# %%
client

# %%
xdf = xr.open_dataset(filename, chunks={})
xdf

# %%
df = xr.open_mfdataset(filename, chunks={}, coords="none", data_vars="minimal",
                                 compat='override', parallel=True)
df

# %%
df = xr.open_mfdataset(f'{work_dir}/wrfxtrm_d01_*', concat_dim='Time', combine='nested', parallel=True, coords="minimal", 
                       data_vars="minimal", compat='override', chunks={})
df

# %%

# %%

# %%

# %%

# %% [markdown] tags=[]
# ## Convert Times bytes to datetime
# This is necessary for the daily model output which is missing the XTIME variable. End result of this cell is a time coordinate variable.

# %% tags=[]
t_str = df['Times'].values
new_times = [datetime.datetime.strptime(tt.tobytes().decode('ascii'), '%Y-%m-%d_%H:%M:%S') for tt in t_str]

df['time'] = (('Time'), new_times)
df = df.rename({'Time': 'time'})
df = df.assign_coords({'time': df.time})
df = df.drop('Times')
df

# %%

# %%

# %%

# %%
df_lat = df['XLAT']
df_lon = df['XLONG']

df_lat = df_lat.squeeze('Time', drop=True)
df_lon = df_lon.squeeze('Time', drop=True)

df['XLAT'] = df_lat
df['XLONG'] = df_lon

# Remove the time dimension from the lat/lon-related variables
# df2 = df.assign_coords(XLAT=df.coords['XLAT'].squeeze('Time'), XLONG=df.coords['XLONG'].squeeze('Time'))


# %%

# %%
df

# %%

# %%

# %%
# %%time
rename_dims = {'south_north': 'y', 'west_east': 'x',
               'south_north_stag': 'y_stag', 'west_east_stag': 'x_stag',
               'Time': 'time'}

rename_vars = {'XLAT': 'lat', 'XLAT_U': 'lat_u', 'XLAT_V': 'lat_v',
               'XLONG': 'lon', 'XLONG_U': 'lon_u', 'XLONG_V': 'lon_v'}

# var_metadata['time'] = dict(axis='T', standard_name='time')

df = df.rename(rename_vars)
df = df.rename(rename_dims)
df = df.assign_coords({'time': df.XTIME})

# Modify the attributes
for cvar in df.variables:
    # Remove unneeded attributes, update the coordinates attribute
    for cattr in list(df[cvar].attrs.keys()):
        if cattr == 'coordinates':
            # Change the coordinates attribute to new lat/lon naming
            orig_coords = df[cvar].attrs[cattr]
            new_coords = []
            for xx in orig_coords.split(' '):
                if xx not in rename_vars:
                    continue
                new_coords.append(rename_vars[xx])
            df[cvar].attrs[cattr] = ' '.join(new_coords)

# %%
# df = df.assign_coords({'time': df.XTIME})
# df = df.rename({'Time': 'time'})

df['time'].attrs['axis'] = 'T'
df['time'].attrs['standard_name'] = 'time'

# %%
df

# %%

# %%

# %% jupyter={"outputs_hidden": true} tags=[]
# %%time
for vv in df.variables:
    cvar = vv
    try:
        minval = df[cvar].min().values
        maxval = df[cvar].max().values
        meanval = df[cvar].mean().values
        print(f'{cvar}: {minval}, {maxval}, {meanval}')
    except TypeError:
        print(f'{cvar}: TypeError')

# %%
df.variables.keys()

# %% [markdown]
# There are some cells where vegation leaf temperature (TV; units=K) is negative

# %%
df['NEG_TV'] = df.TV.where(df.TV < 0.0)

print(f'Count of negative vegetation leaf temperature (TV) values: {df.NEG_TV.count().values}')

# %%
# %%time
df['NEG_TV'].hvplot.quadmesh(x='lon', y='lat', geo=True, rasterize=True, cmap='turbo', tiles='OSM')

# %% [markdown]
# There are some cells where leaf area index (LAI; units=m2 m-2) is negative

# %%
df['NEG_LAI'] = df.LAI.where(df.LAI < 0.0)

print(f'Count of negative leaf area index (LAI) values: {df.NEG_LAI.count().values}')

# %%
# %%time
df['NEG_LAI'].hvplot.quadmesh(x='lon', y='lat', geo=True, rasterize=True, cmap='turbo', tiles='OSM')

# %% [markdown]
# ## WRF constants

# %%
work_dir = '/Volumes/USGS_NHM2/datasets/wrf-conus404/originals'
filename = f'{work_dir}/wrf_constants_conus404_update_v3.nc'

xdf = xr.open_dataset(filename, chunks={})
xdf

# %%

# %% [markdown]
# ### Plot using projection coordinates

# %%
# Plot using projection coordinates
xdf['LANDMASK'].west_east.plot()

# %% [markdown]
# ### Plot using projection coordinates transformed to lat/lon

# %%
varname = 'LANDMASK'
crs_info = xdf.crs

xx = xdf.west_east.values
yy = xdf.south_north.values

globe = ccrs.Globe(ellipse='sphere', semimajor_axis=6370000, semiminor_axis=6370000)
lcc = ccrs.LambertConformal(globe=globe, 
                            central_longitude=crs_info.longitude_of_central_meridian, 
                            central_latitude=crs_info.latitude_of_projection_origin,
                            standard_parallels=crs_info.standard_parallel)
ax = plt.axes(projection=ccrs.PlateCarree())

xdf[varname].plot.contourf(ax=ax, transform=lcc)
ax.coastlines()

ax.add_feature(cartopy.feature.BORDERS, linestyle='-');
ax.set_extent([xx.min(), xx.max(), yy.min(), yy.max()], crs=lcc)
# ax.set_extent([xx.min(), xx.max(), yy.min(), yy.max()], crs=lcc)

# %%

# %%

# %%

# %%

# %%
thresh = 1e9
acswdnb_0 = 8.402427e+08
i_acswdnb_0 = 290.0

acswdnb_1 = 8.420878e+08
i_acswdnb_1 = 290.0

swdnb = (acswdnb_1 + thresh * i_acswdnb_1) - (acswdnb_0 + thresh * i_acswdnb_0)
print(swdnb)

# %%
print((acswdnb_1 + thresh * i_acswdnb_1) / 3600.)

# %%
swdnb / 3600.

# %%
import math


# %%
def ncl_wrf_rh(qv, pres, temp):
    es0 = 0.6112 / 0.001   # [Pa]
    # svp1 = 0.6112   # [kPa] 
    svp2 = 17.67
    svp3 = 29.65
    svpt0 = 273.15
    epsilon = 0.622   # Rd / Rv; []
        
    # Vapor pressure [Pa]
    es = es0 * math.exp((svp2 * (temp - svpt0)) / (temp - svp3))   # 
    
    # Saturation specific humidity
    qs = (epsilon * es) / (pres - es * (1.0 - epsilon))
    
    # Saturation mixing ratio
    # NOTE: Original ncl code computed rh using mixing ratio and saturation specific humidity which
    #       produces slightly different rh values
    qvs = qs / (1 + qs)
    
    # es = 10 * svp1 * math.exp((svp2 * (temp - svpt0)) / (temp - svp3))   # 
    # qvs = (epsilon * es) / (0.01 * pres - (1.0 - epsilon) * es)
    # sh = qv / qvs
    rh = 100 * max(min(qv / qvs, 1.0), 0.0)
    
    print(f'{es=}')
    print(f'{qs=}')
    print(f'{qvs=}')
    # print(f'{sh=}')
    print(f'{rh=}')

    
def wrf_rh(qv, pres, temp):
    # from: https://forum.mmm.ucar.edu/phpBB3/viewtopic.php?t=9134#:~:text=Re%3A%20relative%20humidity&text=The%20Relative%20Humidity%20(RH)%20is,you%20can%20easily%20obtain%20RH.
    
    # real, parameter :: svp1=611.2
    # real, parameter :: svp2=17.67
    # real, parameter :: svp3=29.65
    # real, parameter :: svpt0=273.15
    # real, parameter :: eps = 0.622
    # rh = 1.E2 * (p*q/(q*(1.-eps) + eps))/(svp1*exp(svp2*(t-svpt0)/(T-svp3)))
    # rh = 100 * (p * q / (q * (1.0 - 0.622) + 0.622)) / (611.2 * exp(17.67 * (t - 273.15) / (t - 29.65)
    es0 = 611.2   # Saturation vapor pressure reference value [Pa]
    svp2 = 17.67
    svp3 = 29.65
    svpt0 = 273.15
    epsilon = 0.622
    
    # e = qv * pres / (qv - 0.622 * qv) + 0.622)
    # e = qv * pres / (0.378 * qv + 0.622)
    
    e = pres * qv / (qv + epsilon)
    
    # NOTE: WRF uses the following equation for vapor pressure. I'm not sure where this derivation came from.
    # e = pres * qv / (qv * (1.0 - epsilon) + epsilon)
    
    # Saturation vapor pressure, Bolton formula (1980)
    # es = 611.2 * exp(17.67 * (temp - 273.15) / (temp - 29.65)
    es = es0 * math.exp(svp2 * (temp - svpt0) / (temp - svp3))
    
    rh = 100 * e / es
    # rh = 100 * (pres * qv / (qv * (1.0 - eps) + eps)) / es
    # rh = 100 * (pres * qv / (qv * (1.0 - eps) + eps)) / (svp1 * math.exp(svp2 * (temp - svpt0) / (temp - svp3)))
    
    print(f'{e=}')
    print(f'{es=}')
    print(f'{rh=}')
    return rh

def teten_rh(qv, pres, temp):
    # Teten's eqn
    es0 = 6.113  # [hPa]
    es0 *= 100   # [Pa]
    epsilon = 0.622  # Rd / Rv; []
    
    temp_c = temp - 273.15   # [C]
    
    # Vapor pressure
    e = (qv * pres) / (epsilon + qv)   # [Pa]
    
    # Saturation vapor pressure
    es = es0 * math.exp(17.269 * temp_c / (temp_c + 237.3))   # [Pa]
    # es = es0 * math.exp((17.269 * (temp - 273.15)) / (temp - 35.86))   # [Pa]
    rh = 100.0 * (e / es)   # 0-100%
    
    print(f'{e=}')
    print(f'{es=}')
    print(f'{rh=}')
    
def wrf_diag_rh(qv, pres, temp):
    # from calc_rh() in model_diag_afwa.F
    pq0 = 379.90516   # How was this derived?
    a2 = 17.2693882
    a3 = 273.16
    a4 = 35.86
    
    # Water vapor specific humidity [unitless]
    sh = qv / (1 + qv)
    
    # Saturation water vapor specific humidity
    # es = 379.90516 / pres * exp(17.2693882 * (temp - 273.16) / (temp - 35.86))
    es = pq0 / pres * math.exp(a2 * (temp - a3) / (temp - a4))
    rh = 100 * max(min(sh / es, 1), 0)
    
    print(f'{sh=}')
    print(f'{es=}')
    print(f'{rh=}')
    return rh

def cordex_rh(qv, pres, temp):
    # August-Roche-Magnus formula for saturated water vapor pressure
    
    # Saturated water vapor pressure [Pa]
    es = 610.94 * math.exp(17.625 * (temp - 273.15) / (temp - 30.11))
    
    # Saturated mixing ratio [kg kg-1]
    # NOTE: I'm not sure this isn't specific humidity (incomplete)
    ws = 0.622 * es / (pres - es)
    
    # original was rh = qv / (ws * 1000)
    rh = 100 * qv / ws
    
    print(f'{es=}')
    print(f'{ws=}')
    # print(f'{q=}')
    print(f'{rh=}')
    
    
def magnus_sat_vp(temp):
    # Magnus saturation vapor pressure formula
    # Relative error < 0.4% over -40C <= t <= 50C
    c1 = 610.94   # [Pa]
    a1 = 17.625
    b1 = 243.04   # [C]
    
    temp_c = temp - 273.15
    es = c1 * math.exp(a1 * temp_c / (b1 + temp_c))
    
    print(f'{es=}')


# %%

# %%
wrf_diag_rh(qv=0.001113043, pres=101543.3, temp=260.0828)

# %%
ncl_wrf_rh(qv=0.001113043, pres=101543.3, temp=260.0828)

# %%
teten_rh(qv=0.001113043, pres=101543.3, temp=260.0828)

# %%

# %%
wrf_rh(qv=0.001113043, pres=101543.3, temp=260.0828)

# %%
cordex_rh(qv=0.001113043, pres=101543.3, temp=260.0828)

# %%
qv=0.001113043
pres=101543.3
temp=260.0828

e_v1 = (qv * pres) / (0.622 + qv)   # teten eqn [Pa]

e_v2 = 6.113 * math.exp(17.2694 * (temp - 273.15) / (temp - 35.86))   # Teten's eqn [mb]
e_v2 *= 100   # Convert to [Pa]

e_v3 = qv / (1 + qv)   # wrf_diag eqn [unitless]; specific humidity

e_v4 = qv * pres / math.exp(17.67 * (temp - 273.15) / (temp - 29.65))   # [Pa]

# from wrf_rh [Pa]
e_v5 = pres * qv / (qv * (1.0 - 0.622) + 0.622)

print(f'{e_v1=}')
print(f'{e_v2=}')
print(f'{e_v3=}')
print(f'{e_v4=}')
print(f'{e_v5=}')


# %%
a1 = 0.622 + qv   # denominator of vapor pressure eqn (Teten)
a2 = 0.378 * qv + 0.622   # wrf_rh

print(f'{a1=}')
print(f'{a2=}')

# %%

# %%
magnus_sat_vp(temp)

# %%

# %%
c1 = 610.94   # [Pa]
a1 = 17.625
b1 = 243.04   # [C]
epsilon = 0.622  # Rd / Rv; []

# Vapor pressure
e = (qv * pres) / (epsilon + qv)   # [Pa]

# Dewpoint temperature
td = (b1 * math.log(e / c1)) / (a1 - math.log(e / c1))

print(f'{e=}')
print(f'{td=}')

# %%


def compute_dewpoint_temperature(temperature, vp, sat_vp):
    #  237.3 * X / ( 17.269 - X ) ;  where X = { ln ( E2 / ESAT2 ) + 17.269 * ( T2 - 273.15 ) / ( T2 - 35.85 ) }
    # equation from Milly's DRB spreadsheet
    x = np.log(vp / sat_vp) + 17.269 * (temperature - 273.15) / (temperature - 35.85)
    return 237.3 * x / (17.269 - x)

def teten_saturation_vp(temperature: Union[float, xr.Dataset, xr.DataArray]):
    """Saturation vapor pressure from temperature

    :param temperature: Temperature [K]
    """
    # Teten's eqn
    es0 = 6.113  # Saturation vapor pressure reference value; [hPa]
    es0 *= 100   # [Pa]

    temp_c = temperature - 273.15   # [C]

    # Saturation vapor pressure
    es = es0 * np.exp(17.269 * temp_c / (temp_c + 237.3))   # [Pa]
    return es

def compute_vp(vp_mixing_ratio: Union[float, xr.Dataset, xr.DataArray],
               pressure: Union[float, xr.Dataset, xr.DataArray]):
    """Water vapor pressure from mixing ratio and pressure

    :param vp_mixing_ratio: Vapor pressure mixing ratio [kg kg-1]
    :param pressure: Pressure [Pa]
    """
    epsilon = 0.622  # Rd / Rv; []

    e = vp_mixing_ratio * pressure / (epsilon + vp_mixing_ratio)
    return e


# %%
e = compute_vp(qv, pres)
es = teten_saturation_vp(temp)

td_milly = compute_dewpoint_temperature(temp, e, es)
print(f'{td_milly=}')

# %%

# %%

# %%

# %%

# %%

# %%
