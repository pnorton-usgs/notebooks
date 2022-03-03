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
import dask
import netCDF4 as nc
import numpy as np
import os
import pandas as pd
import xarray as xr

from collections import OrderedDict

# from pyPRMS.prms_helpers import dparse
from pyPRMS.ParameterFile import ParameterFile

# %%
# NHMv1.0
workdir = '/Volumes/USGS_NHM1/calibrations/NHMv10/DAYMET_releases/byHRU'
param_filename = f'{workdir}/NHM-PRMS.param'
filename = f'{workdir}/output/NHM-PRMS_data_release.csv'

poi_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/NHM_v1.0/poi_data'
poi_data = f'{poi_dir}/*_pois.nc'

# NHMv1.1
# workdir = '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHRU/20210624_calib'
# workdir = '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHRU/20210617_gm_PRECAL'
# param_filename = f'{workdir}/myparam.param'
# filename = f'{workdir}/output_variables/stats.csv'
# outfilename = f'{workdir}/output_variables/streamflow.nc'

# poi_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/poi_data'
# poi_data = f'{poi_dir}/*_pois.nc'

# %% [markdown]
# ## Read the model output streamflow
# POIs and segment info is pulled from the header

# %%
def read_csv_header(filename):
    fhdl = open(filename, 'r')

    # First and second rows are headers
    hdr1 = fhdl.readline().strip()

    fhdl.close()
    
    tmp_flds = hdr1.split(' ')
    tmp_flds.remove('Date')

    flds = {nn+3: hh for nn, hh in enumerate(tmp_flds)}

    # poi_flds maps column index to POI and is used to rename the dataframe columns from indices to station IDs
    poi_flds = OrderedDict()

    # poi_seg_flds maps POI to the related segment ID
    poi_seg_flds = OrderedDict()

    for xx, yy in flds.items():
        tfld = yy.split('_')
        segid = int(tfld[2]) - 1  # Change to zero-based indices
        poiid = tfld[4]

        poi_flds[xx] = poiid
        poi_seg_flds[poiid] = segid

    return poi_flds, poi_seg_flds



# %%
def read_streamflow(filename, field_names):
    df = pd.read_csv(filename, sep='\s+', header=None, skiprows=2, parse_dates={'time': [0, 1, 2]}, 
                     index_col='time')

    df.rename(columns=field_names, inplace=True)
    
    return df


# %%

poi_flds, poi_seg_flds = read_csv_header(filename)
model_pois_segments = list(poi_seg_flds.values())

df = read_streamflow(filename, field_names=poi_flds)

# %%
df.info()

# %%

# %%
# Need the list of POIs and date range from the model streamflow output file
model_pois = df.columns.tolist()

st_date = min(df.index.tolist())
en_date = max(df.index.tolist())

print(f'Number of POIs: {len(model_pois)}')
print(f'Start date: {st_date}')
print(f'End date: {en_date}')

# %%

# %% [markdown]
# ## Read the model parameter file
# The parameters needed from this file are poi_gage_id and seg_cum_area. The seg_cum_area parameter is converted to square miles.

# %%
pfile = ParameterFile(param_filename, verbose=True, verify=True)

# Get the list of POI IDs in the parameter file
poi_gage_id = pfile.parameters['poi_gage_id']
poi_gage_id_list = poi_gage_id.data.tolist()

# Get the cumulative drainage area for each of the POI segments
seg_cum_area = pfile.parameters.get_dataframe('seg_cum_area')

# using iloc because model_pois_segments are 0-based local indices
model_pois_area = seg_cum_area.iloc[model_pois_segments].to_numpy(dtype=float)

# Convert model area from acres to square miles
model_pois_area /= 640.

# %%

# %% [markdown]
# ## Verify the order of POIs in the output file and the parameter file are the same

# %%
# Check that the parameter file poi order and output file order are the same
for xx, yy in poi_flds.items():
    if yy != poi_gage_id_list[xx-3]:
        print(f'Out of order: {yy}, {poi_gage_id_list[xx-3]}')

# %%

# %% [markdown]
# ## Read the streamflow observation file(s)
# The streamflow observations are read from the local netcdf files used by Bandit for model extractions.

# %%
# Open the streamflow observation files
poi_df = xr.open_mfdataset(poi_data, chunks={'poi_id': 2000}, combine='nested',
                           concat_dim='poi_id', decode_cf=True, engine='netcdf4')
poi_df

# %%
poi_obs_crap = poi_df['poi_id'].to_pandas()

# %%
poi_obs_list = poi_obs_crap.tolist()

# %%
# set(model_pois) - set(poi_obs_list)
set(poi_obs_list) - set(model_pois)


# %%

# %%

# %%

# %%

# %%

# %%
# Helper function to pull variables from the streamflow observation files
def get_var(var, df, gageids, stdate=None, endate=None):
    if stdate is None or endate is None:
        data = df[var].loc[gageids].to_pandas()
    else:
        data = df[var].loc[gageids, stdate:endate].to_pandas()
    return data


# %%

# %%
# This does a bunch of selects using an array of POI IDs which causes a warning from xarray
# related to performance. The following silences the warning.
# See: https://docs.dask.org/en/latest/array-slicing.html
# dask.config.set({"array.slicing.split_large_chunks": False})

# poiname_list = self.get('poi_name').tolist()
poiname_list = get_var('poi_name', poi_df, model_pois).tolist()

max_poiid_len = len(max(model_pois, key=len))
max_poiname_len = len(max(poiname_list, key=len))

# Create a netCDF file for the CBH data
nco = nc.Dataset(outfilename, 'w', clobber=True)

# Create the dimensions
nco.createDimension('poiid_nchars', max_poiid_len)
nco.createDimension('poi_id', len(model_pois))
nco.createDimension('poiname_nchars', max_poiname_len)
nco.createDimension('time', None)

reference_time = st_date.strftime('%Y-%m-%d %H:%M:%S')
cal_type = 'standard'

# Create the variables
timeo = nco.createVariable('time', 'f4', 'time')
timeo.calendar = cal_type
timeo.units = f'days since {reference_time}'

poiido = nco.createVariable('poi_id', 'S1', ('poi_id', 'poiid_nchars'), zlib=True)
poiido.long_name = 'Point-of-Interest ID'
poiido.cf_role = 'timeseries_id'
poiido._Encoding = 'ascii'

poinameo = nco.createVariable('poi_name', 'S1', ('poi_id', 'poiname_nchars'), zlib=True)
poinameo.long_name = 'Name of POI station'

lato = nco.createVariable('latitude', 'f4', 'poi_id', zlib=True)
lato.long_name = 'Latitude'
lato.units = 'degrees_north'

lono = nco.createVariable('longitude', 'f4', 'poi_id', zlib=True)
lono.long_name = 'Longitude'
lono.units = 'degrees_east'

draino = nco.createVariable('drainage_area', 'f4', 'poi_id',
                            fill_value=nc.default_fillvals['f4'], zlib=True)
draino.long_name = 'Drainage Area'
draino.units = 'mi2'

draineffo = nco.createVariable('drainage_area_contrib', 'f4', 'poi_id',
                               fill_value=nc.default_fillvals['f4'], zlib=True)
draineffo.long_name = 'Effective drainage area'
draineffo.units = 'mi2'

model_da_o = nco.createVariable('drainage_area_model', 'f4', 'poi_id', 
                                 fill_value=nc.default_fillvals['f4'], zlib=True)
model_da_o.long_name = 'Model drainage area'
model_da_o.units = 'mi2'


varo = nco.createVariable('discharge', 'f4', ('poi_id', 'time'),
                          fill_value=nc.default_fillvals['f4'], zlib=True)
varo.long_name = 'discharge'
varo.units = 'ft3 s-1'

model_varo = nco.createVariable('discharge_model', 'f4', ('poi_id', 'time'),
                          fill_value=nc.default_fillvals['f4'], zlib=True)
model_varo.long_name = 'Model discharge'
model_varo.units = 'ft3 s-1'

nco.setncattr('Description', 'POI model output for PRMS')
nco.setncattr('FeatureType', 'timeSeries')
# nco.setncattr('Bandit_version', __version__)
# nco.setncattr('NHM_version', nhmparamdb_revision)

# Write the Streamgage IDs
poiido[:] = nc.stringtochar(np.array(model_pois).astype('S'))

timeo[:] = nc.date2num(pd.to_datetime(df.index).tolist(),
                       units=f'days since {reference_time}',
                       calendar=cal_type)

# Write the streamgage observations
varo[:, :] = get_var('discharge', poi_df, model_pois, st_date, en_date).to_numpy(dtype=float)

# Write the simulated streamflow
model_varo[:, :] = df.to_numpy(dtype=float).T

poinameo[:] = nc.stringtochar(np.array(poiname_list).astype('S'))
lato[:] = get_var('latitude', poi_df, model_pois).to_numpy(dtype=float)
lono[:] = get_var('longitude', poi_df, model_pois).to_numpy(dtype=float)
draino[:] = get_var('drainage_area', poi_df, model_pois).to_numpy(dtype=float)
draineffo[:] = get_var('drainage_area_contrib', poi_df, model_pois).to_numpy(dtype=float)
model_da_o[:] = model_pois_area
nco.close()

# %%
nco.close()

# %%

# %%
