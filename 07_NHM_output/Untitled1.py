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
# import dask
# import netCDF4 as nc
import datetime
import numpy as np
import os
import pandas as pd
import re
import xarray as xr
import sys

from io import StringIO

from urllib.request import urlopen  # , Request
from urllib.error import HTTPError

from collections import OrderedDict

# from pyPRMS.prms_helpers import dparse
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ParamDb import ParamDb

# %%
paramdb_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/paramdb_v10_daymet_CONUS'

# workdir = '/Volumes/USGS_NHM1/calibrations/NHMv11/gridmet_byHRU/20210617_gm_PRECAL'
# param_filename = f'{workdir}/myparam.param'

poi_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/poi_data'
poi_data = f'{poi_dir}/NWIS_pois.nc'

st_date = datetime.datetime(1940, 10, 1)
en_date = datetime.datetime(2020, 9, 30)

# %%
pdb = ParamDb(paramdb_dir=paramdb_dir, verbose=True, verify=True)

# %%
# Get the POI gages from the parameter database
poi_gage_id = pdb.parameters['poi_gage_id']
poi_gage_id_list = poi_gage_id.data.tolist()

# %%
len(poi_gage_id_list)

# %%

# %% [markdown]
# ## Read the observed streamflow values

# %%
obs_df = xr.open_mfdataset(poi_data, chunks={'poi_id': 2000}, combine='nested',
                                    concat_dim='poi_id', decode_cf=True, engine='netcdf4')

# %%
obs_df

# %%
obs_poi_id = obs_df['poi_id'].to_pandas()

# %%
obs_poi_id

# %%
obs_poi_id_list = obs_poi_id.tolist()

# %%
missing_from_obs = set(poi_gage_id_list) - set(obs_poi_id_list)

# %%
len(missing_from_obs)

# %%
missing_from_obs

# %%
# URLs can be generated/tested at: http://waterservices.usgs.gov/rest/Site-Test-Tool.html
base_url = 'http://waterservices.usgs.gov/nwis'

t1 = re.compile('^#.*$\n?', re.MULTILINE)   # remove comment lines
t2 = re.compile('^5s.*$\n?', re.MULTILINE) 

def _retrieve_from_nwis(url):
    response = urlopen(url)
    encoding = response.info().get_param('charset', failobj='utf8')
    streamgage_site_page = response.read().decode(encoding)

    # Strip the comment lines and field length lines from the result
    streamgage_site_page = t1.sub('', streamgage_site_page, 0)
    streamgage_site_page = t2.sub('', streamgage_site_page, 0)

    return streamgage_site_page


def get_nwis_site_fields():
    # Retrieve a single station and pull out the field names and data types

    stn_url = f'{base_url}/site/?format=rdb&sites=01646500&siteOutput=expanded&' + \
              'siteStatus=active&parameterCd=00060&siteType=ST'

    response = urlopen(stn_url)
    encoding = response.info().get_param('charset', failobj='utf8')
    streamgage_site_page = response.read().decode(encoding)

    # Strip the comment lines and field length lines from the result
    streamgage_site_page = t1.sub('', streamgage_site_page, 0)

    # nwis_dtypes = t2.findall(streamgage_site_page)[0].strip('\n').split('\t')
    nwis_fields = StringIO(streamgage_site_page).getvalue().split('\n')[0].split('\t')

    nwis_final = {}
    for fld in nwis_fields:
        code = fld[-2:]
        if code in ['cd', 'no', 'nm', 'dt']:
            nwis_final[fld] = np.str_
        elif code in ['va']:
            nwis_final[fld] = np.float32
        else:
            nwis_final[fld] = np.str_

    return nwis_final


def get_nwis_sites(stdate, endate, sites=None, regions=None):
    cols = get_nwis_site_fields()

    # Columns to include in the final dataframe
    include_cols = ['agency_cd', 'site_no', 'station_nm', 'dec_lat_va', 'dec_long_va', 'dec_coord_datum_cd',
                    'alt_va', 'alt_datum_cd', 'huc_cd', 'drain_area_va', 'contrib_drain_area_va']

    # Start with an empty dataframe
    nwis_sites = pd.DataFrame(columns=include_cols)

    url_pieces = OrderedDict()
    url_pieces['format'] = 'rdb'
    url_pieces['startDT'] = stdate.strftime('%Y-%m-%d')
    url_pieces['endDT'] = endate.strftime('%Y-%m-%d')
    # url_pieces['huc'] = None
    url_pieces['siteOutput'] = 'expanded'
    url_pieces['siteStatus'] = 'all'
    # url_pieces['parameterCd'] = '00060'  # Discharge
    # url_pieces['siteType'] = 'ST'
    url_pieces['hasDataTypeCd'] = 'dv'
    url_pieces['access'] = 3

    # NOTE: If both sites and regions parameters are specified the sites
    #       parameter takes precedence.
    if sites is None:
        # No sites specified so default to HUC02-based retrieval
        url_pieces['huc'] = None

        if regions is None:
            # Default to HUC02 regions 1 thru 18
            regions = list(range(1, 19))
        if isinstance(regions, (list, tuple)):
            pass
        else:
            # Single region
            regions = [regions]
    else:
        # One or more sites with specified
        url_pieces['sites'] = None

        if isinstance(sites, (list, tuple)):
            pass
        else:
            # Single string, convert to list of sites
            sites = [sites]

    if 'huc' in url_pieces:
        # for region in range(19):
        for region in regions:
            sys.stdout.write(f'\r  Region: {region:02}')
            sys.stdout.flush()

            url_pieces['huc'] = f'{region:02}'
            url_final = '&'.join([f'{kk}={vv}' for kk, vv in url_pieces.items()])
            stn_url = f'{base_url}/site/?{url_final}'

            streamgage_site_page = _retrieve_from_nwis(stn_url)

            # Read the rdb file into a dataframe
            df = pd.read_csv(StringIO(streamgage_site_page), sep='\t', dtype=cols, usecols=include_cols)

            nwis_sites = nwis_sites.append(df, ignore_index=True)
            sys.stdout.write('\r                      \r')
    else:
        for site in sites:
            sys.stdout.write(f'\r  Site: {site} ')
            sys.stdout.flush()

            url_pieces['sites'] = site
            url_final = '&'.join([f'{kk}={vv}' for kk, vv in url_pieces.items()])

            stn_url = f'{base_url}/site/?{url_final}'

            try:
                streamgage_site_page = _retrieve_from_nwis(stn_url)

                # Read the rdb file into a dataframe
                df = pd.read_csv(StringIO(streamgage_site_page), sep='\t', dtype=cols, usecols=include_cols)

                nwis_sites = nwis_sites.append(df, ignore_index=True)
            except HTTPError as err:
                if err.code == 404:
                    sys.stdout.write(f'HTTPError: {err.code}, site does not meet criteria - SKIPPED\n')
            sys.stdout.write('\r                      \r')

    field_map = {'agency_cd': 'poi_agency',
                 'site_no': 'poi_id',
                 'station_nm': 'poi_name',
                 'dec_lat_va': 'latitude',
                 'dec_long_va': 'longitude',
                 'alt_va': 'elevation',
                 'drain_area_va': 'drainage_area',
                 'contrib_drain_area_va': 'drainage_area_contrib'}

    nwis_sites.rename(columns=field_map, inplace=True)
    nwis_sites.set_index('poi_id', inplace=True)
    nwis_sites = nwis_sites.sort_index()

    return nwis_sites


# %%
nwis_sites = get_nwis_sites(stdate=st_date, endate=en_date)

# %%
print(f'before: {len(nwis_sites)}')
nwis_sites_sub = nwis_sites[nwis_sites.index.isin(poi_gage_id_list)]
print(f'after: {len(nwis_sites_sub)}')

# %%
nwis_sites.head()

# %%
nwis_sites.loc['10309101'] = ['USGS', 'INVALID poi_id', 38.9459, -119.7796, 'NAD83', 0, 'NGVD29', '16050201', 0.0, 0.0]

# %%
nwis_sites.tail()

# %%
nwis_sites_list = nwis_sites_sub.index.tolist()
# nwis_sites_list
set(poi_gage_id_list) - set(nwis_sites_list)

# %%
# nwis_sites.loc['01111200']

# %%
aa = get_nwis_sites(stdate=st_date, endate=en_date, sites=poi_gage_id_list)

# %%

# %%
cdf = pd.read_csv('/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/NHM_v1.0/crap_error', sep ='\t')

# %%
cdf.info()

# %%
aa = cdf['contrib_drain_area_va'].tolist()

# %%
for ii, xx in enumerate(aa):
    print(ii, float(xx))

# %%
cdf.iloc[3961]

# %%

# %%
