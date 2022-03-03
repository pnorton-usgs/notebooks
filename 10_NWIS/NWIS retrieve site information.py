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
import datetime
import numpy as np
import pandas as pd
import re
from io import StringIO
import sys

from urllib.request import urlopen, Request
from urllib.error import HTTPError

# %%
workdir = '/Users/pnorton/tmp/streamflow_work'
stnfile = f'{workdir}/something'

base_url = 'http://waterservices.usgs.gov/nwis'

t1 = re.compile('^#.*$\n?', re.MULTILINE)   # remove comment lines
t2 = re.compile('^5s.*$\n?', re.MULTILINE)  # remove field length lines


# %%

def nwis_site_fields():
    # Retrieve a single station and pull out the field names and data types
    stn_url = f'{base_url}/site/?format=rdb&sites=01646500&siteOutput=expanded&siteStatus=active&parameterCd=00060&siteType=ST'

    response = urlopen(stn_url)
    encoding = response.info().get_param('charset', failobj='utf8')
    streamgage_site_page = response.read().decode(encoding)

    # Strip the comment lines and field length lines from the result
    streamgage_site_page = t1.sub('', streamgage_site_page, 0)

    nwis_dtypes = t2.findall(streamgage_site_page)[0].strip('\n').split('\t')
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


# %%
# Get fields and datatypes
cols = nwis_site_fields()

# Columns to include in the final dataframe
include_cols = ['agency_cd', 'site_no', 'station_nm', 'dec_lat_va', 'dec_long_va', 'dec_coord_datum_cd',
                'alt_va', 'alt_datum_cd', 'huc_cd', 'drain_area_va', 'contrib_drain_area_va']

# Start with an empty dataframe
nwis_sites = pd.DataFrame(columns=include_cols)
# nwis_sites.info()

# %%
# Retrieve the station information by HUC2

# URLs can be generated/tested at: http://waterservices.usgs.gov/rest/Site-Test-Tool.html
base_url = 'http://waterservices.usgs.gov/nwis'

for region in range(19):
    # region = '01'
    print(f'Region {region+1:02}')
    stn_url = f'{base_url}/site/?format=rdb&huc={region+1:02}&siteOutput=expanded&siteStatus=all&parameterCd=00060&siteType=ST'

    response = urlopen(stn_url)
    encoding = response.info().get_param('charset', failobj='utf8')
    streamgage_site_page = response.read().decode(encoding)

    # Strip the comment lines and field length lines from the result
    streamgage_site_page = t1.sub('', streamgage_site_page, 0)
    streamgage_site_page = t2.sub('', streamgage_site_page, 0)

    # Read the rdb file into a dataframe
    df = pd.read_csv(StringIO(streamgage_site_page), sep='\t', dtype=cols, usecols=include_cols)

    nwis_sites = nwis_sites.append(df, ignore_index=True)

# %%
nwis_sites.info()

# %%

# %%
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
nwis_sites.to_csv('/Users/pnorton/Projects/National_Hydrology_Model/datasets/streamflow/nwis_sites_simple.tab', 
                  sep='\t')

# %%

# %%

# %%

# %%
miss_site = ['11403800', '12158010', '12157250', '04095380', '12158040', '11066460', '06026420',
             '11173500', '06204070', '12202300', '11230200', '06354490', '12113347', '11292700',
             '06307990', '01180000', '06175520', '06024020', '06288400', '05051500', '09343300',
             '06036805', '11065000', '06023100', '11192950', '06027600', '12143700', '06287800',
             '12202420', '11316700']
nwis_missing = nwis_sites[nwis_sites.index.isin(miss_site)]
print(f'Number of missing sites: {len(miss_site)}')
print(nwis_missing)

# %%
set(miss_site) - set(nwis_missing.index.tolist())

# %%
all_sites = nwis_sites.index.tolist()

# %%
nwis_sites.loc['01010500']['drainage_area']

# %% [markdown]
# ## Check daily completeness

# %%
daily_filename = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/streamflow/conus_daily_HUC_03_obs.tab'

# %%
col_names = ['agency_cd', 'site_no', 'parameter_cd', 'ts_id', 'loc_web_ds', 'month_nu', 'day_nu', 'begin_yr',
             'end_yr', 'count_nu', 'mean_va']
col_types = [np.str_, np.str_, np.str_, np.int, np.str_, np.int, np.int, np.int, np.int, np.int, np.float]
cols = dict(zip(col_names, col_types))

# %%
daily_df = pd.read_csv(daily_filename, sep='\t', dtype=cols)
daily_df.info()


# %%
def nwis_load_daily_statistics(src_dir):
    col_names = ['agency_cd', 'site_no', 'parameter_cd', 'ts_id', 'loc_web_ds', 'month_nu', 'day_nu',
                 'begin_yr', 'end_yr', 'count_nu', 'mean_va']
    col_types = [np.str_, np.str_, np.str_, np.int, np.str_, np.int, np.int, np.int, np.int, np.int, np.float]
    cols = dict(zip(col_names, col_types))

    # Start with an empty dataframe
    nwis_daily = pd.DataFrame(columns=col_names)

    for region in range(18):
        # region = '01'
        sys.stdout.write(f'\rRegion: {region+1:02}')
        sys.stdout.flush()
        # print(f'Region {region+1:02}')

        # Read the rdb file into a dataframe
        df = pd.read_csv(f'{src_dir}/conus_daily_HUC_{region+1:02}_obs.tab', sep='\t', dtype=cols)

        nwis_daily = nwis_daily.append(df, ignore_index=True)
    print('')
    return nwis_daily


# %%
daily_df = nwis_load_daily_statistics('/Users/pnorton/Projects/National_Hydrology_Model/datasets/streamflow')

# %%
daily_df[daily_df['site_no'] == '01066000']

# %%
# nwis_cnt = daily_df.groupby(['site_no'])['loc_web_ds'].count()
nwis_cnt = daily_df.groupby(['site_no'])['loc_web_ds'].nunique(dropna=False)

# nwis_cnt[nwis_cnt < 2]
print(nwis_cnt)

# %%
nwis_cnt = nwis_cnt[nwis_cnt[]]

# %%
nwis_cnt['05051500']

# %%
daily_sum[daily_sum == por_days.days]

# %%
daily_sum['09474000']

# %%
nwis_multiple = daily_df[daily_df['loc_web_ds'].notnull()]
aa = nwis_multiple['site_no'].tolist()
print(set(aa))

print(daily_df[~daily_df['site_no'].isin(aa)])

# %%
daily_df[daily_df['loc_web_ds'].notnull()]

# %%

# %%

# %%
nwis_sites.loc['05051500', ]


# %%

# %%
def nwis_site_param_cd(poi_id):
    col_names = ['agency_cd', 'site_no', 'station_nm', 'site_tp_cd', 'dec_lat_va', 'dec_long_va', 
                 'coord_acy_cd', 'dec_coord_datum_cd', 'alt_va', 'alt_acy_va', 'alt_datum_cd', 'huc_cd', 
                 'data_type_cd', 'parm_cd', 'stat_cd', 'ts_id', 'loc_web_ds', 'medium_grp_cd', 
                 'parm_grp_cd', 'srs_id', 'access_cd', 'begin_date', 'end_date', 'count_nu']
    col_types = [np.str_, np.str_, np.str_, np.str_, np.float, np.float,
                 np.str_, np.str_, np.float, np.float, np.str_, np.str_,
                 np.str_, np.str_, np.str_, np.int, np.str_, np.str_,
                 np.str_, np.int, np.str_, np.str_, np.str_, np.int]
    cols = dict(zip(col_names, col_types))

    # Retrieve a single station and pull out the field names and data types
    stn_url = f'{base_url}/site/?format=rdb&sites={poi_id}&seriesCatalogOutput=true&siteStatus=all'

    response = urlopen(stn_url)
    encoding = response.info().get_param('charset', failobj='utf8')
    streamgage_site_page = response.read().decode(encoding)

    # Strip the comment lines and field length lines from the result
    streamgage_site_page = t1.sub('', streamgage_site_page, 0)
    streamgage_site_page = t2.sub('', streamgage_site_page, 0)

    df = pd.read_csv(StringIO(streamgage_site_page), sep='\t', dtype=cols, usecols=['site_no', 'parm_cd'])
    
    return df['parm_cd'].tolist()


# %%
aa = nwis_site_param_cd('11383730')
print(aa)
print('00060' in aa)
print('70331' in aa)

# %%

# %%

# %% [markdown]
# ## Get station period-of-record information from NWIS

# %%
# Retrieve the station period of record information by HUC2
col_names = ['agency_cd', 'site_no', 'station_nm', 'site_tp_cd', 'dec_lat_va', 'dec_long_va', 
             'coord_acy_cd', 'dec_coord_datum_cd', 'alt_va', 'alt_acy_va', 'alt_datum_cd', 'huc_cd', 
             'data_type_cd', 'parm_cd', 'stat_cd', 'ts_id', 'loc_web_ds', 'medium_grp_cd', 
             'parm_grp_cd', 'srs_id', 'access_cd', 'begin_date', 'end_date', 'count_nu']
col_types = [np.str_, np.str_, np.str_, np.str_, np.float, np.float,
             np.str_, np.str_, np.float, np.float, np.str_, np.str_,
             np.str_, np.str_, np.str_, np.int, np.str_, np.str_,
             np.str_, np.int, np.str_, np.str_, np.str_, np.int]
cols = dict(zip(col_names, col_types))

nwis_sites_exp = pd.DataFrame(columns=col_names)


# URLs can be generated/tested at: http://waterservices.usgs.gov/rest/Site-Test-Tool.html
base_url = 'http://waterservices.usgs.gov/nwis'

for region in range(19):
    # region = '01'
    print(f'Region {region+1:02}')
    stn_url = f'{base_url}/site/?format=rdb&huc={region+1:02}&seriesCatalogOutput=true&outputDataTypeCd=dv&siteStatus=all'

    response = urlopen(stn_url)
    encoding = response.info().get_param('charset', failobj='utf8')
    streamgage_site_page = response.read().decode(encoding)

    # Strip the comment lines and field length lines from the result
    streamgage_site_page = t1.sub('', streamgage_site_page, 0)
    streamgage_site_page = t2.sub('', streamgage_site_page, 0)

    # Read the rdb file into a dataframe
    df = pd.read_csv(StringIO(streamgage_site_page), sep='\t', dtype=cols, usecols=col_names)

    nwis_sites_exp = nwis_sites_exp.append(df, ignore_index=True)
    
# Coerce some of the field datatypes
nwis_sites_exp['count_nu'] = pd.to_numeric(nwis_sites_exp['count_nu'], errors='raise', downcast='integer')
nwis_sites_exp['ts_id'] = pd.to_numeric(nwis_sites_exp['ts_id'], errors='raise', downcast='integer')
nwis_sites_exp['srs_id'] = pd.to_numeric(nwis_sites_exp['srs_id'], errors='raise', downcast='integer')

# %%
nwis_sites_exp.to_csv('/Users/pnorton/Projects/National_Hydrology_Model/datasets/streamflow/nwis_sites_por.tab', sep='\t',
                      columns=col_names)

# %% [markdown]
# ## Read station period-of-record information from a file

# %%
col_names = ['agency_cd', 'site_no', 'station_nm', 'site_tp_cd', 'dec_lat_va', 'dec_long_va',
             'coord_acy_cd', 'dec_coord_datum_cd', 'alt_va', 'alt_acy_va', 'alt_datum_cd', 'huc_cd',
             'data_type_cd', 'parm_cd', 'stat_cd', 'ts_id', 'loc_web_ds', 'medium_grp_cd',
             'parm_grp_cd', 'srs_id', 'access_cd', 'begin_date', 'end_date', 'count_nu']
col_types = [np.str_, np.str_, np.str_, np.str_, np.float, np.float,
             np.str_, np.str_, np.float, np.float, np.str_, np.str_,
             np.str_, np.str_, np.str_, np.int, np.str_, np.str_,
             np.str_, np.int, np.str_, np.str_, np.str_, np.int]
cols = dict(zip(col_names, col_types))

nwis_sites_exp = pd.read_csv('/Users/pnorton/Projects/National_Hydrology_Model/datasets/streamflow/nwis_sites_por.tab',
                             sep='\t', dtype=cols, index_col=0)

# %%
nwis_sites_exp.info()

# %%
nwis_sites_exp[(nwis_sites_exp['site_no']=='05051500') & (nwis_sites_exp['parm_cd']=='00060') &
               (nwis_sites_exp['stat_cd']=='00003')][['begin_date', 'end_date']].to_numpy()[0].tolist()


# %%
def has_multiple_timeseries(df, poi_id, param_cd, stat_cd):
    return df[(df['site_no'] == poi_id) & (df['parm_cd'] == param_cd) & 
              (df['stat_cd'] == stat_cd)].shape[0] > 1
             #(df['stat_cd'] == stat_cd) & (df['loc_web_ds'].notnull())].shape[0] > 1


# %%
has_multiple_timeseries(nwis_sites_exp, '01019000', '00060', '00003')


# %%

# %%
def has_param_statistic(df, poi_id, param_cd, stat_cd):
    # Returns true if poi_id  has a single entry for the given param_cd and stat_cd
    return df[(df['site_no'] == poi_id) & (df['parm_cd'] == param_cd) & (df['stat_cd'] == stat_cd)].shape[0] == 1


# %%
has_param_statistic(nwis_sites_exp, '07252406', '00060', '00003')


# %%
def nwis_site_type(poi_id):
    # Retrieve a single station and pull out the field names and data types
    stn_url = f'{base_url}/site/?format=rdb&sites={poi_id}&siteStatus=all'

    response = urlopen(stn_url)
    encoding = response.info().get_param('charset', failobj='utf8')
    streamgage_site_page = response.read().decode(encoding)

    # Strip the comment lines and field length lines from the result
    streamgage_site_page = t1.sub('', streamgage_site_page, 0)
    streamgage_site_page = t2.sub('', streamgage_site_page, 0)

    df = pd.read_csv(StringIO(streamgage_site_page), sep='\t', header=0) #, usecols=['site_type_cd'])

    return df['site_tp_cd'].tolist()


# %%
nwis_site_type('07176460')[0]

# %%
poi_id = '07176460'

# Retrieve a single station and pull out the field names and data types
stn_url = f'{base_url}/site/?format=rdb&sites={poi_id}&siteStatus=all'

response = urlopen(stn_url)
encoding = response.info().get_param('charset', failobj='utf8')
streamgage_site_page = response.read().decode(encoding)

# Strip the comment lines and field length lines from the result
streamgage_site_page = t1.sub('', streamgage_site_page, 0)
streamgage_site_page = t2.sub('', streamgage_site_page, 0)

df = pd.read_csv(StringIO(streamgage_site_page), sep='\t', header=0) #, usecols=['site_type_cd'])

# %%
df

# %%
daily_df[daily_df['site_no'] == '02160000']

# %%
