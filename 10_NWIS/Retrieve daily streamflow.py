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
import socket
import sys

from io import StringIO

from urllib.request import urlopen, Request
from urllib.error import HTTPError

# %%
workdir = '/Users/pnorton/tmp/streamflow_work'
stnfile = f'{workdir}/something'

base_url = 'http://waterservices.usgs.gov/nwis'

t1 = re.compile('^#.*$\n?', re.MULTILINE)   # remove comment lines
t2 = re.compile('^5s.*$\n?', re.MULTILINE)  # remove field length lines

# %%
# Set timeout in seconds - if not set defaults to infinite time for response
timeout = 30
socket.setdefaulttimeout(timeout)


# %%
# class WaterServices:
#     def __init__(self):
#         pass
    
#     def get():
#         # Get something from water services
#         response = urlopen(f'{BASE_NWIS_URL}/dv/{url_final}')
#         encoding = response.info().get_param('charset', failobj='utf8')
#         streamgage_obs_page = response.read().decode(encoding)

# %%
def nwis_load_daily_statistics(src_dir):
    col_names = ['agency_cd', 'site_no', 'parameter_cd', 'ts_id', 'loc_web_ds', 'month_nu', 'day_nu',
                 'begin_yr', 'end_yr', 'count_nu', 'mean_va']
    col_types = [np.str_, np.str_, np.str_, np.int, np.str_, np.int, np.int, np.int, np.int, np.int, np.float]
    cols = dict(zip(col_names, col_types))

    # Start with an empty dataframe
    nwis_daily = pd.DataFrame(columns=col_names)

    for region in range(1, 19):
        # region = '01'
        sys.stdout.write(f'\rRegion: {region:02}')
        sys.stdout.flush()
        # print(f'Region {region+1:02}')

        # Read the rdb file into a dataframe
        df = pd.read_csv(f'{src_dir}/conus_daily_HUC_{region:02}_obs.tab', sep='\t', dtype=cols)

        nwis_daily = nwis_daily.append(df, ignore_index=True)
    print('')
    return nwis_daily


# %%
nwis_daily = nwis_load_daily_statistics('/Users/pnorton/Projects/National_Hydrology_Model/datasets/streamflow')

# %%
nwis_cnt = nwis_daily.groupby(['site_no'])['loc_web_ds'].nunique()
# nwis_cnt = nwis_cnt[nwis_cnt < 2]

# %%
nwis_cnt['07311630']

# %%
nwis_daily[nwis_daily['site_no'] == '07311630']

# %%
nwis_daily

# %%
