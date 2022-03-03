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
#     display_name: Python [conda env:bandit]
#     language: python
#     name: conda-env-bandit-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
from __future__ import (absolute_import, division, print_function)
from future.utils import iteritems, iterkeys

import pandas as pd
import msgpack
import os

from pyPRMS.constants import REGIONS
from pyPRMS.prms_helpers import dparse
from pyPRMS.Cbh import Cbh


# %%
def get_parameter(filename):
    with open(filename, 'rb') as ff:
        return msgpack.load(ff, use_list=False)


# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/daymet'
merged_paramdb_dir = '/Users/pnorton/Projects/National_Hydrology_Model/paramDb/merged_params_3'
var = 'prcp'
# r17_prcp.cbh.gz

CBH_VARNAMES = ['prcp', 'tmin', 'tmax']
CBH_INDEX_COLS = [0, 1, 2, 3, 4, 5]

# %%
first = True

hru_nhm_to_local = get_parameter('{}/hru_nhm_to_local.msgpack'.format(merged_paramdb_dir))
hru_nhm_to_region = get_parameter('{}/hru_nhm_to_region.msgpack'.format(merged_paramdb_dir))

for rr in REGIONS:
# for rr, rvals in iteritems(hru_nhm_to_region):
    rvals = hru_nhm_to_region[rr]
    print('Examining {} ({} to {})'.format(rr, rvals[0], rvals[1]))
    if rvals[0] >= rvals[1]:
        raise ValueError('Lower HRU bound is greater than upper HRU bound.')

    idx_retrieve = {}

    for yy in hru_nhm_to_local.keys():
        if rvals[0] <= yy <= rvals[1]:
            # print('\tMatching region {}, HRU: {} ({})'.format(rr, yy, hru_order_ss[yy]))
            idx_retrieve[yy] = hru_nhm_to_local[yy]

    if len(idx_retrieve) > 0:
        # The current region contains HRUs in the model subset
        # Read in the data for those HRUs
        cbh_file = '{}/{}_{}.cbh.gz'.format(workdir, rr, var)

        if not os.path.isfile(cbh_file):
            # Missing data file for this variable and region
            bandit_log.error('Required CBH file, {}, is missing. Unable to continue'.format(cbh_file))
            raise IOError('Required CBH file, {}, is missing.'.format(cbh_file))

        cc1 = Cbh(filename=cbh_file, indices=idx_retrieve)
        cc1.read_cbh()

        if first:
            outdata = cc1.data.copy()
            first = False
        else:
            outdata = pd.merge(outdata, cc1.data, on=[0, 1, 2, 3, 4, 5])

# %%
outdata.info()

# %%
outdata.head()

# %%
outdata.rename(columns={0:'yy', 1:'mm', 2:'dd', 3:'hr', 4:'min', 5:'sec'}, inplace=True)

# %%
import feather

# %%
feather.write_dataframe(outdata, '{}/daymet_{}.feather'.format(workdir, var))

# %%
import pyarrow as pa
import pyarrow.parquet as pq

# %%
arrow_table = pa.Table.from_pandas(outdata)
pq.write_table(arrow_table, '{}/daymet_{}.pq'.format(workdir, var), use_dictionary=True, compression='snappy')

# %%
# Read a pyarrow table into a dataframe
ndf = pq.read_table('{}/daymet_{}.pq'.format(workdir, var), nthreads=1).to_pandas()

# %%
