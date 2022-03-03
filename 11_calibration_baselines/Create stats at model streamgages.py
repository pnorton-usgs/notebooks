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
import numpy as np
# import netCDF4 as nc
import pandas as pd

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/code_lauren/NHM_CONUS_CALIBRATION/TandM_FIGURES'
filename = f'{workdir}/StreamGageInfo.csv'

calibdir = '/Volumes/USGS_NHM1/calibrations/NHMv10/DAYMET_releases/byHRU'
rel_filename = f'{calibdir}/NHM-PRMS_data_release.csv'

# %%
# fld_names = ['GAGEid', 'Ref_Non', 'DI', 'Latitude', 'Longitude', 
#              'NWISda', 'GFda', 'GageName', 'Comment1', 'Comment2']
# fld_type = [np.str_, np.str_, np.int_, np.float_, np.float_,
#             np.float_, np.float_, np.str_, np.str_, np.str_]
fld_names = ['GAGEid', 'Ref_Non', 'DI', 'GFda']
fld_type = [np.str_, np.str_, np.int_, np.float_]
fld_cols = dict(zip(fld_names, fld_type))

df = pd.read_csv(filename, sep=',', usecols=fld_names, dtype=fld_cols, skipinitialspace=True)

# %%
df.info()

# %%
df.head()

# %%

# %%
# Read the NHM-PRMS_data_release.csv file
fhdl = open(rel_filename, 'r', encoding='ascii')

# Get the column headers
hdr_flds = fhdl.readline().strip('\n').split(' ')
fhdl.close()

hdr_flds.pop(0)    # Remove the 'date' header entry
gage_id_map = {ii + 3: vv.split('_')[4] for ii, vv in enumerate(hdr_flds)}

# Read in the data
sim_df = pd.read_csv(rel_filename, sep=' ', header=None, skiprows=2, skipinitialspace=True,
                     parse_dates={'date': [0, 1, 2]},
                     engine='c', memory_map=True, index_col='date')

# Rename columns to gage ids
sim_df.rename(columns=gage_id_map, inplace=True)

# %%
sim_df.info()

# %%
# Pull observed streamflow from the source netcdf files used by Bandit
# /Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/poi_data

# %%

# %%

# %%
