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
workdir = '/Volumes/USGS_NHM1/calibrations/NHMv10/MAURER_releases/fix_runs/byHRU_musk'
run_out_dir = f'{workdir}/output'


obs_file = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/NHM_v1.0/poi_data/*.nc'
sim_file = f'{run_out_dir}/NHM-PRMS_data_release.csv'

ctl_file = f'{workdir}/NHM-PRMS_data_release.control'

# %%
import os
import xarray as xr

from collections import OrderedDict
from pyPRMS.ControlFile import ControlFile
from pyPRMS.ParameterFile import ParameterFile


# %%
def get_var(var, df, gageids, stdate=None, endate=None):
    if stdate is None or endate is None:
        data = df[var].loc[gageids].to_pandas()
    else:
        data = df[var].loc[gageids, stdate:endate].to_pandas()
    return data

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

# %%
poi_flds, poi_seg_flds = read_csv_header(sim_file)

# %%
poi_flds

# %%
ctl = ControlFile(f'{ctl_file}')

param_file = os.path.normpath(os.path.join(f'{workdir}/{ctl.get("param_file").values}'))
pfile = ParameterFile(param_file, verbose=True, verify=True)

# Get the list of POI IDs in the parameter file
poi_gage_id = pfile.parameters['poi_gage_id']
poi_gage_id_list = poi_gage_id.data.tolist()

# %%
len(poi_gage_id_list)

# %%

# %%
obs_df = xr.open_mfdataset(obs_file, chunks={'poi_id': 2000}, combine='nested',
                           concat_dim='poi_id', decode_cf=True, engine='netcdf4')

# %%
obs_pois = obs_df['poi_id'].to_pandas().tolist()

# %%
set(poi_gage_id_list) - set(obs_pois)

# %%
poiname_list = get_var('poi_name', obs_df, poi_gage_id_list).tolist()

# %%
len(poiname_list)

# %%

# %%

# %%

# %%
