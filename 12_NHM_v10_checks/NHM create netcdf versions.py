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
from collections import OrderedDict
from pyPRMS.ParamDb import ParamDb
# from pyPRMS.ParamDbRegion import ParamDbRegion

# %%
basedir = '/Users/pnorton/Projects/National_Hydrology_Model'

# base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv10/byHRU_test1'
# work_dir = f'{base_dir}/testDb'

# V1.1 gridmet paramdb; master branch
# /Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS
workdir = f'{basedir}/datasets/paramdb_v11'
paramdb = f'{workdir}/paramdb_v11_gridmet_CONUS'
outdir = f'{basedir}/calibrations/NHMv11/gridmet_byHRU'
ncname = 'A_master_20210520.nc'

# V1.1 gridmet byHRU calibration results
# /Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/gridmet_byHRU
# workdir = f'{basedir}/calibrations/NHMv11/gridmet_byHRU'
# paramdb = f'{workdir}/byHRU_pdb'
# outdir = f'{workdir}'
# ncname = 'B_byHRU.nc'


# V1.0 daymet 
# workdir = f'{basedir}/datasets/paramdb_v10'
# paramdb = f'{workdir}/paramdb_v10_daymet_CONUS'
# outdir = f'{workdir}/parameters'
# ncname = 'D_new_byHW_musk_obs.nc'

# v1.0 maurer
# workdir = f'{basedir}/datasets/paramdb_v10'
# paramdb = f'{workdir}/maurer_fix_work/maurer_byHRUdyn_musk_obs'
# outdir = f'{workdir}/maurer_fix_work'
# ncname = 'G_byHRUdyn_musk_obs_new.nc'


# byHRU merge test results
# workdir = f'{basedir}/calibrations/NHMv10/byHRU_test1'
# paramdb = f'{workdir}/testDb'

# Lauren's Daymet repo
# workdir = f'{basedir}/datasets/paramdb_v10'
# paramdb = f'{workdir}/lauren_repo/Daymet'
# paramdb = '/Users/pnorton/tmp/LH_fb8729'


# workdir = f'{basedir}/calibrations/NHMv10/byHRU_test1/'
# paramdb = f'{workdir}/paramdb_v10_daymet_CONUS'
# outdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/parameters'


# workdir = f'{basedir}/datasets/paramdb_v10'
# paramdb = f'{workdir}/nhmparamdb_lhay_daymet_byHRU'
# paramdb = f'{workdir}/paramdb_v10_daymet_CONUS'
# outdir = f'{workdir}/parameters'

# pdb = ParamDbRegion(paramdb_dir=paramdb, verbose=True, verify=True)

pdb = ParamDb(paramdb_dir=paramdb, verbose=True, verify=True)

# %%
pdb.parameters.check()

# %%
# pdb.write_parameter_file(filename=f'{outdir}/master.param')

# %%
pdb.write_netcdf(filename=f'{outdir}/{ncname}')

# %%
# Convert to new-style paramdb
pdb.write_paramdb('/Users/pnorton/tmp/LH_fb8729')

# %%
import pandas as pd

param_stats = []

for pp in pdb.parameters.values():
    param_stats.append(pp.stats())
    
df = pd.DataFrame.from_records(param_stats, columns=['name', 'min', 'max', 'mean', 'median'])

# %%
df

# %%
aa = pdb.parameters.get_dataframe('snarea_curve')

# %%
aa

# %%
bb = aa.loc[7285].tolist()

for xx in bb:
    print(xx)

# %%
aa = pdb.parameters.get_dataframe('jh_coef')
aa.head()

# %%
bb = aa.loc[6918].tolist()
for xx in bb:
    print(xx)

# %%
aa = pdb.parameters.get_dataframe('carea_max')
aa.loc[7285].tolist()



# %%
pdb.parameters.get('cecn_coef').stats()

# %%
aa = pdb.parameters.get_dataframe('mann_n')
# aa.head(10)

for xx in [3705, 4226, 3703]:
    print(aa.loc[xx].tolist())

# %%
