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
#     display_name: Python [conda env:bandit_py3]
#     language: python
#     name: conda-env-bandit_py3-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import datetime

from Bandit.model_output import ModelOutput

# %%
workdir = '/Volumes/parker_rocks/NHM_output/netcdf'
filename = f'{workdir}/tmaxf.nc'

st_date = datetime.datetime(1980, 10, 1)
en_date = datetime.datetime(1980, 10, 31)

the_hrus = [1000, 999, 1001, 1002]


# %%
modout = ModelOutput(filename=filename, varname='tmaxf', startdate=st_date, enddate=en_date, nhm_hrus=the_hrus)

# %%
modout.dataset

# %%
modout.write_csv('/Users/pnorton/Projects/National_Hydrology_Model/tmp')

# %%
modout.write_netcdf('/Users/pnorton/Projects/National_Hydrology_Model/tmp')

# %%
