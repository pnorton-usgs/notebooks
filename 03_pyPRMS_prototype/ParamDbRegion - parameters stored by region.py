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
#     display_name: Python [conda env:gis_38]
#     language: python
#     name: conda-env-gis_38-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
from pyPRMS.ParamDbRegion import ParamDbRegion
from pyPRMS.ParameterFile import ParameterFile

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/crap/Maurer'

pdb = ParamDbRegion(paramdb_dir=workdir, verbose=True, verify=True)

# %%
print(pdb.parameters['hru_deplcrv'])

# %%
snow_index = pdb.parameters['hru_deplcrv'].data[0]
pdb.parameters['snarea_curve'].data.reshape((-1, 11))[snow_index-1, :]

# %% [markdown]
# ## Check snarea_curve in py27 parameter file

# %%
workdir2 = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190418_hw_6000/py27'

# %%
pfile = ParameterFile('{}/myparam.param'.format(workdir2))

# %%
snow_index = pfile.parameters['hru_deplcrv'].data[40]
pfile.parameters['snarea_curve'].data.reshape((-1, 11))[snow_index-1, :]

# %% [markdown]
# ## Check snarea_curve in py37 parameter file

# %%
workdir3 = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190418_hw_6000'
pfile2 = ParameterFile('{}/myparam.param'.format(workdir3))

# %%
snow_index = pfile2.parameters['hru_deplcrv'].data[40]
pfile2.parameters['snarea_curve'].data.reshape((-1, 11))[snow_index-1, :]

# %%
hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp'
hru_layer_name = None
hru_shape_key='hru_id_nat'

# %%
pdb.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)

# %%
pdb.parameters.plot('dday_intcp', output_dir=None, cmap='tab20')

# %%
