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
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
from pyPRMS.ParameterFile import ParameterFile

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20200415_red_river_v2'
filename = f'{workdir}/myparam.param'

outdir = workdir
# outdir = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/v1_plots'
# outdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/SNODAS/unmasked'
# outdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/ERA5'

# %%
pdb = ParameterFile(filename, verbose=True, verify=True)

# %%
pdb.parameters.shapefile_hrus(f'{workdir}/GIS/HRU_subset.shp', layer_name=None, shape_key='nhru_v11')

# pdb.parameters.shapefile_hrus('/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp', 
#                               layer_name=None, shape_key='hru_id_nat')

# pdb.parameters.shapefile_hrus('/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb', 
#                               layer_name='nhruv11_sim30', shape_key='nhru_v11')

# %%
pdb.parameters.plot('hru_deplcrv', output_dir=None)

# %%
crv = pdb.parameters.get_dataframe('snarea_curve')
crv.T.plot(cmap = plt.cm.get_cmap('rainbow', crv.shape[0]))

# %%
pdb.parameters['hru_deplcrv'].stats()

# %%
