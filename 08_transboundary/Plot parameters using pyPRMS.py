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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
from pyPRMS.ParamDb import ParamDb

# %%
workdir = '/Users/pnorton/tmp/check_paramdb_v11'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_CONUS'
# outdir = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/v11_SRTM_plots_b'

# outdir = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/v1_plots'
outdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/SNODAS/unmasked'
# outdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/ERA5'


# Segment lines
seg_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
seg_layer_name = 'nsegment_v11'
seg_shape_key = 'nsegment_v11'

# %%
pdb = ParamDb(paramdb_dir=workdir, verbose=True, verify=True)

# %%
# pdb.parameters.shapefile_hrus('/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp', 
#                               layer_name=None, shape_key='hru_id_nat')

pdb.parameters.shapefile_hrus('/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb', 
                              layer_name='nhruv11_sim30', shape_key='nhru_v11')

pdb.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %%
pdb.parameters.plot('hru_deplcrv', output_dir=None, cmap='tab10')

# %%
# Plot the snow depletion curves
crv = pdb.parameters.get_dataframe('snarea_curve')
crv.T.plot(cmap=plt.cm.get_cmap('tab10', crv.shape[0]))

# %%
pdb.parameters['hru_deplcrv'].stats()

# %%
pdb.parameters.plot('mann_n', output_dir='/Users/pnorton/tmp/check_paramdb_v11', linewidth=6.0, 
                    facecolor='snow', edgecolor='whitesmoke', 
                    vary_color=True, vary_width=True, cmap='tab20')

# %%
