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
from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.plot_helpers import get_projection

# %%
# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb'
# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb_CONUS'
# work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/paramdb_daymet_CONUS'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20200415_red_river_v2'
filename = f'{workdir}/myparam.param'

# shpfile = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
# layer_name = 'nsegment_v11'

# HRU polygons
# hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
# hru_layer_name = 'nhruv11_sim30'
# hru_shape_key = 'nhru_v11'

hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp'
hru_layer_name = None
hru_shape_key='hru_id_nat'

# Segment lines
seg_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
seg_layer_name = 'nsegment_v11'
seg_shape_key = 'nsegment_v11'

# %%
pdb = ParameterFile(filename, verbose=True, verify=True)
# pdb = ParamDb(paramdb_dir=work_dir, verbose=True, verify=True)

# %%
# Load the HRU shapes
pdb.parameters.shapefile_hrus(f'{workdir}/GIS/HRU_subset.shp', layer_name=None, shape_key='nhru_v11')
# pdb.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)

seg_shape_key = 'nsegment_v'
pdb.parameters.shapefile_segments(f'{workdir}/GIS/Segments_subset.shp', layer_name=None, shape_key=seg_shape_key)

# pdb.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %%
# pdb.parameters.plot('hru_area', output_dir=None, linewidth=6.0, 
#                     facecolor='snow', edgecolor='whitesmoke', 
#                     vary_color=True, vary_width=True, cmap='viridis')
pdb.parameters.plot('slowcoef_sq', output_dir=None, cmap='viridis')

# %%
pdb.parameters.plot('seg_width', output_dir=None, linewidth=6.0, 
                    facecolor='snow', edgecolor='whitesmoke', 
                    vary_color=True, vary_width=True, cmap='tab20')

# %%
pdb.parameters['seg_width'].stats()

# %%
