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
import pandas as pd
from pyPRMS.ParamDb import ParamDb
# from pyPRMS.ParameterFile import ParameterFile
# from pyPRMS.ParameterNetCDF import ParameterNetCDF

# %%
src_param = '/Users/pnorton/Downloads/soil_rechr_max_frac_hru_id_nhm_061621.csv'

workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS'

outdir = '/Users/pnorton/Projects/National_Hydrology_Model/FY2021_projects/model_support/NHMv1.1'
newpdb_dir = f'{outdir}/20210617_FIX_soil_rechr_max_frac/paramdb_new'

hru_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_nhrusim_Topology.gdb'
hru_layer_name = 'nhruv11_sim30'
hru_shape_key='nhru_v11'

# Segment lines
seg_geodatabase = '/Users/pnorton/Projects/National_Hydrology_Model/Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
seg_layer_name = 'nsegment_v11'
seg_shape_key = 'nsegment_v11'

# %%
pdb = ParamDb(workdir, verbose=True, verify=True)

# %%
pdb.parameters.check()

# %%
new_vals = pd.read_csv(src_param, sep=',', usecols=[0, 1], index_col=0)

new_vals.head()

# %%
cparam = 'soil_rechr_max_frac'

for row in new_vals.itertuples(index=True, name='Pandas'):
#     print(row.Index, row.soil_rechr_max_frac)
    
    pdb.parameters.update_element(cparam, row.Index, row.soil_rechr_max_frac)

# %%

# %%
pdb.parameters.shapefile_hrus(hru_geodatabase, layer_name=hru_layer_name, shape_key=hru_shape_key)
# pdb.parameters.shapefile_segments(seg_geodatabase, layer_name=seg_layer_name, shape_key=seg_shape_key)

# %%
cparam = 'soil_rechr_max_frac'

# pdb.parameters.plot(cparam, limits='valid', linewidth=0.0, edgecolor='whitesmoke', 
#                     mask_defaults='darkgrey', cmap='gist_earth_r', output_dir='/Users/pnorton/tmp')
pdb.parameters.plot(cparam, limits='valid', linewidth=0.0, edgecolor='whitesmoke', 
                    mask_defaults='darkgrey', cmap='gist_earth_r')

# %%
pdb.write_paramdb(newpdb_dir)

# %%

# %%
