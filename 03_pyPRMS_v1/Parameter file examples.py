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
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import pandas as pd
import xml.dom.minidom as minidom
import xml.etree.ElementTree as xmlET

from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ControlFile import ControlFile
from pyPRMS.ValidParams import ValidParams

# %%
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20220214_gm_pipestem'

output_dir = f'{work_dir}/output'

filename = f'{work_dir}/myparam.param'

gis_dir = f'{work_dir}/GIS'

pfile = ParameterFile(filename, verbose=True, verify=True)

# Set pointers to the parameters and dimensions within the parameterSet
# This simplifies code for access to the members of parameters and dimensions
myparams = pfile.parameters
mydims = pfile.dimensions

# %% [markdown]
# ## Dimensions

# %%
# Print the information for all the dimensions
print('='*10 + ' Dimensions string representation ' + '='*10)
print(mydims)

# %% [markdown]
# ### Dimension names

# %%
print('='*10 + ' Dimension names ' + '='*10)

print(mydims.keys())

# %% [markdown]
# ### Dimension objects

# %%
print('='*10 + ' Dimension objects ' + '='*10)
print(list(mydims.values()))

# %% [markdown]
# ### Iterate over dimensions

# %%
for kk in mydims.values():
    print(f'{kk.name} = {kk.size}')

# %% [markdown]
# ## Parameters

# %% [markdown]
# ### Get information about a parameter

# %%
curr_param = 'gwflow_coef'

print(myparams.get(curr_param))

# %% [markdown]
# ### Accessing different parts of a parameter

# %% [markdown]
# #### Get parameter data as an ndarray

# %%
# Returns an ndarray of the data for a parameter
print('='*10 + f' Data for {curr_param} '+ '='*10)
# print(myparams.get(curr_param).data)
print(myparams[curr_param].data)

# %% [markdown]
# #### Get parameter data as a dataframe

# %%
myparams.get_dataframe(curr_param)

# %% [markdown]
# #### Get the defined dimensions for the parameter

# %%
print('='*10 + f' Dimensions defined for {curr_param} ' + '='*10)
print(myparams[curr_param].dimensions)

# %% [markdown]
# #### Datatype for parameter

# %%
# Datatype for the parameter
print('='*10 + f' Datatype for {curr_param} ' + '='*10)
print(myparams[curr_param].datatype)

# %% [markdown]
# ### 2D Parameters returned as a dataframe
# The column names for a 2D parameter are uniquely labeled with the parameter name and index values from the second dimension

# %%
myparams.get_dataframe('tmax_allsnow').head()

# %% [markdown]
# ### Checking basic consistency of parameters
# <li>Can check a single parameter or all parameters in a ParameterSet</li>
# <li>Checks to make sure the total number of values equals the total declared dimension size</li>
# <li>Outputs either OK or BAD. This could be modified to return error codes instead</li>
# <li>Additional consistency checks could be added</li>

# %% [markdown]
# #### Check parameter for internal consistency

# %%
myparams.get(curr_param).check()

# %% [markdown]
# #### Check consistency of all parameters in the Parameter file

# %%
# Check consistency of all parameters in the ParameterSet
myparams.check()

# %% [markdown]
# #### Get some simple statistics

# %%
myparams[curr_param].stats()

# %%

# %%

# %% [markdown]
# ### Remove parameter

# %%
tst_param = 'albset_snm'

# The current parameter set has the tst_param
print('='*10 + f' Before removing {tst_param} ' + '='*10)
print(f'{tst_param} exists = {myparams.exists(tst_param)}')

# We can remove that parameter
# pfile.parameterset.parameters.remove('albset_snm')
myparams.remove(tst_param)

# And now it's gone
print('='*10 + f' After removing {tst_param} '.format(tst_param) + '='*10)
print(f'{tst_param} exists = {myparams.exists(tst_param)}')

# %% [markdown]
# ### Iterate over parameters

# %%
# Loop through the parameter objects in the ParameterSet and display name, datatype, and dimension info
for vv in list(myparams.values()):
    print(f'{vv.name}: datatype = {vv.datatype}')
    
    for dd in vv.dimensions.values():
        print(f'\t{dd.name} = {dd.size}')

# %% [markdown]
# ### Write new parameter file

# %%
pfile.write_parameter_file(filename=f'{output_dir}/new_param_file', header=['Written by ParameterFile','Updated from calibration'])

# %% [markdown]
# ### Get dimensions in xml format

# %%
xmlstr = minidom.parseString(xmlET.tostring(pfile.xml_global_dimensions)).toprettyxml(indent='    ')
print(xmlstr)

# %% [markdown]
# ### Get parameters in xml format

# %%
xmlstr = minidom.parseString(xmlET.tostring(pfile.xml_global_parameters)).toprettyxml(indent='    ')
print(xmlstr)

# %% [markdown]
# ### Create table with parameter information

# %%
out_list = []
for pp in pfile.parameters.values():
    dims = []
    for dd in pp.dimensions.values():
        dims.append(dd.name)

    out_list.append([pp.name, pp.datatype, pp.units, pp.description, pp.minimum, pp.maximum, pp.default, dims])

col_names = ['parameter', 'datatype', 'units', 'description', 'minimum', 'maximum',
             'default', 'dimensions']
df = pd.DataFrame.from_records(out_list, columns=col_names)
df.to_csv(f'{output_dir}/parameters.csv', sep='\t', index=False)

# %% [markdown]
# ### Plotting parameters for a model extracted by Bandit

# %% [markdown]
# #### Load GIS information first
# The layer_name is always None for extracted model shapefiles
# The HRU shape_key for an extracted domain is different between v1.0 and v1.1 extractions. 
# <li>v1.0 shape_key = 'hru_id_nat'</li>
# <li>v1.1 shape_key = 'nhru_v1_1'</li>
# The segment shape_key is also different between v1.0 and v1.1 extractions.
# <li>v1.0 shape_key = 'seg_id_nat'</li>
# <li>v1.1 shape_key = 'nsegment_v'</li>

# %%
pfile.parameters.shapefile_hrus(f'{gis_dir}/HRU_subset.shp', layer_name=None, shape_key='nhru_v1_1')

# %% [markdown]
# #### Plot a HRU parameter

# %%
# Save plot to a file
# pfile.parameters.plot('snowinfil_max', output_dir=f'{workdir}', cmap='terrain')

# Display plot in notebook
pfile.parameters.plot('snowinfil_max', cmap='terrain')

# %% [markdown]
# ### Plot a segment parameter

# %%
pfile.parameters.shapefile_segments(f'{gis_dir}/Segments_subset.shp', layer_name=None, shape_key='nsegment_v')

# %%
# NOTE: 2022-02-14 PAN - This doesn't work for some reason. I'll look at pyPRMS and see what's going on.
pfile.parameters.plot('seg_width', linewidth=6.0, facecolor='snow', edgecolor='whitesmoke', 
                      vary_color=True, vary_width=True, cmap='tab20')

# %%

# %%

# %%

# %%

# %%

# %%

# %%
