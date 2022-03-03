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
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ControlFile import ControlFile
from pyPRMS.ValidParams import ValidParams
# reload(pf)

# %% [markdown]
# ### Load a parameter file

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/tmp/nvwsc'
# filename = '{}/Desch.params6'.format(workdir)

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/Study_areas/sagehen/input/prms'
# filename = '{}/prms.params'.format(workdir)

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/prms5/pipestem'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/prms5/red_river_of_the_south'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20191017_newnhm'
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/pipestem'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20200512_pipestem'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/extraction_requests/20190626_columbia_plateau'

filename = '{}/myparam.param'.format(workdir)

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/sample_hrus'
# filename = f'{workdir}/ALL.SCE_HRU1'


gis_dir = f'{workdir}/GIS'

# filename = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/regions/subset_testing/crap.param'

pfile = ParameterFile(filename, verbose=True, verify=True)

# Set pointers to the parameters and dimensions within the parameterSet
# This simplifies code for access to the members of parameters and dimensions
myparams = pfile.parameters
mydims = pfile.dimensions

# %%
pfile.dimensions['nsegment'].size

# %% [markdown]
# ### List the dimensions that are included in the parameterSet

# %%
# Print the information for all the dimensions
# Two approaches are shown
# print(pfile.parameterset.dimensions)
print('='*10 + 'Dimensions string representation' + '='*10)
print(mydims)

# List all the dimension names (two approaches)
print('='*10 + 'Dimension names' + '='*10)
# print(pfile.parameterset.dimensions.keys())
print(list(mydims.keys()))

# List all the dimension objects (two approaches)
print('='*10 + 'Dimension objects' + '='*10)
# print(pfile.parameterset.dimensions.values())
print(list(mydims.values()))

# %% [markdown]
# # Iterate through the available dimensions

# %%
for kk in mydims.values():
    print('{} = {}'.format(kk.name, kk.size))

# %% [markdown]
# ### Get information about a parameter

# %%
curr_param = 'gwflow_coef'

# Here's one way to get parameter information (three approaches)
# print(pfile.parameterset.parameters.get('gwflow_coef'))
print(myparams.get(curr_param))
# print(myparams['gwflow_coef'])

# %% [markdown]
# ### Accessing different parts of the parameter
# #### These are things like the associated dimensions, datatype, and data

# %%
# Returns an ndarray of the data for a parameter
print('='*10 + 'Data for {}'.format(curr_param) + '='*10)
# print(myparams.get(curr_param).data)
print(myparams[curr_param].data)

# Dimensions defined for the parameter
print('='*10 + 'Dimensions defined for {}'.format(curr_param) + '='*10)
print(myparams[curr_param].dimensions)

# Datatype for the parameter
print('='*10 + 'Datatype for {}'.format(curr_param) + '='*10)
print(myparams[curr_param].datatype)

# %% [markdown]
# ### Get the parameter data returned as a Pandas DataFrame
# The dataframe will contain the selected data for the parameter. If the parameter is associated with either nhm_id or nhm_seg then that information is included in the dataframe as the index.

# %%
# Get parameter data as a dataframe
# Using .head() to only display the first five rows

# There are multiple approaches to doing this
# pfile.parameterset.parameters.get_DataFrame(curr_param).head()
myparams.get_dataframe(curr_param).head()

# Using the .as_dataframe property of an individual parameter
# Will return just the parameter data with no associated 
# nhm_id or nhm_seg parameters
# myparams[curr_param].as_dataframe.head()

# %%
myparams.get_dataframe('nhm_id')
# myparams.get_dataframe(curr_param).loc[57874].values[0]

# %% [markdown]
# ### 2D Parameters returned as a dataframe
# The column names for a 2D parameter are uniquely labeled with the parameter name and index values from the second dimension

# %%
# .head() is used to limit output to the top five rows

# pfile.parameterset.parameters.get_DataFrame('tmax_allsnow')
myparams.get_dataframe('tmax_allsnow').head()

# %% [markdown]
# ### Get the names of the dimensions for a parameter

# %%
# Get the names of the dimensions for the selected parameter
# pfile.parameterset.parameters[curr_param].dimensions.keys()
list(myparams.get(curr_param).dimensions.keys())

# %% [markdown]
# ### Get a dictionary of parameter objects from the ParameterSet

# %%
# Get the dictionary of parameter objects
# pfile.parameterset.parameters.items()
list(myparams.items())

# %% [markdown]
# ### Access parameters by index

# %%
# Get the dimensions for second parameter (e.g. tstorm_mo)
# print(pfile.parameterset.parameters.items()[1][1].dimensions)

print(list(myparams.items())[1][1].dimensions)
#                               ^---- index 1 of the selected parameter contains the parameter object
#                            ^---- second parameter from parameters.items()
# All indices are zero-based

# %%
# myparams[1][1].dimensions

# %% [markdown]
# ### Checking basic consistency of parameters
# <li>Can check a single parameter or all parameters in a ParameterSet</li>
# <li>Checks to make sure the total number of values equals the total declared dimension size</li>
# <li>Outputs either OK or BAD. This could be modified to return error codes instead</li>
# <li>Additional consistency checks could be added</li>

# %%
# Check a nparameter for internal consistency
# pfile.parameterset.parameters[curr_param].check()
myparams.get(curr_param).check()

# %%
# Check consistency of all parameters in the ParameterSet
# pfile.parameterset.parameters.check()
myparams.check()

# %% [markdown]
# ### You can remove a parameter from the ParameterSet

# %%
tst_param = 'albset_snm'

# The current parameter set has the tst_param
print('='*10 + 'Before removing {}'.format(tst_param) + '='*10)
print('{} exists = {}'.format(tst_param, myparams.exists(tst_param)))

# We can remove that parameter
# pfile.parameterset.parameters.remove('albset_snm')
myparams.remove(tst_param)

# And now it's gone
print('='*10 + 'After removing {}'.format(tst_param) + '='*10)
print('{} exists = {}'.format(tst_param, myparams.exists(tst_param)))

# %% [markdown]
# ### Example of iterating through the ParameterSet

# %%
# Loop through the parameter objects in the ParameterSet and display name, datatype, and dimension info
for vv in list(myparams.values()):
    print('{}: datatype = {}'.format(vv.name, vv.datatype))
    
    for dd in vv.dimensions.values():
        print('\t{} = {}'.format(dd.name, dd.size))

# %% [markdown]
# ### Writing a new parameter file
# Output parameter file will have the same order of dimensions and parameters as the input parameter file

# %%
pfile.write_parameter_file(filename='new_param_file', header=['Written by ParameterFile','Updated from calibration'])

# %%
from future.utils import iteritems

import pyPRMS.ParameterFile as pf
reload(pf)

# %% [markdown]
# ### Load a parameter file

# %%
filename = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/regions/subset_testing/crap.param'

pfile = pf.ParameterFile(filename)

# Set pointers to the parameters and dimensions within the parameterSet
# This simplifies code for access to the members of parameters and dimensions
myparams = pfile.parameterset.parameters
mydims = pfile.parameterset.dimensions

# %%
### Get the position of a dimension for a variable

print(myparams.get('tstorm_mo').dimensions.get_position('nmonths'))
# print(myparams.get('tstorm_mo'))

# %% [markdown]
# ### Write parameters in paramDb format

# %%
pfile.write_paramdb('{}/paramdb'.format(workdir))

# %% [markdown]
# ### Get dimensions in xml format

# %%
import xml.dom.minidom as minidom
import xml.etree.ElementTree as xmlET

xmlstr = minidom.parseString(xmlET.tostring(pfile.xml_global_dimensions)).toprettyxml(indent='    ')
print(xmlstr)
# with open('crap.xml', 'w') as f:
#     f.write(xmlstr.encode('utf-8'))

# %% [markdown]
# ### Get parameters in xml format

# %%
xmlstr = minidom.parseString(xmlET.tostring(pfile.xml_global_parameters)).toprettyxml(indent='    ')
print(xmlstr)

# %%

# %%
# Loop through the parameter objects in the ParameterSet and display name, datatype, and dimension info
for vv in myparams.values():
#     print('{}'.format(vv.name))
    print('{}: datatype = {}'.format(vv.name, vv.datatype))
    
#     for dd in vv.dimensions.values():
#         print('\t{} = {}'.format(dd.name, dd.size))

# %%
# Load master list of valid parameters
vpdb = ValidParams()


# %%
# Read the control file and get a list of requested modules
ctl = ControlFile('{}/control.default'.format(workdir))

modules_used = ctl.modules.values()

# for xx in ctl.modules.keys():
#     if xx == 'precip_module':
#         if ctl.modules[xx] == 'climate_hru':
#             modules_used.append('precipitation_hru')
#     elif xx == 'temp_module':
#         if ctl.modules[xx] == 'climate_hru':
#             modules_used.append('temperature_hru')
#     else:
#         modules_used.append(ctl.modules[xx])
# ctl.write('control.rrs')

# %%
ctl.modules.values()

# %%
# Get the list of parameters that are required for the selected modules
required_params = vpdb.get_params_for_modules(modules=modules_used)
print('Number of required parameters: {}'.format(len(required_params)))

# %%
print(vpdb.parameters['hru_area'])

# %%
len(pfile.parameters.keys())

# %%
pfile.remove_unneeded_parameters(required_params=required_params)

# %%
len(pfile.parameters.keys())

# %%
header=['Written by parameterSet.write_parameter_file()','####']
pfile.write_parameter_file('crap.param', header=header)

# %%
pfile.write_netcdf('rrs_reduced.nc')

# %%
b = myparams['gwflow_coef'].data

# %%
b.size

# %%
b[0] = 95.

# %%
b

# %%
myparams['gwflow_coef'].data = b

# %%
myparams['snarea_curve'].unique().size

# %%
myparams.get_dataframe('snarea_curve').head()

# %%
myparams['hru_deplcrv'].unique().size

# %%
myparams['snarea_curve'].data.reshape((-1, 11)).shape[0]

# %%
myparams['tstorm_mo'].dimensions.get_position('nmonths')

# %%
# Load master list of valid parameters
vpdb = ValidParams()
vpdb.parameters['tstorm_mo'].dimensions.keys()

# %%
print(vpdb.parameters['snarea_curve'])

# %%
aa = vpdb.parameters['tstorm_mo'].dimensions.copy()

# %%
from pyPRMS.Dimensions import Dimensions
isinstance(aa, Dimensions)

# %%
type(aa)

# %%
type(vpdb.parameters['tstorm_mo'].dimensions)

# %%
for kk, vv in iteritems(aa):
    print(kk, vv.size)

# %%
pfile.expand_parameter('tstorm_mo')

# %%
new_sizes = []

for kk, vv in iteritems(aa):
    new_sizes.append(vv.size)
    
print(new_sizes)

# %%
new_sizes = [vv.size for vv in aa.values()]
print(new_sizes)

# %%
print(pfile.parameters['lake_evap_adj'])

# %%
pfile.parameters['lake_evap_adj'].data

# %%
pfile.expand_parameter('lake_evap_adj')

# %%
print(pfile.parameters['lake_evap_adj'])

# %%
pfile.parameters['lake_evap_adj'].data

# %%
pfile.expand_parameter('tstorm_mo')

# %%
pfile.degenerate_params()

# %%
print(vpdb.parameters['jh_coef'])

# %%
pfile.parameters['jh_coef'].data

# %%
pfile.expand_parameter('jh_coef')

# %%
pfile.parameters['jh_coef'].data[127,:]

# %%
pfile.parameters['jh_coef'].data.shape

# %%
pfile.parameters['tmax_adj'].data

# %%
pfile.expand_parameter('tmax_adj')

# %%
pfile.parameters['tmax_adj'].data[:,11]

# %%
pfile.parameters['tmax_adj'].data.shape

# %%
print(pfile.parameters['snarea_curve'])

# %%
pfile.parameters['snarea_curve'].data

# %%
print(pfile.parameters['hru_deplcrv'])

# %%
pfile.expand_parameter('hru_deplcrv')

# %%
print(pfile.parameters['hru_deplcrv'])

# %%
pfile.parameters['hru_deplcrv'].data

# %%
print(pfile.parameters['snarea_curve'])

# %%
pfile.parameters['snarea_curve'].data.reshape((-1, 11))[8,:]

# %%
pfile.parameters['snarea_curve'].data

# %%
pfile.expand_parameter('hru_deplcrv')

# %%
pfile.write_parameter_file('{}/Desch.params6.test'.format(workdir))

# %%
print(pfile.parameters['tstorm_mo'])

# %%
pfile.updated_params

# %%
print(myparams['tmax_allrain_offset'])

# %%
myparams['poi_gage_id'].data

# %%
for xx in myparams['tmax_allrain_offset'].dimensions.items():
    print(xx)

# %%
# pfile.parameters.shapefile_hrus('/Users/pnorton/Projects/National_Hydrology_Model/notebooks/GIS/all_nhru_simple_lcc/nhruNationalIdentifier.shp', 
#                               layer_name=None, shape_key='hru_id_nat')
pfile.parameters.shapefile_hrus(f'{gis_dir}/HRU_subset.shp', 
                              layer_name=None, shape_key='hru_id_nat')

# %%
pfile.parameters.plot('snowinfil_max', output_dir=f'{workdir}', cmap='terrain')

# %%
myparams['snowinfil_max'].stats()

# %%
