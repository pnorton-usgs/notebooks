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
from future.utils import iteritems

from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ControlFile import ControlFile
from pyPRMS.ValidParams import ValidParams
# reload(pf)

# %% [markdown]
# ### Load a parameter file

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/tmp/nvwsc'
filename = '{}/Desch.params6'.format(workdir)

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/Study_areas/sagehen/input/prms'
# filename = '{}/prms.params'.format(workdir)

# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/prms5/pipestem'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/src/prms5/red_river_of_the_south'
# filename = '{}/myparam.param'.format(workdir)

# filename = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/regions/subset_testing/crap.param'

pfile = ParameterFile(filename)

# Set pointers to the parameters and dimensions within the parameterSet
# This simplifies code for access to the members of parameters and dimensions
myparams = pfile.parameters
mydims = pfile.dimensions

# %%
# Load master list of valid parameters
vpdb = ValidParams('/Users/pnorton/PycharmProjects/pyPRMS/pyPRMS/xml/parameters.xml')

# %%
# Read the control file and get a list of requested modules
ctl = ControlFile('{}/control'.format(workdir))

modules_used = ctl.modules.values()

# %%
ctl.modules.values()

# %%
# Get the list of parameters that are required for the selected modules
required_params = vpdb.get_params_for_modules(modules=modules_used)
print('Number of required parameters: {}'.format(len(required_params)))

# %%
type(required_params)

# %%
len(pfile.parameters.keys())

# %%
pfile.reduce_by_modules(ctl)

# %%
len(pfile.parameters.keys())

# %%
header=['Written by parameterSet.write_parameter_file()','####']
pfile.write_parameter_file('crap.param', header=header)

# %%
pfile.write_netcdf('rrs_reduced.nc')

# %%
print(sorted(list(required_params)))

# %%
required_params.difference(set(pfile.parameters.keys()))

# %%
remove_list = set(pfile.parameters.keys()).difference(required_params)
print(remove_list)

# %%
set(pfile.parameters.keys()).difference(required_params)

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
