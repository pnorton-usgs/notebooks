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
#     display_name: Python [conda env:idp_bandit]
#     language: python
#     name: conda-env-idp_bandit-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import xml.dom.minidom as minidom
import xml.etree.ElementTree as xmlET

import pyPRMS.NhmParamDb as nhm
import pandas as pd
# reload(nhm)

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/paramdb_v2'

pdb = nhm.NhmParamDb(workdir)

# %%
pdb.available_parameters

# %%
# Get specified HRUs for given parameter
# for pp in pdb.available_parameters:
#     thedata = pdb.parameters.get(pp)
# thedata = pdb.parameters.get('gwflow_coef').data
# thedata[52094]
print(pdb.parameters.get('tmax_allsnow').dimensions.keys())
for pp in pdb.available_parameters:
    if set(pdb.parameters.get(pp).dimensions.keys()).intersection({'nhru', 'ngw', 'nssr'}):
        pd.DataFrame(pdb.parameters.get_dataframe(pp).loc[(11485, 52095),]).to_csv('/Users/pnorton/tmp/crapper/{}.csv'.format(pp))
#         print(pd.DataFrame(pdb.parameters.get_dataframe(pp).loc[(11485, 52095),]))

# %%
# print(pdb.parameters['tstorm_mo'].dimensions.get_position('nmonths'))
print(pdb.parameters.items())

# %%
print(pdb.dimensions)

# %%
print(pdb.parameters.get('dday_slope'))

# %%
print(pdb.parameters.get('gwflow_coef').data)

# %%
print(pdb.parameters['snarea_curve'])

# %%
print(pdb.parameters['hru_deplcrv'].data)

# %%
print(pdb.parameters.get('tmax_allsnow').modules)

# %%
print(pdb.parameters['tstorm_mo'].dimensions.tostructure())

# %%
bb = []
bb.insert(1, 'test')
bb.insert(0, 'firstpos')

print(bb)

# %%
print(pdb.parameters['tstorm_mo'].tostructure())

# %%
print(pdb.parameters['snarea_curve'].tostructure())

# %%
pdb.parameters['tstorm_mo'].get_dimsize_by_index(1)

# %%
print(pdb.parameters['tstorm_mo'].tolist())

# %%
pdb.parameters.get('dday_slope').as_dataframe

# %%
aa = pdb.parameters['rain_cbh_adj'].maximum

# %%
idx = 64203 * 11
aa[idx:idx+11]

# %%
pdb.parameters['tstorm_mo'].data.shape

# %%
pdb.parameters.check()

# %%
minidom.parseString(xmlET.tostring(pdb.xml_global_parameters)).toprettyxml(indent='    ')

# %%
pdb.write_parameters_xml('/Users/pnorton/tmp')

# %%
tuple(pdb.parameters['tstorm_mo'].dimensions.keys())

# %% [markdown]
# ### Write parameters to netCDF file

# %%
pdb.write_netcdf('crap.nc')

# %%
print(pdb.parameters['hru_area'].minimum)
type(pdb.parameters['hru_area'].maximum)

# %% [markdown]
# #### What's the longest number of characters in the poi_gage_id array?

# %%
len(max(pdb.parameters['poi_gage_id'].data, key=len))

# %%
aa = pdb.parameters['poi_gage_id'].dimensions.keys().extend(['poi_gage_id_chars'])
print(aa)


# %%
aa = pdb.parameters['tmax_cbh_adj'].dimensions.keys()
print(aa)
aa.reverse()
print(aa)

# %%
pdb.parameters['tmax_cbh_adj'].data

# %%
print(pdb.parameters['x_coef'])

# %%
import netCDF4 as nc

# %%
nc.default_fillvals

# %%
from future.utils import iteritems

cparam = {'name': 'test',
          'datatype': 'float',
          'dimensions': {'nhru': {'size': 445,
                                  'position': 1},
                         'nmonth': {'size': 12,
                                    'position': 2}},
          'data': []}

# %%
for xx, yy in iteritems(cparam):
    print(xx, yy)

# %%
for xx, yy in iteritems(cparam['dimensions']):
    print(xx, yy['size'])

# %%
print(pdb.parameters['radmax'].dimensions)

# %%
for xx, yy in iteritems(pdb.parameters['radmax'].dimensions):
    print(xx, yy.size)

# %%
import os.path

# %%
os.path.split('/Users/pnorton/Projects/National_Hydrology_Model/notebooks/03_pyPRMS_prototype/crap.nc')

# %%
os.path.splitext('/Users/pnorton/Projects/National_Hydrology_Model/notebooks/03_pyPRMS_prototype/crap.nc')

# %%
os.path.basename('/Users/pnorton/Projects/National_Hydrology_Model/notebooks/03_pyPRMS_prototype/crap.nc')

# %%
# os.path.dirname('/Users/pnorton/Projects/National_Hydrology_Model/notebooks/03_pyPRMS_prototype/crap.nc')
os.path.dirname('/Users/pnorton/Projects/National_Hydrology_Model/notebooks/03_pyPRMS_prototype')

# %%
os.path.isfile('/Users/pnorton/Projects/National_Hydrology_Model/notebooks/03_pyPRMS_prototype/crap.nc')

# %%
