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
import pyPRMS.NhmParamDb as nhm
import pandas as pd
# reload(nhm)

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/nhmparamdb'

pdb = nhm.NhmParamDb(workdir)

# %%
list(pdb.available_parameters)

# %%
pdb.hru_nhm_to_region

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
print(pdb.parameters)

# %%
print(pdb.dimensions)

# %%
print(pdb.parameters.get('gwflow_coef'))

# %%
print(pdb.parameters.get('gwflow_coef').data)

# %%
print(pdb.parameters['snarea_curve'])

# %%
print(pdb.parameters['hru_deplcrv'].data)

# %%
bb = pdb.parameters.get_DataFrame('hru_deplcrv')
print(bb.head())

# %%
print(pdb.parameters['tstorm_mo'].dimensions.tostructure())

# %%
print(pdb.parameters['tstorm_mo'].tostructure())

# %%
print(pdb.parameters['snarea_curve'].tostructure())

# %%
pdb.parameters['tstorm_mo'].get_dimsize_by_index(1)

# %%
print(pdb.parameters['tstorm_mo'].tolist())

# %%
pdb.parameters.get_DataFrame('snarea_curve')

# %%
aa = pdb.parameters['snarea_curve'].data

# %%
idx = 64203 * 11
aa[idx:idx+11]

# %%
pdb.parameters['tstorm_mo'].data.size / 12

# %%
pdb.parameters.check()

# %%
x = []

# %%
x.insert(0, 'b')
x.insert(1, 'a')


# %%
x

# %%
