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
#     display_name: Python [default]
#     language: python
#     name: python2
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import pyPRMS.ParamDb as pdb
# reload(nhm)

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/paramDb'

params = pdb.ParamDb(workdir)

# %%
params.available_parameters

# %%
# print(pdb.parameters['tstorm_mo'].dimensions.get_position('nmonths'))
print(params.parameters)

# %%
print(params.dimensions)

# %%
print(params.parameters.get('gwflow_coef'))

# %%
print(params.parameters.get('gwflow_coef').data)

# %%
print(params.parameters['snarea_curve'])

# %%
print(params.parameters['hru_deplcrv'].data)

# %%
bb = params.parameters.get_DataFrame('hru_deplcrv')
print(bb.head())

# %%
print(params.parameters['tstorm_mo'].dimensions.tostructure())

# %%
print(params.parameters['tstorm_mo'].tostructure())

# %%
print(params.parameters['snarea_curve'].tostructure())

# %%
params.parameters['tstorm_mo'].get_dimsize_by_index(1)

# %%
print(params.parameters['tstorm_mo'].tolist())

# %%
params.parameters.get_DataFrame('snarea_curve')

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
