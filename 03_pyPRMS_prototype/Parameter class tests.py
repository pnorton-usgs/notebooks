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
from pyPRMS.Parameters import Parameter

# %%
param = Parameter('joebob')

# %%
print(param)

# %%
param.datatype = 1
print(param)

# %%
param.datatype = 5
print(param)

# %%
param.add_dimension('nmonths', 12)
param.add_dimension('nhru', 4)
print(param)

# %%
