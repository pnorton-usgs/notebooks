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
#     display_name: Python [conda env:pyGDP]
#     language: python
#     name: conda-env-pyGDP-py
# ---

# %%

# %%
from owslib.csw import CatalogueServiceWeb

# %%
csw = CatalogueServiceWeb('')

# %%
csw.identification.type

# %%
[op.name for op in csw.operations]

# %%
csw.getdomain('GetRecords.resultType')

# %%
