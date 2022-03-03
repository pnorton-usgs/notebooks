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
import pyPRMS.Dimensions as dmn
from collections import OrderedDict
import xml.etree.ElementTree as xmlET
# reload(dmn)

# %%
dims = dmn.Dimensions()

# %%
print(dims)

# %%
# Invalid dimension name, should produce a warning
dims.add('george', 1)

# %%

# %%
# Valid dimension name, should be added
dims.add('nmonths', 12)

# %%
print(dims)

# %%
dims.tostructure()

# %%
dims['nmonths'].size

# %%
xmlET.dump(dims.xml)

# %%
dims = dmn.Dimensions()

# %%
dims.add('nmonths', 12)
dims.add('one', 1)

# %%
print(dims)

# %%
dims.exists('two')

# %%
xmlET.dump(dims.xml)

# %%
# If nmonths already exists this is silently ignored
dims.add('nmonths', 13)

# %%
print(dims)

# %%
dims.remove('one')
print(dims)

# %%
dims['one'] += 1

# %%
dmn.Dimensions.Dimension.name

# %%
