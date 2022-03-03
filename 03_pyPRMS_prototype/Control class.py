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
from pyPRMS.Control import Control

# %%
ctl = Control()

# %%
ctl.control_variables.keys()

# %%
ctl['start_time'].values

# %%
ctl.has_dynamic_parameters

# %%
ctl.header

# %%
ctl.modules

# %%
