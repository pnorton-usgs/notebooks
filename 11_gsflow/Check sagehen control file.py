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
from pyPRMS.ControlFile import ControlFile

# %%
base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/sagehen/prms6'
control_dir = f'{base_dir}'
control_file = f'{control_dir}/prms6.control'

# %%
ctl = ControlFile(control_file)

# %%
ctl.control_variables.keys()

# %%
print(ctl.get('windspeed_day'))

# %%
modules_used = ctl.modules.values()
print(modules_used)

# %%
ctl.modules.items()

# %%
ctl.get('temp_module').values

# %%

# %%
ctl.write(f'{control_dir}/prms6.control')

# %%
