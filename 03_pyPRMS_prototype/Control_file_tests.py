# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import datetime
import os

from pyPRMS.ControlFile import ControlFile
# reload(pyPRMS.ControlFile)

# %%
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20190418_hw_6000'
# workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20210312_red_river'
# control_file = '{}/control.default.bandit'.format(workdir)
# control_file = '/Users/pnorton/notes/bandit_default_stuff/control.default'
# control_file = '/Users/pnorton/tmp/prms_chk/4551/control.default.bandit'
control_file = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/control.default'

# %%
ctl = ControlFile(control_file, verbose=True)

# %% [markdown] tags=[]
# ### List modules defined in control file

# %%
ctl.modules

# %%
ctl.control_variables['initial_deltat'].values

# %%
print(ctl.get('start_time'))

# %%
ss = ctl.get('start_time').values

# %%
stdate = datetime.datetime(*ss)

# %%
print(stdate)
ww = ctl.get('prms_warmup').values
print(f'prms_warmup: {ww}')


# %%

# %%
ctl.write('crap.control')

# %%

# %%
nhru_prefix = ctl.get('nhruOutBaseFileName').values
print(nhru_prefix)

# %%
os.path.normpath(os.path.join(f'{workdir}/{nhru_prefix}bob.csv'))

# %%
ctl.control_variables.keys()

# %%
# ctl.remove('aniOutVar_names')

# %%
ctl.get('stream_temp_shade_flag').values

# %%
ctl.get('end_time').values

# %%
ctl.write('control.crap')

# %%
ctl.get('snareathresh_dynamic').values

# %%
ctl.get('dyn_transp_flag').default

# %%
ctl.get('prms_warmup').size

# %%
ctl.get('dyn_covden_flag').valid_values

# %%
ctl.get('dyn_covtype_flag').values

# %%
ctl.get('dyn_covtype_flag').value_repr

# %%
ctl.get('prms_warmup').value_repr

# %%
ctl.get('dyn_covtype_flag').associated_values

# %%
ctl.dynamic_parameters

# %%
ctl.get('dyn_dprst_flag').associated_values

# %%
ctl.get('dyn_dprst_flag').valid_values

# %%
ctl.has_dynamic_parameters

# %%
ctl.dynamic_parameters

# %%
ctl.get('dyn_intcp_flag').values

# %%
print(ctl.get('dyn_intcp_flag'))

# %%
print(ctl.get('dyn_covtype_flag'))

# %%
ctl.get('start_time').associated_values

# %%
ctl.get('nsegmentOutVar_names').values

# %%
