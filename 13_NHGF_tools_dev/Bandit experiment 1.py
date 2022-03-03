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
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
from pyPRMS.ControlFile import ControlFile
from pyPRMS.ParamDb import ParamDb

from Bandit.bandit_helpers import parse_gages, set_date, subset_stream_network

# %%
ctl_filename = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/control.default'
paramdb_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS'

# %%
# Read a control file
ctl = ControlFile(ctl_filename)

# %%
# %%time
pdb = ParamDb(paramdb_dir)

# %%
params_before = list(pdb.parameters.keys())
print(f'Number of parameters: {len(params_before)}')

# %%
# Wrap in conditional based on no_filter_params command line argument
pdb.control = ctl
pdb.reduce_by_modules()

# %%
params_after = list(pdb.parameters.keys())
print(f'Number of parameters: {len(params_after)}')

# %%
dag_ds = pdb.parameters.stream_network(tosegment='tosegment_nhm', seg_id='nhm_seg')

# %%
dsmost_seg = [30119]
uscutoff_seg = []
dag_ds_subset = subset_stream_network(dag_ds, uscutoff_seg, dsmost_seg)

# %%
new_nhm_seg = [ee[0] for ee in dag_ds_subset.edges]

# %%
new_nhm_seg

# %%
# get_subset(self, name: str, global_ids: List[int])
pdb.parameters.get_subset('K_coef', new_nhm_seg)

# %%
aa = pdb.parameters.hru_to_seg

# %%
# newdict = {k: testdict[k] for k in keep}
bb = {kk: aa[kk] for kk in new_nhm_seg}

# %%
bb

# %%
noroute = [30120, 30122]
cc = {kk: aa[kk] for kk in noroute if aa[kk] == 0}

# %%
cc

# %%

# %%
