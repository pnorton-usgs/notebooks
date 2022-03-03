# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
# ---

# %%
import prms_control_lib as ctl

# %%
workdir = "/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t1/06267400/runs/TST1/-0001"
ctlfile = "%s/control/default.control" % workdir

# %%
reload(ctl)
cobj = ctl.control(ctlfile)

# %%
cobj.get_var('data_file')

# %%
