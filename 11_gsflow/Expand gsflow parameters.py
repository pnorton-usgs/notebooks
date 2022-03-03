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

# %% [markdown]
# GSFLOW parameter files sometimes contain 'collapsed' parameter values (e.g. where a duplicated value is specified as N*value). This routine will expand those values to a single value per line.

# %%
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/sagehen/input/prms'
filename = f'{work_dir}/prms_orig.params'

new_file = f'{work_dir}/prms_orig_v2.params'

# %%
fhdl = open(filename, 'r')
rawdata = fhdl.read().splitlines()
fhdl.close()
it = iter(rawdata)

out_fhdl = open(new_file, 'w')

for line in it:
    if line in ['** Dimensions **', '** Parameters **', ]:
        out_fhdl.write(f'{line}\n')
    else:
        aa = line.split('*')
        
        if len(aa) == 2:
            for ex in range(int(aa[0])):
                out_fhdl.write(f'{aa[1]}\n')
        else:
            out_fhdl.write(f'{line}\n')

out_fhdl.close()

# %%

# %%
