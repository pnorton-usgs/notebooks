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
import pymake

# %% [markdown]
# The pymake library is the library that Chris and Joe created to build the MODFLOW6 code. One of the functions provides the dependency graph generation. I believe pymake requires `graphviz`, `pydot`, and `pydotplus`. I did the following to install pymake in one of my conda environments.
#
# ```bash
# conda activate <some_env>
#
# conda install pydot pydotplus graphviz
#
# cd <dir_to_clone_pymake_into>
# git clone https://github.com/modflowpy/pymake.git
#
# pip install pymake
# ```
#
# Below are examples of generating dependency graphs from source code. The output_dir should exist prior to running the make_plots routine.

# %%
sourcecode_dir = '/Users/pnorton/Projects/National_Hydrology_Model/src/coretran/src'
output_dir = '/Users/pnorton/tmp/coretran_dot'

pymake.make_plots(sourcecode_dir, output_dir, include_subdir=True, extension='.pdf')

# %%
sourcecode_dir = '/Users/pnorton/src/modflow6/src'
output_dir = '/Users/pnorton/tmp/modflow_dot'

pymake.make_plots(sourcecode_dir, output_dir, include_subdir=True, extension='.pdf')

# %%
sourcecode_dir = '/Users/pnorton/Projects/National_Hydrology_Model/code/prms5'
output_dir = '/Users/pnorton/tmp/prms5_dot'

pymake.make_plots(sourcecode_dir, output_dir, include_subdir=True, extension='.pdf')

# %%
sourcecode_dir = '/Users/pnorton/Projects/National_Hydrology_Model/src/prms'
output_dir = '/Users/pnorton/tmp/prms6_dot'

pymake.make_plots(sourcecode_dir, output_dir, include_subdir=True, extension='.pdf')

# %%
sourcecode_dir = '/Users/pnorton/Projects/National_Hydrology_Model/src/gsflow_v2/GSFLOW/src'
output_dir = '/Users/pnorton/tmp/gsflow2_dot'

pymake.make_plots(sourcecode_dir, output_dir, include_subdir=True, extension='.pdf')

# %%
