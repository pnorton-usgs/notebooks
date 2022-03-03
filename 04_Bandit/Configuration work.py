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
import Bandit.bandit_cfg as bc

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit'
config_file = f'{workdir}/bandit_v11.cfg'

# %%
config = bc.Cfg(config_file)

# %%
config.cbh_dir

# %%
bb = config.cbh_var_map

# %%
bb['tmax']

# %%
list(bb.keys())

# %%
config.exists('cbh_var_map')

# %%
for cbhvar, cfv in config.cbh_var_map.items():
    print(cbhvar, cfv)

# %%
for cfv in config.cbh_var_map.values():
    print(cfv)

# %%
print(config)

# %%

# %%
import sys
from ruamel.yaml import YAML

# %%
inp = """\
- &CENTER {x: 1, y: 2}
- &LEFT {x: 0, y: 2}
- &BIG {r: 10}
- &SMALL {r: 1}
# All the following maps are equal:
# Explicit keys
- x: 1
  y: 2
  r: 10
  label: center/big
# Merge one map
- <<: *CENTER
  r: 10
  label: center/big
# Merge multiple maps
- <<: [*CENTER, *BIG]
  label: center/big
# Override
- <<: [*BIG, *LEFT, *SMALL]
  x: 1
  label: center/big
"""

yaml = YAML()
data = yaml.load(inp)
assert data[7]['y'] == 2

# %%
for xx in data:
    print(xx)

# %%
yaml.dump(data, sys.stdout)

# %%

# %%
