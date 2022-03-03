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
#     display_name: Python 2
#     language: python
#     name: python2
# ---

# %%
from __future__ import (absolute_import, division, print_function)
from future.utils import iteritems, itervalues

import prms_cfg
reload(prms_cfg)

# %%
workdir = "/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/pipestem_testfiles"

config = prms_cfg.cfg('%s/basin.cfg' % workdir, expand_vars=False)

# %%
filelist = []

for vv in itervalues(config.get_value('objfcn')):
    print(vv['obs_file'], vv['sd_file'])
    if vv['obs_file']:
        filelist.append(vv['obs_file'])
    if vv['sd_file']:
        filelist.append(vv['sd_file'])
        
print(filelist)
    

# %%
ff = [vv['obs_file'] for vv in itervalues(config.get_value('objfcn'))]
print(ff)

# %%
print config.param_range_file

# %%
config.write_config('crap.cfg')

# %%
