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
#     display_name: Python [conda env:bandit_py3]
#     language: python
#     name: conda-env-bandit_py3-py
# ---

# %%
from future.utils import iteritems

from pyPRMS.ValidParams import ValidParams
from pyPRMS.NhmParamDb_AK import NhmParamDb_AK

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit_AK_test/akhmparamdb'

src_db = NhmParamDb_AK(workdir)
mpdb = ValidParams()

# %%
mpdb.keys()

# %%
for kk, vv in iteritems(mpdb):
#     print(kk, vv.dimensions.keys())
#     print(src_db.parameters[kk].dimensions)
    
    if src_db.parameters.exists(kk):
        if set(src_db.parameters[kk].dimensions.keys()).issubset(set(vv.dimensions.keys())):
            pass
#             print(kk, 'OK')
        else:
            print('{} dimensions should be {}, NOT {}'.format(kk, list(vv.dimensions.keys()), list(src_db.parameters[kk].dimensions.keys())))
            
            if 'one' in vv.dimensions.keys():
                print('  --reset to scalar')
                for kx in list(src_db.parameters[kk].dimensions.keys()):
                    src_db.parameters[kk].dimensions.remove(kx)
                src_db.parameters[kk].dimensions.add('one', size=1)
                crap = src_db.parameters[kk].data[0]
                src_db.parameters[kk].data = [crap]
                
            if 'nsegment' in vv.dimensions.keys():
                print('  --reset to nsegment')
                orig_size = src_db.parameters[kk].dimensions['nhru'].size
                src_db.parameters[kk].dimensions.remove('nhru')
                src_db.parameters[kk].dimensions.add('nsegment', size=orig_size)
                
            if 'ndeplval' in vv.dimensions.keys():
                print(' --reset to ndeplval')
                orig_nhru = src_db.parameters[kk].dimensions['nhru'].size
                for kx in list(src_db.parameters[kk].dimensions.keys()):
                    src_db.parameters[kk].dimensions.remove(kx)
                
                src_db.parameters[kk].dimensions.add('ndeplval', size=orig_nhru*11)
                
#             print(kk, vv.dimensions.keys())
#             print(src_db.parameters[kk].dimensions.keys())

# %%
print(src_db.parameters['snarea_curve'])

# %%
src_db.parameters.check()

# %%
print(src_db.parameters['hru_deplcrv'].data)

# %%
print(src_db.parameters['snarea_curve'].as_dataframe.values.reshape((-1, 11)).shape[0])

# %%
