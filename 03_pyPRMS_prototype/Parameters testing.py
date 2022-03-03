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
#     display_name: Python [default]
#     language: python
#     name: python2
# ---

# %%
import pyPRMS.Parameters as prm
reload(prm)

import xml.etree.ElementTree as xmlET

# %%
params = prm.Parameters()

# %%
print(params)

# %%
params.add('joe')

# %%
print(params['joe'])

# %%
params.get('joe').datatype = 1
params.get('joe').dimensions.add('nmonths', 12)
print(params['joe'])

# %%
# params.get_DataFrame('joe')
params['joe'].as_dataframe

# %%
params.get('joe').data = [1,2,3,4,5,6,7,8,9,10,11,12]
print(params['joe'])

# %%
# params.get_DataFrame('joe')
params['joe'].as_dataframe

# %%
params.exists('joe')

# %%

# %%
# Get list of dimensions defined for a parameter
params.get('joe').dimensions.keys()

# %%
xmlET.dump(params['joe'].xml)

# %%
params['joe'].paramDb

# %%
with open('crap.csv', 'w') as f:
    f.write(params['joe'].paramDb)

# %%
params.write_paramdb('testout')

# %%
