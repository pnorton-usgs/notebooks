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
import datetime

# %%
nobs = 8270
dummy_data = '-999 ' * nobs

hdl = open('sf_data', 'w')

hdl.write('Created by Bandit\n')
hdl.write('/////////////////////////////////////////////////////////////////////////\n')
hdl.write('// Station IDs for runoff:\n')
hdl.write('// ID\n')
hdl.write('/////////////////////////////////////////////////////////////////////////\n')
hdl.write('// Unit: runoff = cfs\n')
hdl.write('/////////////////////////////////////////////////////////////////////////\n')
hdl.write('runoff {}\n'.format(nobs))
hdl.write('#########################################################\n')

start = datetime.datetime(1930,1,1)
end = datetime.datetime(2018,12,31)
date_generated = [start + datetime.timedelta(days=x) for x in range(0, (end-start).days+1)]

for date in date_generated:
    hdl.write(date.strftime("%Y %m %d") + ' 0 0 0 ' + dummy_data + '\n')
    
hdl.close()

# %%
start = datetime.datetime.strptime("21-06-2014", "%d-%m-%Y")
end = datetime.datetime.strptime("07-07-2014", "%d-%m-%Y")
date_generated = [start + datetime.timedelta(days=x) for x in range(0, (end-start).days)]

for date in date_generated:
    print(date.strftime("%d-%m-%Y"))

# %%
