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
from sets import Set
import numpy as np
import prms_lib as prms
import os
reload(prms)

#workdir = "/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t3/06191500/runs/2015-03-03_1821/03250"
workdir = "/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/pipestem_testfiles"
statvarDir = workdir
paramDir = workdir
controlDir = workdir

# %% [markdown]
# #### Statvar files

# %%
sv = prms.statvar("%s/daymet.statvar" % statvarDir)

# Get the filename and the vars
print sv.filename
df = sv.data
df.drop(['rec'], axis=1, inplace=True)

print sv.vars

# %%
#df.index.to_pydatetime()
#df['runoff']

plotvars = df.columns.tolist()
print plotvars
plotvars.remove('runoff')
print plotvars

# %% [markdown]
# #### Control file

# %%
# reload(prms)

controlfile = 'daymet.control'

ctl = prms.control('%s/%s' % (workdir, controlfile))

ctl.clear_parameter_group('statVar')
ctl.clear_parameter_group('aniOut')
ctl.clear_parameter_group('dispVar')
ctl.clear_parameter_group('mapOut')
ctl.clear_parameter_group('nhruOut')


# strip path information from file parameters
chg_paths = ['model_output_file', 'data_file', 'ani_output_file', 'var_init_file',
             'stat_var_file', 'var_save_file', 'stats_output_file', 'tmax_day', 
             'tmin_day', 'precip_day', 'param_file', 'csv_output_file']

for pp in chg_paths:
    ctl.replace_values(pp, os.path.basename(ctl.get_values(pp)))
# --------------------------------------------

ctl.modules

# ctl.write_control_file('%s/control/daymet.control.new' % controlDir)

# %% [markdown]
# #### Load database of all available parameters by module

# %%
pdb = prms.param_db('/Users/pnorton/Projects/National_Hydrology_Model/regions/all_modules')

# Get parameters by module name
mypdb = pdb.module_params(ctl.modules)

# Get subset of parameters required by supplied modules
mypdb_subset = pdb.get_param_subset(ctl.modules)

# %%
param.check_all_module_vars(mypdb)

# %% [markdown]
# ### Input parameter files

# %%
reload(prms)
#paramDir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t3/06191500/runs/2015-03-02_0813/00001'
# paramDir = '/Users/pnorton/Projects/National_Hydrology_Model/regions/r10U'
param = prms.parameters('%s/daymet.params' % workdir)

# %%
# Expand input parameters, adding, expanding, and/or converting parameters as needed
param.expand_params(mypdb_subset)

# %%
# Check that all required parameters are included
param.check_all_module_vars(mypdb)

# %%
# basin_*_frost are only required when model_mode=FROST
param.remove_param('basin_fall_frost')
param.remove_param('basin_spring_frost')
param.remove_param('fall_frost')
param.remove_param('spring_frost')

# prms_warmup was moved into control file
param.remove_param('prms_warmup')

# %%
srcdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/Snow_Depletion_Curves/20160331'
files = ['r10U_SDCs.param']

for ff in files:
    param.add_param_from_file('%s/%s' % (srcdir, ff))

# %%
param.write_param_file('%s/daymet.params.expanded' % workdir)

# %%
print param.dimensions
param.add_dimension('joe', 900)
print param.dimensions
param.resize_dim('joe', '1000')
print param.dimensions

# %%
aa = [1,2,3,4,5,6,7,8,9,10]
bb = np.array(aa).reshape((2,5))
cc = np.array(aa).reshape([2,5], order='F')
print aa
print bb
print cc
print bb[0,:]
print cc[0,:]

# %%

# %%

# %%
tst = param.get_var('snarea_curve')

vals = tst['values'].reshape((11,-1), order='F')
print vals.shape

# %%
param.pull_hru2(2280, 'crap.param')

# %%
tst = param.get_var('snarea_curve')
print tst

# %% [markdown]
# ### Load full database of available parameters by module

# %%
# Adding parameters from file(s)
srcdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/regions/r10U'
files = ['carea_min.param', 'dprst_frac.param', 'imperv_frac.param', 
         'sro_to_dprst_imperv.param', 'sro_to_dprst_perv.param', 
         'tmax_allsnow.param']

for ff in files:
    param.add_param_from_file('%s/%s' % (srcdir, ff))


# %%
# Write out the new input parameter file
param.write_param_file('%s/daymet.params.expanded' % paramDir)

# %%

# %% [markdown]
# ### Streamflow data file

# %%
import re
reload(prms)
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/r10U_streamflow'
streamflowData = '%s/1960-2010.data' % workdir

sf = prms.streamflow(streamflowData)

# %%
print 'metaheader:', sf.metaheader
print 'units:', sf.units
print 'types:', sf.types

# %%
sf.selectByStation(['06191000','06191500','06192500'])
#sf.selectByStation(['06477500'])
#sf.clearSelectedStations()


# %%
sf.writeSelectedStations('crap.data')

# %%
import pandas as pd
b = sf.data.reset_index()
b.head()

# %%
c = pd.DatetimeIndex(b['thedate'])
b['year'] = c.year
b['month'] = c.month
b['day'] = c.day
b['hour'] = c.hour
b['minute'] = c.minute
b['second'] = c.second
b.head()

# %%
import csv
b.to_csv('test.txt', index=False, header=False, date_format='%Y %m %d %H %M %S', sep=' ', quotechar='\0')


# %% [markdown]
# ### Test read_gdp()

# %%
reload(prms)

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets'
filename = 'snow_only_SNODASPRCP_Daily_r09_byHRU_2004-01-01_2014-12-31.csv'

gdp_test = prms.read_gdp('%s/%s' % (workdir, filename), missing_val=400.)

# %%
gdp_test.head(10)

# %% [markdown]
# ### Test read_cbh()

# %%
reload(prms)
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets'
filename = 'daymet_1980_2010_prcp.cbh'

cbh_test = prms.read_cbh('%s/%s' % (workdir, filename))



# %%
cbh_test.info()

# %%
print '%7.2f' % 234.4

# %%
aa = ['aa', 'bbb', 'c', 'dd']
len(max(aa, key=len))

# %%
