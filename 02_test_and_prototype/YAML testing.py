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
import collections
import yaml
import pprint

# %%
workdir = "/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t4"
cfgfile = '%s/basin_yaml.cfg' % workdir
outfile = '%s/basin_yaml_out.cfg' % workdir

# %%

for project in yaml.load_all(open(cfgfile)):
    print project
    #pprint.pprint(project)

# %%
#collections.OrderedDict(cfg)
cfg = {'start_date_model': '1980-10-1',
       'start_date': '1982-10-1',
       'end_date': '2010-9-30',
       'nstart': 200, 
       'nsets': 28,
       'nparam': 14, 
       'ntests': 2,
       'base_dir': '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing',
       'base_calib_dir': '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t4', 
       'template_dir': '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/gage_master', 
       'runs_sub': 'runs',
       'prms_control_sub': None, 
       'prms_input_sub': None, 
       'prms_output_sub': None, 
       'basins_xref_file': '%(template_dir)s/lookup.csv',
       'basins_file': '%(base_calib_dir)s/basin_list.txt',
       'log_file': 'optim_log.log', 
       'default_param_list_file': 'default_param_ranges.txt', 
       'param_range_file': 'param_limits.txt', 
       'test_func_margin_file': 'test_func_margins.txt', 
       'prms_control_file': 'daymet.control', 
       'prms_input_file': 'daymet.params',
       'cmd_cp': '/bin/cp', 
       'cmd_ln': '/bin/ln', 
       'cmd_mkdir': '/bin/mkdir', 
       'cmd_prms': '%(base_dir)s/bin/prms',
       'cmd_mocom': '%(base_dir)s/bin/mocom', 
       'calib_run': '/Users/pnorton/Projects/National_Hydrology_Model/code/prms_calib/prms_post_mocom.py', 
       'stats_script': '/Users/pnorton/Projects/National_Hydrology_Model/code/prms_calib/prms_objfcn.py', 
       'plot_script': 'plot_day.ALL.csh',

       'objfcn': {'OF1': {'of_stat': 'NRMSE', 'of_intv': 'daily', 
                          'obs_type': 'value', 'obs_intv': 'daily', 'obs_file': 'filename', 
                          'sdval':6, 'sdfile': 'filename', 
                          'obs_var': 'runoff', 'sim_var': 'basin_cfs'},
                  'OF2': {'of_stat': 'NRMSE', 'of_intv': 'daily', 
                          'obs_type': 'value', 'obs_intv': 'daily', 'obs_file': 'filename', 
                          'sdval':7, 'sdfile': 'filename', 
                          'obs_var': 'runoff', 'sim_var': 'basin_cfs'},
                  'OF3': {'of_stat': 'NRMSE', 'of_intv': 'daily', 
                          'obs_type': 'value', 'obs_intv': 'daily', 'obs_file': 'filename', 
                          'sdval': None, 'sdfile': None, 
                          'obs_var': 'runoff', 'sim_var': 'basin_cfs'}},

       'of_link': {'OFa': {'of_names': ['OF1', 'OF2'], 'of_wgts': [0.8, 0.2]},
                   'OFb': {'of_names': ['OF3'], 'of_wgts': [1.0]}}}
print cfg

# %%
print yaml.dump(cfg, open(outfile, 'w'))

# %%
for project in yaml.load_all(open(outfile)):
    print project

incfg = yaml.load(open(outfile))

# %%
incfg['of_link']['OFa']

# %%
incfg.of_link

# %%
