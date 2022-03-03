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
from pyPRMS.ParameterFile import ParameterFile
from pyPRMS.ControlFile import ControlFile
from pyPRMS.ValidParams import ValidParams

# %%
base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms6/sagehen'
param_dir = f'{base_dir}/prms5'
param_file = f'{param_dir}/prms_grid_v2_prms5.param'

control_dir = f'{base_dir}/prms6'
control_file = f'{control_dir}/prms6.control'


# %%
# Load master list of valid parameters
vpdb = ValidParams()

# %%
ctl = ControlFile(control_file)

# %%
pf = ParameterFile(param_file, verbose=True)

# %%
# print(pf.parameters.keys())

# %%
pf.parameters.check()

# %%
modules_used = ctl.modules.values()

# %%
# Get the list of parameters that are required for the selected modules
required_params = vpdb.get_params_for_modules(modules=modules_used)
print('Number of required parameters: {}'.format(len(required_params)))

# %%
print(f'Number of parameters in file: {len(pf.parameters.keys())}')

# %%
print('Parameters required by modules but missing from parameter file')
print(set(required_params).difference(set(pf.parameters.keys())))

# %%
print('Parameters in parameter file that are not needed')
print(set(pf.parameters.keys()).difference(set(required_params)))

# %%
print(required_params)

# %%
for rp in required_params:
    if not pf.parameters.exists(rp):
        print(f'{rp} is missing')

# %%
pf.degenerate_parameters()

# %%
# Remove the parameters that are not needed for the selected modules in the control file
pf.reduce_by_modules(ctl)

# %%

# %%
expand_list = ['fastcoef_lin', 'soil_rechr_init_frac', 'smidx_coef', 
               'transp_tmax', 'slowcoef_lin', 'radadj_slope', 'radmax',
               'fastcoef_sq', 'ssr2gw_exp', 'freeh2o_cap', 'tmax_allsnow',
               'dday_slope', 'smidx_exp', 'dday_intcp', 'transp_end',
               'transp_beg', 'carea_max', 'adjmix_rain', 'potet_sublim', 
               'ssstor_init_frac', 'radadj_intcp', 'cecn_coef',
               'ppt_rad_adj', 'soil2gw_max', 'epan_coef', 'melt_look', 
               'snowinfil_max', 'pref_flow_den', 'jh_coef',
               'tmax_index', 'melt_force', 'slowcoef_sq',
               'radj_wppt', 'radj_sppt', 'tstorm_mo', 'gwstor_min',
               'gwflow_coef', 'gwsink_coef', 'gwstor_init', 'emis_noppt', 
               'snowpack_init', 'snow_cbh_adj', 'rain_cbh_adj',
               'tmax_cbh_adj', 'tmin_cbh_adj', 'tmax_allrain_offset', 'segment_type',
               'obsin_segment', 'obsout_segment', 'segment_flow_init']
expand_list.sort

for pp in expand_list:
    print(pp)
    pf.expand_parameter(pp)

# %%
pf.expand_parameter('hru_deplcrv')

# %%
pf.write_netcdf(f'{base_dir}/prms6/prms_grid_v3.nc')
# pf.write_netcdf(f'{param_dir}/prms_grid_v3.nc')

# %%
pf.write_parameter_file(f'{param_dir}/crap.param')

# %%
print(pf.parameters['snarea_curve'])

# %%
pf.parameters['snarea_curve'].data.shape

# %%
print(pf.dimensions)

# %%
pf.dimensions.get('ndepl').size = 1
pf.dimensions.get('ndeplval').size = 1408

# %%
type(pf.parameters['obsin_segment'].maximum)

# %%
print(pf.parameters['hru_elev'])

# %%
pf.parameters['hru_elev'].stats()

# %%
pf.parameters['hru_elev'].data[95]

# %%
