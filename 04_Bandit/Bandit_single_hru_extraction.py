# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %%
import numpy as np

from pyPRMS.ControlFile import ControlFile
from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterSet import ParameterSet
from pyPRMS.ValidParams import ValidParams

from pyPRMS.constants import HRU_DIMS

# %% [markdown]
# ### Exploratory work for extracting a single HRU from the National Hydrologic Model

# %%
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS'

control_filename = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHRU_sample/HRU3323_run2/control_single'

out_param_file = '/Users/pnorton/tmp/single_hru_test.param'

# %%
# %%time
# Load the control file
ctl = ControlFile(control_filename)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load master list of valid parameters
vpdb = ValidParams()

# Build list of parameters required for the selected control file modules
required_params = vpdb.get_params_for_modules(modules=list(ctl.modules.values()))

pdb = ParamDb(paramdb_dir=work_dir, verbose=True, verify=True)

nhm_params = pdb.parameters
nhm_global_dimensions = pdb.dimensions

# %%
nhm_id = nhm_params.get('nhm_id').tolist()
nhm_id_to_idx = nhm_params.get('nhm_id').index_map

# %%
# HRU to extract: 83897

# %%
hru_order_subset = [83897]

# %%
print(hru_order_subset)
# nhm_params.get('hru_segment_nhm').tolist()[hru_order_subset]

# %%
params = list(nhm_params.keys())
params.sort()

# %%
# Use hru_order_subset to pull selected indices for parameters with nhru dimensions
# hru_order_subset contains the in-order indices for the subset of hru_segments
# toseg_idx contains the in-order indices for the subset of tosegments
# --------------------------------------------------------------------------

# ==========================================================================
# ==========================================================================
# Get subset of hru_deplcrv using hru_order
# A single snarea_curve can be referenced by multiple HRUs
hru_deplcrv_subset = nhm_params.get_subset('hru_deplcrv', hru_order_subset)

uniq_deplcrv = list(set(hru_deplcrv_subset))
uniq_deplcrv0 = [xx - 1 for xx in uniq_deplcrv]

uniq_dict = {}
for ii, xx in enumerate(uniq_deplcrv):
    uniq_dict[xx] = ii + 1

# Create new hru_deplcrv and renumber
new_hru_deplcrv = [uniq_dict[xx] for xx in hru_deplcrv_subset]

# %%

# %%

# %%

# %%
dims = {}
for kk in nhm_global_dimensions.values():
    dims[kk.name] = kk.size
    
# Remove dimensions that are not used in a single-HRU model
del dims['nsegment']
del dims['npoigages']

# Resize dimensions to the model subset
crap_dims = dims.copy()  # need a copy since we modify dims
for dd, dv in crap_dims.items():
    # dimensions 'nmonths' and 'one' are never changed
    if dd in HRU_DIMS:
        dims[dd] = 1
    elif dd == 'ndeplval':
        dims[dd] = 11
        dims['ndepl'] = 1

# %%
dims

# %%
new_hru_deplcrv = [1]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Build a ParameterSet for output
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
new_ps = ParameterSet()

for dd, dv in dims.items():
    new_ps.dimensions.add(dd, dv)

# Add nobs (missing from paramdb but needed)
new_ps.dimensions.add('nobs', 1)

new_params = list(required_params)

# WARNING: 2019-04-23 PAN
#          Very hacky way to remove parameters that shouldn't always get
#          included. Need to figure out a better way.
check_list = ['basin_solsta', 'gvr_hru_id', 'hru_solsta', 'humidity_percent',
              'irr_type', 'obsout_segment', 'rad_conv', 'rain_code', 'hru_lon']

for xx in check_list:
    if xx in new_params:
        if xx in ['basin_solsta', 'hru_solsta', 'rad_conv']:
            if not new_ps.dimensions.exists('nsol'):
                new_params.remove(xx)
            elif new_ps.dimensions.get('nsol') == 0:
                new_params.remove(xx)
        elif xx == 'humidity_percent':
            if not new_ps.dimensions.exists('nhumid'):
                new_params.remove(xx)
            elif new_ps.dimensions.get('nhumid') == 0:
                new_params.remove(xx)
        elif xx == 'irr_type':
            if not new_ps.dimensions.exists('nwateruse'):
                new_params.remove(xx)
            elif new_ps.dimensions.get('nwateruse') == 0:
                new_params.remove(xx)
        elif xx == 'gvr_hru_id':
            if ctl.get('mapOutON_OFF').values == 0:
                new_params.remove(xx)
        elif xx in ['hru_lat', 'hru_lon', ]:
            if not nhm_params.exists(xx):
                new_params.remove(xx)

new_params.sort()
for pp in params:
    if pp in new_params:
        src_param = nhm_params.get(pp)

        if src_param.name in ['nhm_seg', 'poi_gage_id', 'poi_type', 'poi_gage_segment']:
            # Skip segment-related parameters
            continue
            
        new_ps.parameters.add(src_param.name)

        ndims = src_param.ndims

        # if args.verbose:
        #     sys.stdout.write('\r                                       ')
        #     sys.stdout.write(f'\rProcessing {src_param.name} ')
        #     sys.stdout.flush()

        dim_order = [dd for dd in src_param.dimensions.keys()]
            
        if 'nsegment' in src_param.dimensions.keys():
            print(f'WARNING: {src_param.name} should not be included in single-HRU extraction')
        else:
            for dd in src_param.dimensions.keys():
                new_ps.parameters.get(src_param.name).dimensions.add(dd, new_ps.dimensions.get(dd).size)
                new_ps.parameters.get(src_param.name).datatype = src_param.datatype

        first_dimension = dim_order[0]
        outdata = None

        # Write out the data for the parameter
        if ndims == 1:
            # 1D Parameters
            if first_dimension == 'one':
                outdata = src_param.data
            elif first_dimension == 'nsegment':
                print('WARNING: Streamflow routing not available in single-HRU extraction')
            elif first_dimension == 'ndeplval':
                # This is really a 2D in disguise, however, it is stored in C-order unlike
                # other 2D arrays
                outdata = src_param.data.reshape((-1, 11))[tuple(uniq_deplcrv0), :].reshape((-1))
            elif first_dimension in HRU_DIMS:
                if pp == 'hru_deplcrv':
                    outdata = np.array(new_hru_deplcrv)
                elif pp == 'hru_segment':
                    # hru_segment is not needed for single HRU
                    pass
                else:
                    outdata = nhm_params.get_subset(pp, hru_order_subset)
            else:
                bandit_log.error(f'No rules to handle dimension {first_dimension}')
        elif ndims == 2:
            # 2D Parameters
            if first_dimension == 'nsegment':
                print('WARNING: nsegment-dimensioned parameters not allowed in single-HRU extraction')
                # outdata = nhm_params.get_subset(pp, new_nhm_seg)
            elif first_dimension in HRU_DIMS:
                outdata = nhm_params.get_subset(pp, hru_order_subset)
            else:
                err_txt = f'No rules to handle 2D parameter, {pp}, which contains dimension {first_dimension}'

        new_ps.parameters.get(src_param.name).data = outdata

# %%
# Write the new parameter file
header = [f'Written by Bandit_single_hru version TEST',
          f'ParamDb revision: TEST']

new_ps.write_parameter_file(out_param_file, header=header, prms_version=5)

# %%

# %%

# %%
