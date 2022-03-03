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
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import numpy as np
import pandas as pd

from collections import OrderedDict

from pyPRMS.ParameterFile import ParameterFile

# %%
step = 'SOM'    # one of AET, RUN, RCH, SCA, SOM, ALL
hru_id = '66119'
workdir = f'/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11/byHRU_sample/HRU{hru_id}'

calib_file = f'{workdir}/parameter_info'
param_file = f'{workdir}/params_byHRU'
ofs_file = f'{workdir}/RESULTS/{step}_OFS_HRU{hru_id}'
mean_param_file = f'{workdir}/RESULTS/{step}_PARAMS_HRU{hru_id}'

cal_sca = False

# %%
# 1 RUN    0.34418
# 2 AET    0.27683
# 3 SCA    0.16787
# 4 SOM    0.13448
# 5 RCH    0.07663
# 6 ALL    1.00000

# %%
calib_df = pd.read_csv(calib_file, sep='\s+', index_col=['Parameter'])
calib_df.head()

# %%
# lines: 1103-1212
opt_params = 0
param_names = OrderedDict()

for row in calib_df.index:
    if not cal_sca and calib_df['Snow?'][row] == 1:
        continue

    opt_params += 1
    param_names[row] = row
    
    if calib_df['byMONTH'][row] == 1:
        for pp in range(3):
            param_names[f'{row}.{pp+1}'] = row
            
        opt_params += 3
        
print(f'Number of optimization parameters: {opt_params}')

# %%
# Read parameter file to get initial mean parameter values
params = ParameterFile(param_file, verbose=True) 

init_mean = OrderedDict()

for op, pp in param_names.items():
    init_mean[pp] = params.parameters.get(pp).data.mean()
#     print(op, params.parameters.get(pp).data.mean())

# %%
# Print of the initial mean parameter values
out = []

for op, pp in param_names.items():
    out.append(init_mean[pp])
    
print(out)

# %%

# %%
p = 'carea_max'

smin = calib_df['RANGEmin'][p]
smax = calib_df['RANGEmax'][p]
pcnt = calib_df['%'][p]
init_val = init_mean[p]

bound = smax - smin
lbnd = max(init_val - (pcnt * bound), smin)
ubnd = min(lbnd + (2.0 * pcnt * bound), smax)

print(f'smin: {smin}  smax: {smax}')
print(f'pcnt: {pcnt}  init_val: {init_val}')
print(f'bound: {bound}')
print(f'lbnd: {lbnd}')
print(f'ubnd: {ubnd}')

# %%

# %%

# %%

# %%

# %%

# %%
# line: 528

# lines: 571-603
#       lower_bnd = rng_min
#       upper_bnd = rng_max
#       bound = rng_max - rng_min

#       where (pcnt_or_rng == 'R')
#         lower_bnd = max(lower_bnd, rng_min)
#         upper_bnd = min(upper_bnd, rng_max)
#       elsewhere
#         lower_bnd = max(A - (pcnt * bound), rng_min)
#         upper_bnd = min(lower_bnd + (2.0 * pcnt * bound), rng_max)
#       end where
orig_diff = OrderedDict()
lower_bnd = OrderedDict()
upper_bnd = OrderedDict()

for op, param in param_names.items():
    orig_diff[op] = 0.01 * (calib_df['RANGEmax'][param] - calib_df['RANGEmin'][param])
    
    if calib_df['%orRANGE'][param] == 'R':
        lower_bnd[op] = calib_df['RANGEmin'][param]
        upper_bnd[op] = calib_df['RANGEmax'][param]
    else:
        bound = calib_df['RANGEmax'][param] - calib_df['RANGEmin'][param]
        lower_bnd[op] = max(init_mean[param] - (calib_df['%'][param] * bound), calib_df['RANGEmin'][param])
        upper_bnd[op] = min(lower_bnd[op] + (2.0 * calib_df['%'][param] * bound), calib_df['RANGEmax'][param])
        
print(orig_diff)

# %%
lower_out = []
upper_out = []

for op, param in param_names.items():
    lower_out.append(lower_bnd[op])
    upper_out.append(upper_bnd[op])
    
print(lower_out)
print(upper_out)

# %%
# Make sure init_mean is within the bounds
# Only done once per STEP
out = []

for op, param in param_names.items():
    init_mean[param] = max(init_mean[param], lower_bnd[op])
    init_mean[param] = min(init_mean[param], upper_bnd[op])
    
    out.append(init_mean[param])
    
print(out)

# %%
ofs_df = pd.read_csv(ofs_file, sep='\s+')
ofs_df.drop(columns=['HRU'], inplace=True)
ofs_df.set_index('RUN', inplace=True)
ofs_df.head()

# %%
param_df = pd.read_csv(mean_param_file, sep='\s+')
param_df.drop(columns=['HRU'], inplace=True)
param_df.set_index('RUN', inplace=True)
param_df.head()

# %%
res_df = pd.merge(ofs_df, param_df, how='left', left_index=True, right_index=True)
res_df.head()

# %%
res_df.sort_values(by=['prmsOF'], axis=0, inplace=True)

# %%
#       n = int(0.25 * real(num_ofs))
#       xmax = OFsort(n)
top25 = int(0.25 * len(res_df))
max_of = res_df.iloc[top25, 0]

print(max_of)

# %%
# line: 518-525
res_df.drop(columns=['prmsOF', 'ofRUN', 'ofAET', 'ofSCA', 'ofRCH', 'ofSOM'], inplace=True)
new_lower_bnd = res_df.iloc[0:top25].min()
new_upper_bnd = res_df.iloc[0:top25].max()

# %%
new_lower_bnd

# %%
diff = new_upper_bnd - new_lower_bnd
print(diff)

# %%
# lines: 532-554

# adjust_params()
# if (diff(c_par) < RNGE1pcnt(c_par)) then
#     new_upper_bnd(c_par) = pcal(num_ofs, c_par) + (0.5 * RNGE1pcnt(c_par))
#     new_lower_bnd(c_par) = pcal(num_ofs, c_par) - (0.5 * RNGE1pcnt(c_par))
#     new_upper_bnd(c_par) = min(new_upper_bnd(c_par), upper_bnd(c_par))
#     new_lower_bnd(c_par) = max(new_lower_bnd(c_par), lower_bnd(c_par))

# lower_bnd(c_par) = new_lower_bnd(c_par)
# upper_bnd(c_par) = new_upper_bnd(c_par)

for xx, yy in orig_diff.items():
    if diff[xx] < yy:
        print('\t **** new less than orig')
        new_upper_bnd[xx] = min(param_df.iloc[-1][xx] + (0.5 * yy), upper_bnd[xx])
        new_lower_bnd[xx] = max(param_df.iloc[-1][xx] - (0.5 * yy), lower_bnd[xx])

    lower_bnd[xx] = new_lower_bnd[xx]
    upper_bnd[xx] = new_upper_bnd[xx]
    print(xx, diff[xx], yy)

# %%
lower_out = []
upper_out = []

for op, param in param_names.items():
    lower_out.append(lower_bnd[op])
    upper_out.append(upper_bnd[op])
    
print(lower_out)
print(upper_out)

# %%
print(param_df.iloc[-1].tolist())

# %%

# %%
p = 'carea_max'

smin = calib_df['RANGEmin'][p]
smax = calib_df['RANGEmax'][p]
pcnt = calib_df['%'][p]
init_val = param_df.iloc[-1].tolist()[0]

bound = smax - smin
lbnd = lower_bnd[p]
ubnd = upper_bnd[p]

print(f'smin: {smin}  smax: {smax}')
print(f'pcnt: {pcnt}  init_val: {init_val}')
print(f'bound: {bound}')
print(f'lbnd: {lbnd}')
print(f'ubnd: {ubnd}')

# %%
