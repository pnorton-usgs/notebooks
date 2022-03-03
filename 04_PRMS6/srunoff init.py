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
#     display_name: Python [conda env:idp_bandit]
#     language: python
#     name: conda-env-idp_bandit-py
# ---

# %%
import math

# %%
# Parameters
dprst_depth_avg = 57.9203
dprst_frac = 0.051539
dprst_frac_init = 0.093112
dprst_frac_open = 0.9
hru_area = 149514.25399
hru_percent_imperv = 0.000018
op_flow_thres = 0.9264
va_open_exp = 0.001
va_clos_exp = 0.001

# other variables
dprst_area_max = 0.0
dprst_area_clos = 0.0
dprst_area_clos_max = 0.0
dprst_area_open = 0.0
dprst_area_open_max = 0.0
dprst_stor_hru = 0.0
dprst_vol_thres_open = 0.0
dprst_vol_clos = 0.0
dprst_vol_clos_frac = 0.0
dprst_vol_clos_max = 0.0
dprst_vol_open = 0.0
dprst_vol_open_frac = 0.0
dprst_vol_open_max = 0.0
imperv_stor = 0.0

has_closed_dprst = False
has_open_dprst = False

# Init stuff
dprst_area_max = dprst_frac * hru_area

hru_area_imperv = hru_percent_imperv * hru_area
hru_area_perv = hru_area - hru_area_imperv - dprst_area_max
hru_frac_perv = hru_area_perv / hru_area

# dprst_frac_clos = 1.0 - dprst_frac_open
# dprst_area_open_max = dprst_area_max * dprst_frac_open
# dprst_area_clos_max = dprst_area_max - dprst_area_open_max


print('hru_area = {}'.format(hru_area))
print('hru_area_imperv = {}'.format(hru_area_imperv))
print('hru_area_perv = {}'.format(hru_area_perv))
print('hru_frac_perv = {}'.format(hru_frac_perv))
print('dprst_area_max = {}'.format(dprst_area_max))

if dprst_area_max > 0.0:
    dprst_frac_clos = 1.0 - dprst_frac_open
    dprst_area_open_max = dprst_area_max * dprst_frac_open
    dprst_area_clos_max = dprst_area_max - dprst_area_open_max
    
    if hru_percent_imperv + dprst_frac > 0.999:
        print('ERROR: Impervious plus depression fraction > 0.999 for HRU')

print('dprst_area_clos_max = {}'.format(dprst_area_clos_max))
print('dprst_area_open_max = {}'.format(dprst_area_open_max))
print('dprst_frac_clos = {}'.format(dprst_frac_clos))


has_closed_dprst = dprst_area_clos_max > 0.0
has_open_dprst = dprst_area_open_max > 0.0

print('-'*40)
print('has_closed_dprst = {}'.format(has_closed_dprst))
print('has_open_dprst = {}'.format(has_open_dprst))
print('-'*40)

if dprst_area_max > 0.0:
    if has_closed_dprst:
        dprst_vol_clos_max = dprst_area_clos_max * dprst_depth_avg
    if has_open_dprst:
        dprst_vol_open_max = dprst_area_open_max * dprst_depth_avg

print('dprst_vol_clos_max = {}'.format(dprst_vol_clos_max))
print('dprst_vol_open_max = {}'.format(dprst_vol_open_max))

# init_vars_from_file stuff
if has_closed_dprst:
    dprst_vol_clos = dprst_frac_init * dprst_vol_clos_max
    
if has_open_dprst:
    dprst_vol_open = dprst_frac_init * dprst_vol_open_max

    
dprst_vol_thres_open = op_flow_thres * dprst_vol_open_max

print('dprst_vol_open = {}'.format(dprst_vol_open))
print('dprst_vol_clos = {}'.format(dprst_vol_clos))
print('dprst_vol_thres_open = {}'.format(dprst_vol_thres_open))

if dprst_vol_open > 0.0:
    open_vol_r = dprst_vol_open / dprst_vol_open_max
    
    if open_vol_r == 0.0:
        frac_op_ar = 0.0
    elif open_vol_r > 1.0:
        frac_op_ar = 1.0
    else:
        frac_op_ar = math.exp(va_open_exp * math.log(open_vol_r))
        
    dprst_area_open = dprst_area_open_max * frac_op_ar
    
    if dprst_area_open > dprst_area_open_max:
        dprst_area_open = dprst_area_open_max

if dprst_vol_clos > 0.0:
    clos_vol_r = dprst_vol_clos / dprst_vol_clos_max
    
    if clos_vol_r == 0:
        frac_cl_ar = 0.0
    elif clos_vol_r > 1.0:
        frac_cl_ar = 1.0
    else:
        frac_cl_ar = math.exp(va_clos_exp * math.log(clos_vol_r))
        
    dprst_area_clos = dprst_area_clos_max * frac_cl_ar
    
    if dprst_area_clos > dprst_area_clos_max:
        dprst_area_clos = dprst_area_clos_max
        
dprst_stor_hru = (dprst_vol_open + dprst_vol_clos) / hru_area

print('dprst_area_clos = {}'.format(dprst_area_clos))
print('dprst_area_open = {}'.format(dprst_area_open))
print('dprst_stor_hru = {}'.format(dprst_stor_hru))

if dprst_vol_open_max > 0.0:
    dprst_vol_open_frac = dprst_vol_open / (dprst_vol_open_max + dprst_vol_clos_max)
#     dprst_vol_open_frac = dprst_vol_open / dprst_vol_open_max
    
if dprst_vol_clos_max > 0.0:
    dprst_vol_clos_frac = dprst_vol_clos / (dprst_vol_open_max + dprst_vol_clos_max)
#     dprst_vol_clos_frac = dprst_vol_clos / dprst_vol_clos_max
    
dprst_vol_frac = (dprst_vol_open + dprst_vol_clos) / (dprst_vol_open_max + dprst_vol_clos_max)

print('dprst_vol_open_frac = {}'.format(dprst_vol_open_frac))
print('dprst_vol_clos_frac = {}'.format(dprst_vol_clos_frac))
print('dprst_vol_frac = {}'.format(dprst_vol_frac))

# %%
dprst_vol_open_frac + dprst_vol_clos_frac

# %%
