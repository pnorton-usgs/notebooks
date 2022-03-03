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
# Parameters
dprst_depth_avg = 57.9203
dprst_frac = 0.051539
dprst_frac_init = 0.093112
dprst_frac_open = 1.0
hru_area = 149514.25399
hru_percent_imperv = 0.000018
op_flow_thres = 0.9264
va_open_exp = 0.001
va_clos_exp = 0.001

soil_moist_init = 1.0
soil_moist_init_frac = 0.370547
soil_moist_max = 2.904892
soil_rechr_init = 0.5
soil_rechr_init_frac = 0.371024
soil_rechr_max_frac = 0.994335

# The variables...
imperv_frac_flag = 1
dprst_frac_flag = 1
dprst_depth_flag = 1

check_imperv = False
check_dprst_frac = False
dprst_flag = 1
dprst_clos_flag = True
dprst_open_flag = True

soil_moist = soil_moist_init_frac * soil_moist_max
soil_rechr = soil_rechr_max_frac * soil_moist_max

soil_moist_tmp = 0.0
soil_rechr_tmp = 0.0

dprst_area_max = 7705.81513639
dprst_vol_clos = 0.0
dprst_vol_open = 41558.0387633

hru_area_imperv = 2.69125657182
hru_area_perv = 141805.747597
hru_frac_perv = 0.948443
imperv_stor = 0.1

hru_area_perv_NEW = 0.0
hru_percent_imperv_new = 0.000017
dprst_frac_NEW = 0.06
dprst_area_max_tmp = 0.0


# %% [markdown]
# ### Deal with hru_percent_imperv first

# %%

hru_area_imperv = hru_percent_imperv * hru_area

# NOTE: dprst_area_max is not handled yet
hru_area_perv = hru_area - hru_area_imperv

hru_frac_perv = hru_area_perv / hru_area

# NOTE: need to deal with lakes 

# %%
# Some scratch junk
hru_percent_imperv_OLD = 0.000018
hru_percent_imperv_NEW = 0.0
hru_percent_imperv_CHG = 0.0

if hru_percent_imperv_NEW > 0.0:
    hru_percent_imperv_CHG = hru_percent_imperv_OLD / hru_percent_imperv_NEW

hru_percent_imperv_EXIST_CHG = hru_percent_imperv_OLD / hru_frac_perv

if hru_percent_imperv_NEW > 0.999:
    print('ERROR')
elif imperv_stor > 0.0:
    if hru_percent_imperv_NEW > 0.0:
        imperv_stor *= hru_percent_imperv_CHG
    else:
        print('WARNING')
        tmp_stor = imperv_stor * hru_percent_imperv_EXIST_CHG
        print('tmp_stor = {}'.format(tmp_stor))
        
print('imperv_stor = {}'.format(imperv_stor))
print('_OLD = {}'.format(hru_percent_imperv_OLD))
print('_NEW = {}'.format(hru_percent_imperv_NEW))
print('_CHG = {}'.format(hru_percent_imperv_CHG))
print('_EXI = {}'.format(hru_percent_imperv_EXIST_CHG))

# %% [markdown]
# ### Handle changes to imperv_stor_max

# %%
# Read in imperv_stor_max updates for current timestep


# %% [markdown]
# ### Handle changes to dprst_depth_avg and dprst_frac

# %%
dprst_area_max = dprst_frac * hru_area
hru_area_perv -= dprst_area_max
hru_frac_perv = hru_area_perv / hru_area

# where active and dprst_area_max > 0
dprst_area_open_max = dprst_area_max * dprst_frac_open
dprst_area_clos_max = dprst_area_max - dprst_area_open_max
dprst_frac_clos = 1.0 - dprst_frac_open

if hru_percent_imperv + dprst_frac > 0.999:
    print('ERROR: impervious plus depression fraction > 0.999')
    
# Dynread is wrong
has_closed_dprst = dprst_area_clos_max > 0.0
has_open_dprst = dprst_area_open_max > 0.0

# For all the active HRUs...
if has_closed_dprst and dprst_area_max > 0.0:
    dprst_vol_clos_max = dprst_area_clos_max * dprst_depth_avg
if has_open_dprst and dprst_area_max > 0.0:
    dprst_vol_open_max = dprst_area_open_max * dprst_depth_avg
    
if has_closed_dprst and dprst_area_max > 0.0:
    dprst_vol_clos = dprst_frac_init * dprst_vol_clos_max
    
if has_open_dprst and dprst_area_max > 0.0:
    dprst_vol_open = dprst_frac_init * dprst_vol_open_max
    
if dprst_area_max > 0.0:
    dprst_vol_thres_open = op_flow_thres * dprst_vol_open_max
    
if dprst_vol_open > 0.0 and dprst_area_max > 0.0:
    pass
    # dprst_area_open = depression_surface_area(dprst_vol_open, dprst_vol_open_max, dprst_area_open_max, va_open_exp)

if dprst_vol_clos > 0.0 and dprst_area_max > 0.0:
    pass
    # dprst_area_clos = depression_surface_area(dprst_vol_clos, dprst_vol_clos_max, dprst_area_clos_max, va_clos_exp)


if dprst_vol_clos_max > 0.0:
    dprst_vol_clos_frac = dprst_vol_clos / dprst_vol_clos_max
if dprst_vol_open_max > 0.0:
    dprst_vol_open_frac = dprst_vol_open / dprst_vol_open_max
    
if dprst_area_max > 0.0:
    dprst_vol_frac = (dprst_vol_open + dprst_vol_clos) / (dprst_vol_open_max + dprst_vol_clos_max)
    dprst_stor_hru = (dprst_vol_open + dprst_vol_clos) / hru_area
    

# %%

# %%

# %%

# %%
soil_moist_tmp = soil_moist
soil_rechr_tmp = soil_rechr
print('------- BEFORE --------')
print('dprst_area_max = {}'.format(dprst_area_max))
print('soil_moist = {}'.format(soil_moist))
print('soil_rechr = {}'.format(soil_rechr))
print('hru_area_imperv = {}'.format(hru_area_imperv))
print('hru_area_perv = {}'.format(hru_area_perv))
print('hru_frac_perv = {}'.format(hru_frac_perv))
print('hru_percent_imperv = {}'.format(hru_percent_imperv))
print('imperv_stor = {}'. format(imperv_stor))


dprst_frac_flag = 1
check_imperv = True

if dprst_frac_flag == 1:
    # ??? Why is this done? We already have dprst_area_max from the prior timestep.
    dprst_area_max = dprst_frac * hru_area
    print('UPDATE: dprst_area_max = {}'.format(dprst_area_max))
    
if check_imperv:
    # frac_imperv could be a new or updated value
    frac_imperv = hru_percent_imperv_new
    print('frac_imperv = {}'.format(frac_imperv))
    
    if frac_imperv > 0.999:
        print('ERROR: Dynamic value of the hru_percent_imperv > 0.999')
    elif imperv_stor > 0.0:
        if frac_imperv > 0.0:
            imperv_stor *= hru_percent_imperv / frac_imperv
        else:
            print('WARNING: Dynamic impervious changed to 0 when impervious storage > 0')
            tmp = imperv_stor * hru_percent_imperv / hru_frac_perv
            soil_moist_tmp += tmp
            soil_rechr_tmp += tmp
            imperv_stor = 0.0
            
    hru_percent_imperv = frac_imperv
    hru_area_imperv = hru_area * frac_imperv

print('------- AFTER --------')
print('dprst_area_max = {}'.format(dprst_area_max))
print('soil_moist = {}'.format(soil_moist))
print('soil_rechr = {}'.format(soil_rechr))
print('soil_moist_tmp = {}'.format(soil_moist_tmp))
print('soil_rechr_tmp = {}'.format(soil_rechr_tmp))
print('hru_area_imperv = {}'.format(hru_area_imperv))
print('hru_area_perv = {}'.format(hru_area_perv))
print('hru_frac_perv = {}'.format(hru_frac_perv))
print('hru_percent_imperv = {}'.format(hru_percent_imperv))
print('imperv_stor = {}'. format(imperv_stor))

# %%
check_dprst_frac = False

if dprst_frac_flag == 1:
    dprst_area_max_tmp = dprst_area_max
    
    if check_dprst_frac:
        dprst_area_max_tmp = dprst_frac_NEW * hru_area
        
    if dprst_area_max_tmp > 0.0:
        if hru_percent_imperv + dprst_area_max_tmp > 0.999:
            if 0.999 - hru_percent_imperv < 0.001:
                print('ERROR: Fraction impervious + fraction dprst > 0.999 for HRU')
                
    tmp = dprst_vol_open + dprst_vol_clos

    if dprst_area_max_tmp == 0.0 and tmp > 0.0:
        print('WARNING: dprst_area reduced to 0 with storage > 0')
        new_storage = tmp / dprst_area_max_tmp / hru_frac_perv
        soil_moist_tmp += new_storage
        soil_rechr_tmp += new_storage
        dprst_vol_open = 0.0
        dprst_vol_clos = 0.0
        print('UPDATE: soil_moist_tmp = {}'.format(soil_moist_tmp))
        print('UPDATE: soil_rechr_tmp = {}'.format(soil_rechr_tmp))
        print('UPDATE: dprst_vol_clos = {}'.format(dprst_vol_clos))
        print('UPDATE: dprst_vol_open = {}'.format(dprst_vol_open))
    
    print('   OLD: dprst_area_max = {}'. format(dprst_area_max))
    dprst_area_open_max = dprst_area_max_tmp * dprst_frac_open
    dprst_area_clos_max = dprst_area_max_tmp - dprst_area_open_max
    dprst_area_max = dprst_area_max_tmp
    print('UPDATE: dprst_area_open_max = {}'.format(dprst_area_open_max))
    print('UPDATE: dprst_area_clos_max = {}'.format(dprst_area_clos_max))
    print('UPDATE: dprst_area_max = {}'. format(dprst_area_max))
    
    if dprst_area_clos_max > 0.0:
        dprst_clos_flag = False
    if dprst_area_open_max > 0.0:
        dprst_open_flag = False
     
    print('   OLD: dprst_frac = {}'.format(dprst_frac))
    dprst_frac = dprst_area_max / hru_area
    dprst_vol_clos_max = dprst_area_clos_max * dprst_depth_avg
    dprst_vol_open_max = dprst_area_open_max * dprst_depth_avg
    dprst_vol_thres_open = dprst_vol_open_max * op_flow_thres
    print('UPDATE: dprst_frac = {}'.format(dprst_frac))
    print('UPDATE: dprst_vol_clos_max = {}'.format(dprst_vol_clos_max))
    print('UPDATE: dprst_vol_open_max = {}'.format(dprst_vol_open_max))
    print('UPDATE: dprst_vol_thres_open = {}'.format(dprst_vol_thres_open))
    
if check_imperv or dprst_frac_flag == 1:
    hru_area_perv_NEW = hru_area - hru_area_imperv
    print('UPDATE: hru_area_perv_NEW = {}'.format(hru_area_perv_NEW))
    
    if dprst_flag == 1:
        hru_area_perv_NEW -= dprst_area_max
        print('UPDATE: hru_area_perv_NEW = {}'.format(hru_area_perv_NEW))
        
        if dprst_area_max + hru_area_imperv > 0.999 * hru_area:
            print('ERROR: Impervious + depression area > 0.99 * hru_area')
    
    if hru_area_perv != hru_area_perv_NEW:
        print('========= hru_area_perv != hru_area_perv_NEW =============')
        if hru_area_perv_NEW < 0.0000001:
            print('ERROR: Pervious area error for dynamic parameter')
        
        tmp = hru_area_perv / hru_area_perv_NEW
        print('  TEMP: tmp = {}'.format(tmp))
        soil_moist_tmp *= tmp
        soil_rechr_tmp *= tmp
        print('UPDATE: soil_moist_tmp = {}'.format(soil_moist_tmp))
        print('UPDATE: soil_rechr_tmp = {}'.format(soil_rechr_tmp))
        
        print('   OLD: hru_area_perv = {}'.format(hru_area_perv))
        print('   OLD: hru_frac_perv = {}'.format(hru_frac_perv))
        hru_area_perv = hru_area_perv_NEW
        hru_frac_perv = hru_area_perv / hru_area
        print('UPDATE: hru_area_perv = {}'.format(hru_area_perv))
        print('UPDATE: hru_frac_perv = {}'.format(hru_frac_perv))

# %%
if dprst_depth_flag == 1:
    # Update dprst_depth_avg with any new values for timestep
    if dprst_area_max > 0.0:
        dprst_vol_clos_max = dprst_area_clos_max * dprst_depth_avg
        dprst_vol_open_max = dprst_area_open_max * dprst_depth_avg
        
        if dprst_vol_open_max > 0.0:
            dprst_vol_open_frac = dprst_vol_open / dprst_vol_open_max
            
        if dprst_vol_clos_max > 0.0:
            dprst_vol_clos_frac =dprst_vol_clos / dprst_vol_clos_max
            
        dprst_vol_frac = (dprst_vol_open + dprst_vol_clos) / (dprst_vol_open_max + dprst_vol_clos_max)
        

# %%
soil_moist = 566.0
soil_adj1 = soil_moist
soil_adj2 = 1.0
imperv_stor = 0.1
hru_area = 149514.25399
hru_percent_imperv = 0.000018
hru_area_imperv = hru_percent_imperv * hru_area
hru_area_perv = hru_area - hru_area_imperv
hru_frac_perv = 1.0 - hru_percent_imperv

hru_percent_imperv_OLD = 0.000017
hru_area_imperv_OLD = hru_percent_imperv_OLD * hru_area
hru_area_perv_OLD = hru_area - hru_area_imperv_OLD
hru_frac_perv_OLD = 1.0 - hru_percent_imperv_OLD

adj = imperv_stor * hru_percent_imperv_OLD / hru_frac_perv_OLD
print(adj)
print('soil_adj1 = {}'.format(soil_adj1))
print('soil_adj2 = {}'.format(soil_adj2))
soil_adj1 += adj
soil_adj2 += adj
print('UPDATE: soil_adj1 = {}'.format(soil_adj1))
print('UPDATE: soil_adj2 = {}'.format(soil_adj2))
#               adj = this%imperv_stor(chru) * hru_percent_imperv_old(chru) / hru_frac_perv_old(chru)
#               soil_moist_chg(chru) = soil_moist_chg(chru) + adj
#               soil_rechr_chg(chru) = soil_rechr_chg(chru) + adj

# %%


adj2 = hru_area_perv_OLD / hru_area_perv
print(adj2)
soil_adj1 *= adj2
soil_adj2 *= adj2
print('UPDATE: soil_adj1 = {}'.format(soil_adj1))
print('UPDATE: soil_adj2 = {}'.format(soil_adj2))
#             adj = hru_area_perv_old(chru) / this%hru_area_perv(chru)
#             soil_moist_chg(chru) = soil_moist_chg(chru) * adj
#             soil_rechr_chg(chru) = soil_rechr_chg(chru) * adj

# %%
dprst_depth_avg = 57.9203
dprst_frac = 0.051539
# dprst_frac = 0.0
dprst_frac_init = 0.093112
dprst_frac_open = 0.9

dprst_area_max = dprst_frac * hru_area

dprst_area_open_max = dprst_area_max * dprst_frac_open
dprst_area_clos_max = dprst_area_max - dprst_area_open_max

dprst_vol_clos_max = dprst_area_clos_max * dprst_depth_avg
dprst_vol_open_max = dprst_area_open_max * dprst_depth_avg

dprst_vol_clos = dprst_frac_init * dprst_vol_clos_max
dprst_vol_open = dprst_frac_init * dprst_vol_open_max

tmp = dprst_vol_open + dprst_vol_clos
print(hru_frac_perv_OLD)

adj3 = tmp / hru_frac_perv_OLD
print(adj3)
soil_adj1 += adj3
soil_adj2 += adj3
print('UPDATE: soil_adj1 = {}'.format(soil_adj1))
print('UPDATE: soil_adj2 = {}'.format(soil_adj2))
#             tmp = this%dprst_vol_open(chru) + this%dprst_vol_clos(chru)

#             if (tmp > 0.0) then
#               write(output_unit, *) 'WARNING: dprst_area_max reduced to 0 with storage > 0'
#               write(output_unit, *) '         Storage was added to soil_moist and soil_rechr'

#               adj = tmp / hru_frac_perv_old(chru)
#               soil_moist_chg(chru) = soil_moist_chg(chru) + adj
#               soil_rechr_chg(chru) = soil_rechr_chg(chru) + adj

# %%
soil_moist + soil_adj2

# %%
