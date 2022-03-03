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

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import math
import numpy as np


# %%
def perv_comp(soil_moist, pptp, ptc, carea_max, smidx_coef, smidx_exp, contrib_frac, infil, srp):
    smidx = 0.0
    srpp = 0.0
    
    smidx = soil_moist + (0.5 * ptc)
    contrib_frac = min(carea_max, smidx_coef * 10**(smidx_exp * smidx))
    
    srpp = contrib_frac * pptp
    infil -= srpp
    srp += srpp
#     print(f'infil: {infil}')
#     print(f'srp: {srp}')
#     print(f'contrib_frac: {contrib_frac}')
    
    return (infil, srp, contrib_frac)


# %%
def check_capacity(soil_moist, soil_moist_max, snowinfil_max, infil, srp):
    capacity = soil_moist_max - soil_moist
    excess = infil - capacity
    
    if excess > snowinfil_max:
        srp += excess - snowinfil_max
        infil = snowinfil_max + capacity
        
#     print(f'infil: {infil}')
#     print(f'srp: {srp}')
    
    return (infil, srp)


# %%
infil = 0.0
avail_water = 0.0
NEARZERO = 1.1920929E-07

pptmix_nopack = False

imperv_stor = 0.0
imperv_stor_max = 0.05

snowinfil_max = 2.0000000
# snowmelt = 0.0264271
snowmelt = 0.0
# pkwater_equiv = 0.0687990
pkwater_equiv = 0.0

soil_moist = 0.6814164
soil_moist_max = 2.9048920

net_ppt = 0.1460000
net_rain = 0.0407064
net_snow = 0.1052936
# net_snow = 0.0

carea_max = 0.571828
smidx_coef = 0.020766
smidx_exp = 0.362845


contrib_frac = 0.0
srp = 0.0

print(net_rain + net_snow)
print(net_ppt - net_snow)

# %%
# pptmix_nopack
if pptmix_nopack:
    infil += net_rain
    avail_water += net_rain
    
    infil, srp, contrib_frac = perv_comp(soil_moist, net_rain, net_rain, carea_max, smidx_coef, smidx_exp, contrib_frac, infil, srp)
   
    print(f'infil: {infil}')
    print(f'srp: {srp}')
    print(f'contrib_frac: {contrib_frac}')
    print(f'avail_water: {avail_water}')
    
# infil = 0.0391870
# srp = 0.0015193
# contrib_frac = 0.0373238

# %%
# snowmelt > 0
    # pkwater_equiv, net_snow, soil_moist_max, snowinfil_max, imperv_stor, imperv_stor_max, avail_water
    # 0.0687990     0.1052936       2.9048920      2.0000000    0.0000000  0.0500000        0.0671335

    # soil_moist net_rain    net_ppt    snowmelt  contrib_frac infil     srp
    # 0.6814164  0.0407064  0.1460000  0.0264271  0.0373238  0.0656142  0.0015193
    
    
if snowmelt > 0.0:
    avail_water += snowmelt
    infil += snowmelt
    print('==== snowmelt > 0.0 ====')

    if pkwater_equiv > 0.0 or net_ppt - net_snow == 0.0:
        # In addition to snowmelt, there is a snowpack and all of the precipitation fell as snow
        print('    -- call check_capacity()')
        # Pervious area computation
        infil, srp = check_capacity(soil_moist, soil_moist_max, snowinfil_max, infil, srp)
    else:
        # Snowmelt occurred which depleted the snowpack
        print('    -- call perv_comp()')
        infil, srp, contrib_frac = perv_comp(soil_moist, snowmelt, net_ppt, carea_max, smidx_coef, smidx_exp, contrib_frac, infil, srp)
        
#     print(f'contrib_frac: {contrib_frac}')
elif pkwater_equiv < NEARZERO:
    # There is NO snowmelt and NO snowpack
    print('==== pkwater_equiv < NEARZERO')
    if net_snow < NEARZERO and net_rain > 0.0:
        # All or almost all the precipitation fell as rain (any snow fall was lost to sublimation)
        avail_water += net_rain
        infil += net_rain
        print('    ==== net_snow < NEARZERO and net_rain > 0.0')
        
        infil, srp, contrib_frac = perv_comp(soil_moist, net_rain, net_rain, carea_max, smidx_coef, smidx_exp, contrib_frac, infil, srp)    
elif infil > 0.0:
    # Some amount of snowpack exists; check if infil exceeds max snow infiltration rate.
    # The infiltation results from rain/snow mix on a snow-free surface.
    print('==== infil > 0.0 ====')
    infil, srp = check_capacity(soil_moist, soil_moist_max, snowinfil_max, infil, srp)
#     print(f'contrib_frac: {contrib_frac}')

print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')
print(f'avail_water: {avail_water}')

# %%
srp = 0.0
# infil = net_rain + snowmelt
infil = snowmelt

infil, srp, contrib_frac = perv_comp(soil_moist, net_rain, net_rain, carea_max, smidx_coef, smidx_exp, contrib_frac, infil, srp)

print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')


# %%
def contrib_fraction_smidx(carea_max, smidx_coef, smidx_exp, soil_moist, precip):
    # Compute the antecedent soil moisture
    smidx = soil_moist + (0.5 * precip)

    # Return the contribution fraction
    return min(carea_max, smidx_coef * 10.0**(smidx_exp * smidx))


# %%
def compute_infil_srp(contrib_frac, precip, infil, perv_runoff):
    adjusted_precip = contrib_frac * precip
    infil = infil - adjusted_precip
#     perv_runoff = perv_runoff + dble(adjusted_precip)    
    return (infil, perv_runoff + adjusted_precip)


# %%

# Example where pptmix_nopack is false
# We have snowmelt, but no snowpack and no precipitation (either net_ppt or net_snow)

#carea_max,  smidx_coef, smidx_exp, soil_moist, net_ppt, soil_moist_max, snowinfil_max
# 0.5718280  0.0207660  0.3628450  2.7847483  0.0000000   2.9048920  2.0000000
carea_max = 0.5718280
smidx_coef = 0.0207660
smidx_exp = 0.3628450
soil_moist = 2.7847483
soil_moist_max = 2.9048920

net_rain = 0.0
net_ppt = 0.0
net_snow = 0.0
upslope_hortonian = 0.0

pkwater_equiv = 0.0
snowinfil_max = 2.0
snowmelt = 0.0443047

infil = snowmelt
srp = 0.0

use_cascades = False
pptmix_nopack = False



# %%
cfrac = contrib_fraction_smidx(carea_max, smidx_coef, smidx_exp, soil_moist, net_ppt)
compute_infil_srp(cfrac, snowmelt, infil, srp)

# %%
compute_infil_srp(cfrac, snowmelt, infil, srp)

# %%
infil, srp = check_capacity(soil_moist, soil_moist_max, snowinfil_max, infil, srp)
print(f'infil: {infil}')
print(f'srp: {srp}')


# %%
def compute_infil_test(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, snowmelt, pkwater_equiv,
                       carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, infil, srp):
    if use_cascades:
        infil += max(0.0, upslope_hortonian)
        
    if pptmix_nopack:
        # Mixed event with no antecedent snowpack
        infil += net_rain
        
    if snowmelt > 0.0:
        infil += snowmelt
    elif pkwater_equiv == 0.0 and net_snow == 0.0 and net_rain > 0.0:
        infil += net_rain
    elif infil > 0.0:
        pass
    
    cfrac = contrib_fraction_smidx(carea_max, smidx_coef, smidx_exp, soil_moist, net_ppt)
    infil, srp = compute_infil_srp(cfrac, net_ppt, infil, srp)
    
    infil, srp = check_capacity(soil_moist, soil_moist_max, snowinfil_max, infil, srp)
    
    return (infil, srp, cfrac)


# %%
srp = 0.0
infil = 0.0

compute_infil_test(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, snowmelt, pkwater_equiv, 
                   carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, infil, srp)


# %%
def compute_infil(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, snowmelt, pkwater_equiv,
                  carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, infil, srp):
    
    avail_water = 0.0
    contrib_frac = 0.0
    
    if use_cascades:
        avail_water = max(0.0, upslope_hortonian)
        infil += max(0.0, upslope_hortonian)
        
        if infil > 0.0:
            infil, srp, contrib_frac = perv_comp(soil_moist, upslope_hortonian, upslope_hortonian, carea_max, 
                                                 smidx_coef, smidx_exp, contrib_frac, infil, srp)
        
    if pptmix_nopack:
        infil += net_rain
        avail_water += net_rain

        infil, srp, contrib_frac = perv_comp(soil_moist, net_rain, net_rain, carea_max, smidx_coef, smidx_exp, contrib_frac, infil, srp)

    if snowmelt > 0.0:
        avail_water += snowmelt
        infil += snowmelt
#         print('==== snowmelt > 0.0 ====')

        if pkwater_equiv > 0.0 or net_ppt - net_snow < NEARZERO:
#             print('    -- call check_capacity()')
            infil, srp = check_capacity(soil_moist, soil_moist_max, snowinfil_max, infil, srp)
        else:
#             print('    -- call perv_comp()')
            infil, srp, contrib_frac = perv_comp(soil_moist, snowmelt, net_ppt, carea_max, smidx_coef, smidx_exp, contrib_frac, infil, srp)
    elif pkwater_equiv < NEARZERO:
#         print('==== pkwater_equiv < NEARZERO')
        if net_snow < NEARZERO and net_rain > 0.0:
            avail_water += net_rain
            infil += net_rain
#             print('    ==== net_snow < NEARZERO and net_rain > 0.0')

            infil, srp, contrib_frac = perv_comp(soil_moist, net_rain, net_rain, carea_max, smidx_coef, smidx_exp, contrib_frac, infil, srp)    
    elif infil > 0.0:
#         print('==== infil > 0.0 ====')
        infil, srp = check_capacity(soil_moist, soil_moist_max, snowinfil_max, infil, srp)
    
    return (infil, srp, contrib_frac)


# %%
srp = 0.0
infil = 0.0

compute_infil(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, snowmelt, pkwater_equiv, 
                   carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, infil, srp)

# %%
# This example triggers for snowmelt AND (pkwater_equiv > 0.0_dp .or. net_ppt - net_snow < NEARZERO)
#          infil, pkwater_equiv, snowmelt
# infil:   0.0868653  0.5148234  0.0868653

# infil, net_rain, net_ppt, net_snow, carea_max, smidx_coef, smidx_exp, soil_moist, net_ppt, soil_moist_max, snowinfil_max, srp
#     melt_chk_infil:   0.0868653  0.0000000  0.0000000  0.0000000 0.5718280  0.0207660  0.3628450  0.0000000  0.0000000 ---   2.9048920  2.0000000  0.0000000

crap2, pkwater_equiv, snowmelt, crap3, net_rain, net_ppt, net_snow, carea_max, smidx_coef, smidx_exp, soil_moist, crap, soil_moist_max, snowinfil_max, srp = [0.0868653, 0.5148234, 0.0868653, 0.0868653, 
                                                                                                                       0.0000000, 0.0000000, 0.0000000, 0.5718280, 
                                                                                                                       0.0207660, 0.3628450, 0.0000000, 0.0000000, 
                                                                                                                       2.9048920, 2.0000000, 0.0000000]

# %%
srp = 0.0
infil = 0.0

compute_infil_test(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, snowmelt, pkwater_equiv, 
                   carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, infil, srp)

# %%
srp = 0.0
infil = 0.0

compute_infil(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, snowmelt, pkwater_equiv, 
                   carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, infil, srp)

# %%

# %%
# Triggers melt_comp
#   there is snowmelt AND no current snowpack AND (either a mixed-event or all-rain event)
#   Because there was an antecedent snowpack the rain component was added to pkwater_equiv which was
#   then completely depleted (it melted) during snowcomp. So ....
# 
# FINAL_infil:   0.2355339  0.0270727  0.1030922
# net_rain, net_ppt, net_snow, carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, snowinfil_max, pkwater_equiv, snowmelt, srp
vstr = '0.0670000  0.0670000  0.0000000  0.5718280  0.0207660  0.3628450  1.8843244  2.9048920  2.0000000  0.0000000  0.2626066  0.0270727'
vals = [float(xx) for xx in vstr.split()]

net_rain, net_ppt, net_snow, carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, snowinfil_max, pkwater_equiv, snowmelt, srp = vals

pptmix_nopack = False

srp = 0.0
infil = 0.0

infil, srp, contrib_frac = compute_infil_test(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, 
                                              snowmelt, pkwater_equiv, carea_max, smidx_coef, smidx_exp, soil_moist, 
                                              soil_moist_max, infil, srp)

print('---- new style ----')
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')

srp = 0.0
infil = 0.0

infil, srp, contrib_frac = compute_infil(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, 
                                         snowmelt, pkwater_equiv, carea_max, smidx_coef, smidx_exp, soil_moist, 
                                         soil_moist_max, infil, srp)

print('---- old style ----')
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')

# %%
infil = snowmelt
srp = 0.0

print(f'pkwater_equiv: {pkwater_equiv}')
print(f'snowmelt: {snowmelt}')
print(f'net_ppt: {net_ppt}')
print(f'net_rain: {net_rain}')
print(f'net_snow: {net_snow}')

cfrac = contrib_fraction_smidx(carea_max, smidx_coef, smidx_exp, soil_moist, snowmelt)
infil, srp = compute_infil_srp(cfrac, snowmelt, infil, srp)

print('-'*40)
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {cfrac}')
print(f'infil + srp: {(infil + srp)}')

infil = snowmelt
srp = 0.0

cfrac = contrib_fraction_smidx(carea_max, smidx_coef, smidx_exp, soil_moist, net_ppt)
infil, srp = compute_infil_srp(cfrac, net_ppt, infil, srp)

print('-'*40)
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {cfrac}')
print(f'infil + srp: {(infil + srp)}')

infil = snowmelt
srp = 0.0

cfrac = contrib_fraction_smidx(carea_max, smidx_coef, smidx_exp, soil_moist, net_ppt)
infil, srp = compute_infil_srp(cfrac, snowmelt, infil, srp)

print('-'*40)
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {cfrac}')
print(f'infil + srp: {(infil + srp)}')

# infil = snowmelt
# srp = 0.0
# cfrac = 0.0
# infil, srp = check_capacity(soil_moist, soil_moist_max, snowinfil_max, infil, srp)

# print('-'*40)
# print(f'infil: {infil}')
# print(f'srp: {srp}')
# print(f'contrib_frac: {cfrac}')

# %%
# Triggers melt_chk_infil

# infil:   0.2304425
#     melt_chk_infil:   0.2304425
vstr = '0.0160000  0.0160000  0.0000000  0.5718280  0.0207660  0.3628450  0.7463949  2.9048920  2.0000000  0.3767069  0.2304425  0.0000000'
#   FINAL_infil:   0.2304425  0.0000000  0.0000000

vals = [float(xx) for xx in vstr.split()]

net_rain, net_ppt, net_snow, carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, snowinfil_max, pkwater_equiv, snowmelt, srp = vals

srp = 0.0
infil = 0.0

infil, srp, contrib_frac = compute_infil_test(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, 
                                              snowmelt, pkwater_equiv, carea_max, smidx_coef, smidx_exp, soil_moist, 
                                              soil_moist_max, infil, srp)

print('---- new style ----')
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')

srp = 0.0
infil = 0.0

infil, srp, contrib_frac = compute_infil(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, 
                                         snowmelt, pkwater_equiv, carea_max, smidx_coef, smidx_exp, soil_moist, 
                                         soil_moist_max, infil, srp)

print('---- old style ----')
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')

# %%
# Triggers nopack_infil and melt_chk_infil
# infil:   0.0407064
#     nopack_infil:   0.0391870
# infil:   0.0656142
#     melt_chk_infil:   0.0656142
vstr = '0.0407064  0.1460000  0.1052936  0.5718280  0.0207660  0.3628450  0.6814164  2.9048920  2.0000000  0.0687990  0.0264271  0.0015193'
#   FINAL_infil:   0.0656142  0.0015193  0.0373238

vals = [float(xx) for xx in vstr.split()]

net_rain, net_ppt, net_snow, carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, snowinfil_max, pkwater_equiv, snowmelt, srp = vals


pptmix_nopack = True

srp = 0.0
infil = 0.0

infil, srp, contrib_frac = compute_infil_test(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, 
                                              snowmelt, pkwater_equiv, carea_max, smidx_coef, smidx_exp, soil_moist, 
                                              soil_moist_max, infil, srp)

print('---- new style ----')
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')

srp = 0.0
infil = 0.0

infil, srp, contrib_frac = compute_infil(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, 
                                         snowmelt, pkwater_equiv, carea_max, smidx_coef, smidx_exp, soil_moist, 
                                         soil_moist_max, infil, srp)

print('---- old style ----')
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')

# %%
# Triggers pkw_infil
# infil:   0.1119039
#     pkw_infil:   0.1059145
vstr = '0.1119039  0.1119039  0.0000000  0.5718280  0.0207660  0.3628450  1.0772618  2.9048920  2.0000000  0.0000000  0.0000000  0.0059893'
#   FINAL_infil:   0.1059145  0.0059893  0.0535222

vals = [float(xx) for xx in vstr.split()]

net_rain, net_ppt, net_snow, carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, snowinfil_max, pkwater_equiv, snowmelt, srp = vals

pptmix_nopack = False

srp = 0.0
infil = 0.0

infil, srp, contrib_frac = compute_infil_test(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, 
                                              snowmelt, pkwater_equiv, carea_max, smidx_coef, smidx_exp, soil_moist, 
                                              soil_moist_max, infil, srp)

print('---- new style ----')
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')

srp = 0.0
infil = 0.0

infil, srp, contrib_frac = compute_infil(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, 
                                         snowmelt, pkwater_equiv, carea_max, smidx_coef, smidx_exp, soil_moist, 
                                         soil_moist_max, infil, srp)

print('---- old style ----')
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')

# %%
# Triggers nopack_infil and infil_chk_infil
# infil:   0.0198876
#     nopack_infil:   0.0181654
# infil:   0.0181654
#     infil_chk_infil:   0.0181654
vstr = '0.0198876  0.1230000  0.1031124  0.5718280  0.0207660  0.3628450  1.6991829  2.9048920  2.0000000  0.0896566  0.0000000  0.0017222'
#   FINAL_infil:   0.0181654  0.0017222  0.0865966

vals = [float(xx) for xx in vstr.split()]

net_rain, net_ppt, net_snow, carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, snowinfil_max, pkwater_equiv, snowmelt, srp = vals

pptmix_nopack = True

srp = 0.0
infil = 0.0

infil, srp, contrib_frac = compute_infil_test(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, 
                                              snowmelt, pkwater_equiv, carea_max, smidx_coef, smidx_exp, soil_moist, 
                                              soil_moist_max, infil, srp)

print('---- new style ----')
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')

srp = 0.0
infil = 0.0

infil, srp, contrib_frac = compute_infil(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, 
                                         snowmelt, pkwater_equiv, carea_max, smidx_coef, smidx_exp, soil_moist, 
                                         soil_moist_max, infil, srp)

print('---- old style ----')
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')

# %%
# Triggers nopack_infil and melt_comp_infil
# infil:   0.4799967
#     nopack_infil:   0.3844725
# infil:   0.4814759
#     melt_comp_infil:   0.4613729
vstr = '0.4799967  0.5770000  0.0970033  0.5718280  0.0207660  0.3628450  2.4650743  2.9048920  2.0000000  0.0000000  0.0970033  0.1156271'
#   FINAL_infil:   0.4613729  0.1156271  0.2072400

vals = [float(xx) for xx in vstr.split()]

net_rain, net_ppt, net_snow, carea_max, smidx_coef, smidx_exp, soil_moist, soil_moist_max, snowinfil_max, pkwater_equiv, snowmelt, srp = vals

pptmix_nopack = True

srp = 0.0
infil = 0.0

infil, srp, contrib_frac = compute_infil_test(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, 
                                              snowmelt, pkwater_equiv, carea_max, smidx_coef, smidx_exp, soil_moist, 
                                              soil_moist_max, infil, srp)

print('---- new style ----')
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')

srp = 0.0
infil = 0.0

infil, srp, contrib_frac = compute_infil(use_cascades, pptmix_nopack, upslope_hortonian, net_ppt, net_rain, 
                                         snowmelt, pkwater_equiv, carea_max, smidx_coef, smidx_exp, soil_moist, 
                                         soil_moist_max, infil, srp)

print('---- old style ----')
print(f'infil: {infil}')
print(f'srp: {srp}')
print(f'contrib_frac: {contrib_frac}')

# %%
