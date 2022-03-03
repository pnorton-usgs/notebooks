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

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
# pkwater_equiv       freeh2o  snow_evap  pk_ice       pss                  potet_sublim  potet            hru_intcpevap

# PRMS5
# 0.67442185478284955 0.0      0.0        0.674421847  0.67442185478284955  0.50          1.42531260E-03   0.0

# PRMS6
# 0.67442185478284955 0.0      0.0        0.674421847  0.67442185478284955  0.50          1.42531248E-03   0.0

pk_ice = 0.674421847
potet_sublim = 0.5
potet = .00142531260
hru_intcpevap = 0.0
snowcov_area = 1.0
snow_evap = 0.0
pkwater_equiv = 0.67442185478284955

# %%
# ez = potet_sublim * potet * snowcov_area - hru_intcpevap
ez = potet_sublim * potet * snowcov_area - hru_intcpevap
ez

# %%
pk_ice -= ez
pk_ice

# %%
# cal = pk_temp * ez * 1.27
# pk_def = pk_def + cal

# %%
pkwater_equiv -= ez
pkwater_equiv

# %%
snow_evap = ez
snow_evap

# %%
# pkwater_eqv                freeh2o          snow_evap       pk_ice           pss
# PRMS5
# 0.67370919848326594        0.00000000       7.12656300E-04  0.673709214      0.67442185478284955
#                                             0.000712656300  
# PRMS6
# 0.67370919854147360        0.00000000       7.12656241E-04  0.673709214      0.67442185478284955
#                                             0.000712656241
# Python
# 0.6737091985428495                          0.00071265624   0.673709191

# %%
#         elh = (597.3 - (0.5653 * tavg(chru))) * 2.54
#         this%potet(chru) = jh_coef_2d(chru, curr_month) * (this%tavg_f(chru) - &
#                               jh_coef_hru(chru)) * swrad(chru) / elh

elh = (597.3 - (0.5653 * -10.4694)) * 2.54
elh

# %%
# tavgc            tavgf        swrad          potet            elh
# PRMS5
# -21.7611103      -7.16999960  43.7490540     1.42531260E-03   1548.38794
# PRMS6
# -21.7611122      -7.17000198  43.7490540     1.42531248E-03   1548.38794

# (temperature - 32.0_dp) / 1.8_dp
# -21.7611122      -16.5000000      -27.0222225
tmaxf = 2.3
tminf = -16.64
tavgf = (tmaxf + tminf) / 2.0
tavgc = (tavgf - 32.0) / 1.8

tavgc = -21.7611103
tavgf = -7.1699996

jh_coef = 0.00166199997
# jh_coef = 0.00486
jh_coef_hru = -37.5222015
swrad = 43.7490540

elh = (597.3 - (0.5653 * tavgc)) * 2.54
potet_calc = jh_coef * (tavgf - jh_coef_hru) * swrad / elh
print(tavgf)
print(tavgc)
print(elh)
print(potet_calc)

# %%
# PRMS5
#      4.00000019E-03   0.00000000       0.00000000       2.01490000E-02   5.66500006E-03
#  *   3.97734018E-03   0.00000000       4.00000019E-03   2.01490000E-02   5.66500006E-03

# PRMS6
#      4.00000019E-03   0.00000000       0.00000000       2.01490000E-02   5.66500006E-03
#  *   4.00000019E-03   0.00000000       0.00000000       2.01490000E-02   5.66500006E-03

netrain = 4.00000019E-03
netsnow = 0.0
intcpstor = 0.0
stor = 2.01490000E-02
cov = 5.66500006E-03

# %%
#       Net_precip = Precip*(1.0-Cov)

#       Intcp_stor = Intcp_stor + Precip

#       IF ( Intcp_stor>Stor_max ) THEN
#         Net_precip = Net_precip + (Intcp_stor-Stor_max)*Cov
#         Intcp_stor = Stor_max
#       ENDIF

net_precip = netrain * (1.0 - cov)
intcpstor += netrain
print(net_precip)
print(netrain)

# %%
# pkwater_equiv             freeh2o         temp         esv          sw    trd           cals  hru_ppt   canopy_covden
# PRMS5
# 1.7932620721694548E-006   3.04461736E-07  -3.12500072  0.757000029  0.0   0.579653502   0.0   
# PRMS6
# 1.7932620721694548E-006   3.04461736E-07  -3.12500072  0.757000029  0.0   0.579653502   0.0   0.0       5.66500006E-03

cst = 6.25514030
pk_temp = 0.0
freeh2o = 3.04461736E-07
niteda = 1
sw = 0.0
cecsub = 0.0
temp = -3.12500072
hru_ppt = 0.0
canopy_covden = 5.66500006E-03
tstorm_mo = 0
esv = 0.757000029
emis_noppt = 0.757

air = 0.585e-7 * ((temp + 273.16)**4.0)
ts = 0.0
sno = 325.7

if temp < 0.0:
    ts = temp
    sno = air
    

emis = esv
if hru_ppt > 0 and tstorm_mo == 1:
    if niteda == 1:
        emis = 0.85
        if trd > 1./3.:
            emis = emis_noppt
    else:
        # ! Daytime
        if trd > 1./3.:
            emis = 1.29 - (0.882 * trd)
        elif trd >= 0.5:
            emis = 0.95 - (0.2 * trd)
            
can = canopy_covden * (air - sno)
sky = (1.0 - canopy_covden) * ((emis * air) - sno)

cal = sky + can + cecsub + sw

print(ts)
print(can)
print(sky)
print(cal)

# %%
print('hi')

# %%
qcond = cst * (ts - pk_temp)
print(qcond)

# %%
calnd = freeh2o * 203.2
dif = cal + calnd
print(calnd)
print(dif)

# %%
#      srp              net_rain         snowmelt
# $1   0.00000000       7.90000036E-02   8.64102319E-02

#        soil_moist       smidx_coef      smidx_exp
#        1.16176057       6.52149990E-02  0.276773006

# out
# srp
# 1.21168327E-02    1.21168364E-02
# .0121168327       .0121168364

net_rain = 7.90000036E-02
snowmelt = 8.64102319E-02
soil_moist = 1.16176057
smidx_coef = 6.52149990E-02
smidx_exp = 0.276773006
carea_max = 0.59518

#           smidx = soil_moist(idx) + (0.5 * ptc)
#           ca_fraction = smidx_coef(idx) * 10.0**(smidx_exp(idx) * smidx)

smidx = soil_moist + (0.5 * snowmelt)
ca_fraction = smidx_coef * 10.0**(smidx_exp * smidx)
print(smidx)
print(ca_fraction)

# srpp = ca_fraction * pptp
srpp = ca_fraction * net_rain
print(srpp)

srp = srpp
print(srp)

# %%
#       if (avail_potet < NEARZERO) then
#         this%et_type = 1
#         avail_potet = 0.0
#       elseif (transp_on == 0) then
#         if (snow_free < 0.01) then
#           this%et_type = 1
#         else
#           this%et_type = 2
#         endif
#       elseif (cov_type > 0) then
#         this%et_type = 3
#       elseif (snow_free < 0.01) then
#         this%et_type = 1
#       else
#         this%et_type = 2
#       endif

if transp_on == 0:
    if snow_free < 0.01:
        et_type = 1
    else:
        et_type = 2
elif cov_type > 0:
    et_type = 3
elif snow_free < 0.01:
    et_type = 1
else:
    et_type = 2
    
print(et_type)

# %% [markdown]
# ### szactet work

# %%
#      soil_moist       soil_rechr        avail_potet      soil_moist_max   soil_rechr_max   snow_free
# prms5
# (^)  0.830617666      0.826172769       1.22456090E-03   2.46563506       2.44951963       1.00000000
# prms6
# (^)  0.830617666      0.826172769       1.22456090E-03   2.46563506       2.44951963       1.00000000

soil_moist = 0.830617666
soil_rechr = 0.826172769
avail_potet = 1.22456090E-03
soil_moist_max = 2.46563506
soil_rechr_max = 2.44951963
snow_free = 1
transp_on = 1
cov_type = 1
soil_type = 2

if transp_on == 0:
    if snow_free < 0.01:
        et_type = 1
    else:
        et_type = 2
elif cov_type > 0:
    et_type = 3
elif snow_free < 0.01:
    et_type = 1
else:
    et_type = 2
    
print(et_type)

# Assuming et_type > 1
pcts = soil_moist / soil_moist_max
pctr = soil_rechr / soil_rechr_max
potet_lower = avail_potet
potet_rechr = avail_potet

print(potet_lower, potet_rechr)

if soil_type == 1:
    if pcts < 0.25:
        potet_lower = 0.5 * pcts * avail_potet
    if pctr < 0.25:
        potet_rechr = 0.5 * pctr * avail_potet
elif soil_type == 2:
    if pcts < 0.5:
        potet_lower = pcts * avail_potet
    if pctr < 0.5:
        potet_rechr = pctr * avail_potet
elif soil_type == 3:
    if pcts < 2./3. and pcts > 1./3.:
        potet_lower = pcts * avail_potet
    elif pcts <= 1./3.:
        potet_lower = 0.5 * pcts * avail_potet
    
    if pctr < 2./3. and pctr > 1./3.:
        potet_rechr = pctr * avail_potet
    elif pctr <= 1./3.:
        potet_rechr = 0.5 * pctr * avail_potet

print(pcts, pctr)
print(potet_lower, potet_rechr)

# %%
print(potet_rechr, soil_rechr, soil_moist)
if et_type == 2:
    potet_rechr *= snow_free

if potet_rechr > soil_rechr:
    potet_rechr = soil_rechr
    soil_rechr = 0.0
else:
    soil_rechr = soil_rechr - potet_rechr
    
if et_type == 2 or potet_rechr >= potet_lower:
    if potet_rechr > soil_moist:
        potet_rechr = soil_moist
        soil_moist = 0.0
    else:
        soil_moist -= potet_rechr
    et = potet_rechr
elif potet_lower > soil_moist:
    et = soil_moist
    soil_moist = 0.0
else:
    soil_moist -= potet_lower
    et = potet_lower

if soil_rechr > soil_moist:
    soil_rechr = soil_moist
    
print(potet_rechr, soil_rechr, soil_moist, et)


# %%

# %% [markdown]
# ### precip form

# %%
tmax_f = 34.02
tmax_c = (tmax_f - 32.0) / 1.8
tmin_f = 20.04
tmin_c = (tmin_f - 32.0) / 1.8

print(tmax_f, tmax_c)
print(tmin_f, tmin_c)

# %%
tdiff_f = tmax_f - tmin_f
tdiff_c = tmax_c - tmin_c

print(tdiff_f, tdiff_c)

# %%
print(tdiff_f / 1.8)
print(tdiff_f * (5./9.))

# %%
# tmax_allsnow_c = -0.330614537
tmax_allsnow_f = 31.404894
tmax_allsnow_c = (tmax_allsnow_f - 32.0) / 1.8

print(tmax_allsnow_c, tmax_allsnow_f)

# %%
prmx_f = (tmax_f - tmax_allsnow_f) / tdiff_f * 1.0
prmx_c = (tmax_c - tmax_allsnow_c) / tdiff_c * 1.0
print(prmx_f, prmx_c)

# %%
import sys
sys.float_info

# %%
tmax_c = 1.5000000000000016
tmax_f = tmax_c * 1.8 + 32.0

print(tmax_c, tmax_f)

# %%
tmax_allsnow_c = -1.2080521053738065
tmax_allsnow_f = tmax_allsnow_c * 1.8 + 32.0

tmax_allrain_c = 1.5000004238552518
tmax_allrain_f = tmax_allrain_c * 1.8 + 32.0

print(tmax_allsnow_c, tmax_allsnow_f)

print(tmax_allrain_c, tmax_allrain_f)

# %%
tmax_allsnow = 29.825507
tmax_allrain_offset = 4.874493

print(tmax_allsnow + tmax_allrain_offset)

# %%
