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
# RADIANS=PI/180.0D0
# real(r64), parameter :: ECCENTRICY = 0.01671_dp
# real(r64), parameter :: DAYSYR = 365.242_dp
# ! 0.016723401  daily change -1.115E-09, eccen = 0.016723401 + (julhour-julhour(1966,1,0,18))+dmin/60)/24*-1.115E-09
# # ! julday(1966,1,0.75 UT) = 2439126.25
# # ! eccen = 0.01675104-0.00004180*T-0.000000126*T^2  T is julian centuries (days time from epoch, is GMT from Jan 0.0
# real(r64), parameter :: DEGDAY = 360.0_dp / DAYSYR
# real(r64), parameter :: DEGDAYRAD = DEGDAY * RADIANS ! about 0.00143356672

# integer(i32) :: i
# real(dp), parameter :: JD(366) = [(dble(i), i=1,366)]
# real(dp), parameter :: Y(366) = (JD - 1.0_dp) * DEGDAYRAD
# real(dp), parameter :: JD3(366) = (JD - 3.0_dp) * DEGDAYRAD
# real(dp), parameter :: Y2(366) = Y * 2.0_dp
# real(dp), parameter :: Y3(366) = Y * 3.0_dp
# real(dp), parameter :: obliquity(366) = 1.0_dp - (ECCENTRICY * cos(JD3))

# this%solar_declination = 0.006918_dp - &
#                        0.399912_dp * COS(Y) + 0.070257_dp * SIN(Y) - &
#                        0.006758_dp * COS(Y2) + 0.000907_dp * SIN(Y2) - &
#                        0.002697_dp * COS(Y3) + 0.00148_dp * SIN(Y3)

PI_12 = 12.0 / np.pi
RADIANS = np.pi / 180.0
ECCENTRICITY = 0.016171
DAYSYR = 365.242
DEGDAY = 360.0 / DAYSYR
DEGDAYRAD = DEGDAY * RADIANS

JD = np.arange(1.0, 12.0, 1)
Y = (JD - 1.0) * DEGDAYRAD
JD3 = (JD - 3.0) * DEGDAYRAD
Y2 = Y * 2.0
Y3 = Y * 3.0
obliquity = 1.0 - (ECCENTRICITY * np.cos(JD3))

solar_declination = 0.006918 - 0.399912 * np.cos(Y) + 0.070257 * np.sin(Y) - 0.006758 * np.cos(Y2) + 0.000907 * np.sin(Y2) - 0.002697 * np.cos(Y3) + 0.00148 * np.sin(Y3)

print('PI_12: {}'.format(PI_12))
print('RADIANS: {}'.format(RADIANS))
print('eccentricity: {}'.format(ECCENTRICITY))
print('daysyr: {}'.format(DAYSYR))
print('degday: {}'.format(DEGDAY))
print('degdayrad: {}'.format(DEGDAYRAD))
# print('solar_declination: {}'.format(solar_declination))
# print('jd: {}'.format(JD))
# print('y: {}'.format(Y))

# %%
def compute_t(lat, solar_declination):
    # lat - latitude in radians
    # solar_declination - array of 366 solar declination values
    
    tx = -np.tan(lat) * np.tan(solar_declination)

#     if tx[tx < -1.0].any():
#         print('Some values < -1.0')
#     if tx[tx > 1.0].any():
#         print('Some values > 1.0')
        
    tx[tx < -1.0] = math.pi
    tx[tx > 1.0] = 0.0
    tx[(tx >= -1.0) & (tx <= 1.0)] = np.arccos(tx)
    return tx


#     if (tx < -1.0):
#         return math.pi
#     elif (tx > 1.0):
#         return 0.0
#     else:
#         return math.acos(tx)

# %%
for xx in range(-180, 180):
    bb = compute_t(float(xx), solar_declination)


# %%
# call this%compute_soltab(this%hru_cossl(chru), this%soltab_horad_potsw(:, chru), &
#                          this%soltab_sunhrs(:, chru), obliquity, &
#                          this%solar_declination, 0.0, 0.0, &
#                          hru_lat(chru), hru_type(chru), chru)

# %%
def compute_soltab(latitude, slope, aspect):
    aspect_radians = aspect * RADIANS
    x0_radians = latitude * RADIANS
    
    # Latitude of the equivalent slope
    x1_radians = np.arcsin(np.cos(slope) * np.sin(x0_radians) + np.sin(slope) * np.cos(x0_radians) * np.cos(aspect_radians))
    
    # Denominator of equation 12 (Lee, 1963)
    d1_radians = np.cos(slope) * np.cos(x0_radians) - np.sin(slope) * np.sin(x0_radians) * np.cos(aspect_radians)
    
    # Difference in longitude between the location of the HRU and equivalent
    # horizontal surface expressed in angle hours. Equation 12 (Lee, 1963).
    x2_radians = np.arctan(np.sin(slope) * np.sin(aspect_radians) / d1_radians)
    if d1_radians < 0.0:
        x2_radians = x2_radians + np.pi
        
    r0 = 2.0
    r1 = (60.0 * 2.0) / np.power(obliquity, 2.0)
   
    print(f'latitude: {latitude}')
    print(f'slope: {slope}')
    print(f'aspect: {aspect}')
    print(f'x0: {x0_radians}')
    print(f'x1: {x1_radians}')
    print(f'd1: {d1_radians}')
    print(f'x2: {x2_radians}')
    print(f'r0: {r0}')
    print(f'r1: {r1}')
    
    t = compute_t(x1_radians, solar_declination)
    t7 = t - x2_radians
    t6 = -t - x2_radians
    print(f't7: {t7}')
    print(f't6: {t6}')
    
    t = compute_t(x0_radians, solar_declination)
    t1 = t
    t0 = -t
    print(f't1: {t1}')
    print(f't0: {t0}')
    
    t3 = np.where(t7 > t1, t1, t7)
    t2 = np.where(t6 < t0, t0, t6)
    print(f't3: {t3}')
    print(f't2: {t2}')
    
#      2.220446049250313E-016
    if np.abs(np.arctan(slope) < 2.220446049250313e-16):
        print('atan of slope near zero')
        solt = func3(0.0, x0_radians, t1, t0, r1, solar_declination)
        sunh = (t1 - t0) * PI_12
    else:
        m = t3 < t2
        print(m)
        t2 = np.where(m, 0.0, t2)
        t3 = np.where(m, 0.0, t3)
        t6 += 2.0 * np.pi
        
        m = t6 < t1
        print(m)
        
        t7 -= 2.0 * np.pi
        print(t7 > t0)


# %%
def func3(v, w, x, y, r1, solar_declination):
    return r1 * PI_12 * (np.sin(solar_declination) * np.sin(w) * (x - y) + np.cos(solar_declination) * np.cos(w) * (np.sin(x+v) - np.sin(y+v)))


# %%
compute_soltab(60.0, 1., 0.6)

# %%
np.arccos(1.557)

# %%
math.cos(1.557)

# %%
