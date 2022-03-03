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
# %matplotlib inline

import matplotlib.pyplot as plt
import numpy as np

# %%

# %%
# fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(17,11))
# ax = axes
# ax.plot(sim_mon_mn.index.to_pydatetime(), sim_mon_mn, linewidth=.25, color='#cc3333', label='MOCOM', alpha=0.5)

seg_outflow = 100.0

width_alpha = 10.0
width_m = np.arange(0.0001, 9.0, 0.001)
# width_alpha = np.arange(0.0001, 2.0, .001)
# width_m = 0.015

seg_width_flow = width_alpha * seg_outflow**width_m

# Seg_width_flow(i) = width_alpha(i) * sngl(Seg_outflow(i)) ** width_m(i)
plt.plot(width_m, seg_width_flow)
plt.ylabel('seg_width_flow')
plt.xlabel('width_m')
plt.show()

# %%
print(width_alpha)

# %%
width_alpha.size

# %%
seg_width_flow.size

# %%
#             tdiff = max(tmax(chru) - tmin(chru), 0.0001)

#             this%prmx(chru) = ((tmax(chru) - tmax_allsnow(chru, month)) / tdiff) * rainmix_adj(chru, month)

# #             ! Make sure prmx is between 0.0 and 1.0
#             this%prmx(chru) = min(max(0.0, this%prmx(chru)), 1.0)

# %%

# %% [markdown]
# ### prmx tests 

# %%
tmax_as = 31.
tmax_ar = tmax_as + 3.5
tm = 20.
adj = 1.0

for tx in range(20, 50):
    tdiff = max(float(tx) - tm, 0.0001)
    prmx = ((float(tx) - tmax_as) / tdiff) * adj
    
    if float(tx) < tmax_as:
        print('{} {:.4f} {:0.3f} ALL-SNOW'.format(tx, tdiff, prmx))
    elif tm > tmax_as or float(tx) >= tmax_ar:
        print('{} {:.4f} {:0.3f} ALL-rain'.format(tx, tdiff, prmx))
#         print(tdiff, prmx, 'tmin > tmax_as OR tmax >= tmax_ar')
    else:
        if prmx < 1.0:
            print('{} {:.4f} {:0.3f} mixed'.format(tx, tdiff, prmx))
        else:
            print('{} {:.4f} {:0.3f} ALL-rain2'.format(tx, tdiff, prmx))


# %%
