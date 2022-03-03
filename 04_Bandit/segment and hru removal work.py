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

# %%
import numpy as np
from collections import OrderedDict

# %% [markdown]
# ### Removing stream segments

# %%
nhm_seg = [30113, 30114, 30115, 30116, 30117, 30118, 30119]
tosegment = [6, 4, 7, 5, 1, 3, 0]
tosegment_nhm = [30118, 30116, 30119, 30117, 30113, 30115, 30120]
poi_gage_segment = [7]


# %%
# create nhm_seg version of poi_gage_segment
poi_gage_segment_nhm = [nhm_seg[ss-1] for ss in poi_gage_segment]
print(poi_gage_segment_nhm)

# %%
# create map connecting nhm_seg ids with tosegment_nhm ids
aa = OrderedDict((seg, toseg) for seg, toseg in zip(nhm_seg, tosegment_nhm))
print(aa)

# %%
new_nhm_seg = [kk for kk in aa.keys()]
print(new_nhm_seg)

new_tosegment = [new_nhm_seg.index(kk)+1 if kk in new_nhm_seg else 0 for kk in aa.values()]
print(new_tosegment)

new_tosegment_nhm = [kk for kk in aa.values()]
print(new_tosegment_nhm)

# %%
# Remove a segment
bb = aa.copy()
del bb[30116]
print(bb)

# %%
# Recreate arrays with new segment lists
new_nhm_seg = [kk for kk in bb.keys()]
print(new_nhm_seg)

new_tosegment = [new_nhm_seg.index(kk)+1 if kk in new_nhm_seg else 0 for kk in bb.values()]
print(new_tosegment)

new_tosegment_nhm = [kk for kk in bb.values()]
print(new_tosegment_nhm)

# %%
# Recreate arrays with new segment lists
new_nhm_seg = [kk for kk in bb.keys()]
print(new_nhm_seg)

new_nhm_seg_dict = {}
for ii, ss in enumerate(new_nhm_seg):
    new_nhm_seg_dict[ss] = ii+1
    
new_tosegment_nhm = [kk for kk in bb.values()]
print(new_tosegment_nhm)

new_tosegment = [new_nhm_seg_dict[kk] if kk in new_nhm_seg_dict else 0 for kk in bb.values()]
# new_tosegment = [new_nhm_seg.index(kk)+1 if kk in new_nhm_seg else 0 for kk in bb.values()]
print(new_tosegment)



# %%
# Update the poi_gage_segment with new indices
new_poi_gage_segment = [new_nhm_seg.index(xx)+1 for xx in poi_gage_segment_nhm]
print(new_poi_gage_segment)

# %% [markdown]
# ### Removing HRUs

# %%
nhm_id = [57874, 57875, 57878, 57881, 57868, 57873, 57879, 57880, 57882, 57883, 57869, 57870, 57864, 57865, 88888, 999]
# hru_segment = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7]
hru_segment = [1, 1, 2, 2, 3, 3, 4, 4, 0, 5, 6, 6, 7, 7]
hru_segment_nhm = [30113, 30113, 30114, 30114, 30115, 30115, 30116, 30116, 0, 30117, 30118, 30118, 30119, 30119, 0, 0]
# hru_segment_nhm = [30113, 30113, 30114, 30114, 30115, 30115, 30116, 30116, 30117, 30117, 30118, 30118, 30119, 30119]
hru_deplcrv = [7, 8, 9, 12, 3, 6, 10, 11, 13, 14, 4, 5, 1, 2]
nonrouted = [88888]

# %%
# Create map connecting nhm_id to hru_segment_nhm ids
cc = OrderedDict((nhm, hru_seg) for nhm, hru_seg in zip(nhm_id, hru_segment_nhm))
print(cc)
print(cc.values())

# %%
# Can we re-construct original arrays/order from the ordered dictionary? Yes
new_nhm_id = [xx for xx in cc.keys()]
print(new_nhm_id)

new_hru_segment_nhm = [vv if vv in new_nhm_seg else vv if kk in nonrouted else 0 if kk == 0 else -1 for kk, vv in cc.items()]
# new_hru_segment_nhm = [kk if kk in new_nhm_seg else 0 if kk == 0 else -1 for kk in cc.values()]
print(new_hru_segment_nhm)

new_hru_segment = [new_nhm_seg.index(kk)+1 if kk in new_nhm_seg else 0 if kk == 0 else -1 for kk in cc.values()]
print(new_hru_segment)

# %%
new_hru_segment_nhm = []

for kk in cc.values():
    if kk in new_nhm_seg:
        new_hru_segment_nhm.append(kk)
    elif kk == 0:
        new_hru_segment_nhm.append(kk)
    else:
        new_hru_segment_nhm.append(-1)

# %%
# Remove a few HRUs
dd = cc.copy()
del dd[57879]
del dd[57880]
print(dd)

# %%
# Re-construct the arrays
new_nhm_id = [xx for xx in dd.keys()]
print(new_nhm_id)

new_hru_segment_nhm = [kk if kk in new_nhm_seg else 0 if kk == 0 else -1 for kk in dd.values()]
print(new_hru_segment_nhm)

new_hru_segment = [new_nhm_seg.index(kk)+1 if kk in new_nhm_seg else 0 if kk == 0 else -1 for kk in dd.values()]
print(new_hru_segment)

# %%
type(new_nhm_id)

# %%
nhm_idx = OrderedDict((id, ii) for ii, id in enumerate(nhm_id))
print(nhm_idx)

# %%
for xx in list(nhm_idx.keys()):
    if xx in [57879, 57880]:
        del nhm_idx[xx]
        
print(nhm_idx)

# %%
xx = [hru_segment_nhm[yy] for yy in nhm_idx.values()]
print(xx)

# %%
# [nhm_seg.index(kk)+1 if kk in nhm_seg else 0 if kk==0 else -1
#                                             for kk in self.get('hru_segment_nhm').data.tolist()]
yy = [new_nhm_seg.index(kk)+1 if kk in new_nhm_seg else 0 if kk==0 else -1 for kk in xx]
print(yy)

# %%
type(nhm_idx.values())

# %%
hru_deplcrv = np.array([7, 8, 9, 12, 3, 6, 10, 11, 13, 14, 4, 5, 1, 2])
snarea_curve_tmp = [0.047029, 0.225739, 0.376232, 0.498508, 0.611377, 0.705436, 0.771276, 0.827711, 0.87474, 0.931175, 0.940581, 0.048659, 0.233563, 0.389271, 0.515784, 0.632565, 0.729883, 0.798006, 0.856396, 0.905055, 0.963446, 0.973177, 0.048554, 0.233059, 0.388432, 0.514673, 0.631202, 0.72831, 0.796286, 0.85455, 0.903105, 0.961369, 0.97108, 0.046454, 0.222979, 0.371632, 0.492412, 0.603902, 0.696809, 0.761845, 0.81759, 0.864044, 0.919789, 0.929079, 0.046565, 0.223512, 0.37252, 0.493588, 0.605344, 0.698474, 0.763665, 0.819543, 0.866108, 0.921986, 0.931299, 0.0467, 0.224161, 0.373602, 0.495022, 0.607103, 0.700503, 0.765884, 0.821924, 0.868624, 0.924664, 0.934004, 0.047971, 0.230259, 0.383765, 0.508489, 0.623618, 0.71956, 0.786719, 0.844284, 0.892254, 0.949819, 0.959413, 0.048331, 0.231991, 0.386651, 0.512313, 0.628308, 0.724971, 0.792635, 0.850632, 0.898964, 0.956962, 0.966628, 0.04808, 0.230786, 0.384643, 0.509652, 0.625045, 0.721206, 0.788518, 0.846215, 0.894295, 0.951991, 0.961608, 0.04734, 0.227234, 0.378723, 0.501808, 0.615425, 0.710106, 0.776383, 0.833191, 0.880531, 0.93734, 0.946808, 0.047468, 0.227847, 0.379745, 0.503162, 0.617086, 0.712022, 0.778477, 0.835439, 0.882907, 0.939869, 0.949363, 0.047643, 0.228685, 0.381142, 0.505013, 0.619355, 0.71464, 0.78134, 0.838512, 0.886154, 0.943326, 0.952854, 0.04858, 0.233184, 0.38864, 0.514948, 0.63154, 0.728699, 0.796711, 0.855007, 0.903587, 0.961883, 0.971599, 0.047474, 0.227873, 0.379788, 0.503219, 0.617156, 0.712103, 0.778566, 0.835534, 0.883008, 0.939976, 0.949471]

uniq_deplcrv = list(set(hru_deplcrv))
print(uniq_deplcrv)

# %%
new_hru_deplcrv = [uniq_deplcrv.index(cc)+1 for cc in hru_deplcrv]
print(new_hru_deplcrv)

# %%
type(hru_deplcrv)

# %%
