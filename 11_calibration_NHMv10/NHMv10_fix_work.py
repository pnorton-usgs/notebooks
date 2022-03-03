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
import numpy as np

# from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile
# from pyPRMS.plot_helpers import get_projection

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/tmp/20200831_model_check'
outdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/parameters'



byHRU_file = f'{workdir}/myparam.param.byHRU'
byHW_file = f'{workdir}/myparam.param.byHW'
byOBS_file = f'{workdir}/myparam.param.byOBS'

# %%
byHRU_pf = ParameterFile(byHRU_file, verbose=True, verify=True)

# %%
bh_radmax = byHRU_pf.parameters.get('radmax').tolist()
bh_radmax.sort()

# %%
byHW_pf = ParameterFile(byHW_file, verbose=True, verify=True)
bhw_radmax = byHW_pf.parameters.get('radmax').tolist()
bhw_radmax.sort()

# %%
byOBS_pf = ParameterFile(byOBS_file, verbose=True, verify=True)
bo_radmax = byOBS_pf.parameters.get('radmax').tolist()
bo_radmax.sort()

# %%
for xx, yy, zz in zip(bh_radmax, bhw_radmax, bo_radmax):
    print(xx, yy, zz)

# %%
len(set(bh_radmax))

# %%
len(set(bhw_radmax))

# %%
len(set(bo_radmax))

# %%
len(set(bh_radmax).difference(set(bhw_radmax)))

# %%
len(set(bhw_radmax).difference(set(bo_radmax)))

# %%
len(set(bh_radmax).difference(set(bo_radmax)))

# %%
byHRU_pf.write_netcdf(f'{outdir}/B_byHRU_from_src.nc')

# %%
byHW_pf.write_netcdf(f'{outdir}/C_byHW_from_src.nc')

# %%
byOBS_pf.write_netcdf(f'{outdir}/D_byOBS_from_src.nc')

# %%

# %%

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/paramdb_v10/maurer_fix_work/byHRU_musk_final_params_compare/FINALparams'
pfile = f'{workdir}/FINALparams_0505'

# %%
pf = ParameterFile(filename=pfile, verbose=False, verify=True)

# %%
nhm_id = pf.parameters.get('nhm_id').data
nhm_seg = pf.parameters.get('nhm_seg').data

cparam = 'carea_max'

# %%
for cid, vals in zip(nhm_seg, pf.parameters.get(cparam).data):
    print(cid, vals)
#     pdb_orig.parameters.update_element(cparam, cid, vals)

# %%
src_data = pf.parameters.get(cparam).data

# %%
np.where(nhm_id == 6918)[0]

# %%
new_vals = [0.5, 0.5, 0.6, 0.5, 0.691814, 0.685925,
            0.685925, 0.5, 0.5, 0.5, 0.5, 0.5]

# %%
np.array_equal(src_data[0], new_vals)

# %%
cc = np.where(src_data == 0.252069)[0]
# cc = np.where(nhm_id == 6922)

# %%
src_data.shape

# %%
type(cc)

# %%
nhm_id

# %%
aa = pf.parameters.get('rad_trncf')

# %%
aa.data

# %%
aa.data[0]

# %%
np.array_equal(aa.data[0], 0.586339)

# %%
