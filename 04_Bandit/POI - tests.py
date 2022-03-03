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

from datetime import datetime
from Bandit.points_of_interest import POI

# %%
src_pois = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/poi_data/*_pois.nc'

# gage_list = ['01AD001', '11AA005']
# gage_list = ['01010070', '01018000']
gage_list = ['01010070', '11AA005']
st_date = datetime(2000, 1, 1)
en_date = datetime(2001, 12, 31)

# %%
pp = POI(src_path=src_pois, st_date=st_date, en_date=en_date, gage_ids=gage_list, verbose=True)

# %%
# pp.write_ascii('crap.csv')
pp.write_netcdf('crap.nc')

# %%
# pp.data['discharge'].loc[:, '11AA005']

# pp.data['discharge'].loc[dict(poi_id=b'11AF005')]

# pp.data.sel(poi_id=b'11AF005')
dd = pp.get('poi_id') # .tolist()

# %%
dd.tolist()
# np.array(dd).astype('S')

# %%
pp.data['poi_id']

# %%
'time' in pp.data.dims

# %%
pp.data['discharge'].loc[['01AD001', '11AF005'], st_date:en_date]

# %%
pp.data.info()


# %%
pp.data['poi_id'].dtype

# %%
pp.data['poi_name']

# %%
pp.data.info()

# %%
xx = ['011A500', 'fff3234']

# %%
yy = np.array(xx, dtype='S10')

# %%
yy

# %%
yy.shape

# %%
