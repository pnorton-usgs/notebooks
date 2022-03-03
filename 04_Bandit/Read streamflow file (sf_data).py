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

# %%
from pyPRMS.Streamflow import Streamflow
import datetime
import pandas as pd

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/notebooks/04_Bandit'
filename = '{}/sf_data'.format(workdir)

st_date = datetime.datetime(2015,9,25)
en_date = datetime.datetime(2015,10,2)

# %%
sf = Streamflow(filename)

# %%

# %%
df = sf.data[st_date:en_date]

# %%
df.head()

# %%
df.info()

# %%
df2 = df.iloc[:,128]

# %%
df2

# %%
# df[df['A'].str.contains("Hello|Britain")==True]
df2[df2.str.contains('_Dis')].head(1).keys()[0].strftime('%Y-%m-%d')

# %%
df2.str.contains('_Dis').sum()

# %%
df2.replace('_Dis', '', regex=True)

# %%
df3 = pd.DataFrame(df2)
df3

# %%
df3['02450000'].str.contains('_Dis').index[0].strftime('%Y-%m-%d')

# %%
