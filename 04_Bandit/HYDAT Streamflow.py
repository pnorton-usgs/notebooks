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
#     display_name: Python [conda env:gis_38]
#     language: python
#     name: conda-env-gis_38-py
# ---

# %%

# %%
import calendar
import datetime
import pandas as pd
import numpy as np
import sqlite3

# %%
db_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/HYDAT'
db_name = 'Hydat.sqlite3'

con = sqlite3.connect(f'{db_dir}/{db_name}')
cur = con.cursor()

# %%
# df = pd.read_sql_query('SELECT * FROM ANNUAL_STATISTICS where DATA_TYPE=="Q";', con)
df = pd.read_sql_query('SELECT * FROM STATIONS;', con)

# %%
df.head()

# %%
