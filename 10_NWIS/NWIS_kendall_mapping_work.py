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
#     display_name: Python [conda env:bandit]
#     language: python
#     name: conda-env-bandit-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import numpy as np
import pandas as pd

from pyPRMS.ParameterFile import ParameterFile

# %%
workdir = '/Users/pnorton/Projects/Streamflow_CONUS/datapull_CONUS_20200123'
filename = f'{workdir}/conus_annual_HUC_ALL_wy1960-2018_kendall.csv'

# Columns for streamgage kendall file
stn_col_names = ['site_no', 'pval', 'tau', 'trend', 'station_nm', 'dec_lat_va', 'dec_long_va',
                 'drain_area_va', 'contrib_drain_area_va']
stn_col_types = [np.str_, np.float_, np.float_, np.int, np.str_, np.float_, np.float_, np.float_, np.float_]

stn_cols = dict(zip(stn_col_names, stn_col_types))

# %%
# Load the kendall results
df_stns = pd.read_csv(filename, sep='\t', usecols=stn_col_names, dtype=stn_cols)
df_stns.set_index('site_no', inplace=True)
df_stns.head()

# %%
jobs_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs'
model_dir = f'{jobs_dir}/20200205_NWIS_wy1960-2018'

# This isn't quite right. The station id directory needs to be included
nhm_id_filename = f'{jobs_dir}/{model_dir}/something'

out_dir = '/Users/pnorton/tmp/nwis_crap'

# %%
# Write streamgages to the model job directory
# (comment out df_stns.set_index('site_no', inplace=True) above first)
df_stns.iloc[:,0].to_csv(f'{model_dir}/streamgages.csv', header=False, index=False)

# %%
stn_list = df_stns.index.tolist()
len(stn_list)

# %%

# %%
# For each station: 
#   open the model
#   read nhm_id
#   create csv file mapping nhm_id to the station trend

# %%
# Upward trends

df_up = df_stns[df_stns['trend'] > 0.0]

print('Upward trends')
for ss in stn_list:
    curr_file = f'{model_dir}/{ss}/myparam.param'
    
    try:
        pfile = ParameterFile(curr_file)
    except FileNotFoundError:
        print(f'  No model exists for {ss}: skipping.')
        continue

    print(ss)
    
    ids = pfile.parameters['nhm_id'].tolist()
    
    try:
        tt = df_up.loc[ss, 'trend']
        ofile = open(f'{out_dir}/{ss}_upward.csv', 'w')

        for ii in ids:
            ofile.write(f'{ii},{tt}\n')

        ofile.close()
    except KeyError:
        continue


# Downward trends
print('Downward trends')
df_down = df_stns[df_stns['trend'] < 0.0]

for ss in stn_list:   
    curr_file = f'{model_dir}/{ss}/myparam.param'
    
    try:
        pfile = ParameterFile(curr_file)
    except FileNotFoundError:
        print(f'  No model exists for {ss}: skipping.')
        continue

    print(ss)
    
    ids = pfile.parameters['nhm_id'].tolist()
    
    try:
        tt = df_down.loc[ss, 'trend']
        ofile = open(f'{out_dir}/{ss}_downward.csv', 'w')

        for ii in ids:
            ofile.write(f'{ii},{tt}\n')

        ofile.close()
    except KeyError:
        continue

# %%
df = pd.read_csv(f'{out_dir}/all_upward_wy1960-2018.csv')
df.head()

# %%
df.info()

# %%
dfg = df.groupby(['nhm_id']).mean()

# %%
dfg.to_csv(f'{out_dir}/all_upward_mean_wy1960-2018.csv', index=True)

# %%
dfg.head()

# %%
df_stns.loc[ss, 'trend']

# %%
df_stns[df_stns['trend'] < 0.0]

# %%
