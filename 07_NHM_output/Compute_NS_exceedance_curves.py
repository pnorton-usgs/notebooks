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
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
from collections import OrderedDict
import geopandas as gpd
import numpy as np
import pandas as pd

# %%
calibrations = ['byHRU', 'byHRU_musk', 'byHRU_musk_obs', 'byHRU_noroute', 'byHRU_noroute_obs']
output_suffix = {'Ref': 'ref', 'Non-ref': 'nonref'}

workdir = '/Volumes/USGS_NHM1/calibrations/NHMv10/DAYMET_releases/model_stats'

falc_cls = 'Non-ref'  # either: Ref or Non-ref

out_filename = f'{workdir}/ns_curves_{output_suffix[falc_cls]}.csv'


# %%
def get_exceedance_curve(ts, cal):
    stat_name = f'{ts.name}_{cal}'
    
    ranked_stat = sorted(ts[ts.notnull()])
    prob = np.arange(len(ranked_stat), dtype=float) + 1.0
    prob = (prob / (len(ranked_stat) + 1.0))
    
    # Return dataframe of exceedence curve values
    df = pd.DataFrame({'exceedance': prob, stat_name: ranked_stat}, columns=['exceedance', stat_name])
    df.set_index('exceedance', inplace=True)
    return df


# %% [markdown]
# ## Read the gage statistics files

# %%
df_dict = OrderedDict()

for cc in calibrations:
    df_dict[cc] = pd.read_csv(f'{workdir}/gage_stats_{cc}.csv', sep=',')
    df_dict[cc] = df_dict[cc][df_dict[cc]['falcone_class'] == falc_cls]
    
    # Remove rows where NS is NaN
    df_dict[cc] = df_dict[cc][df_dict[cc]['ns'].notna()]

# %% [markdown]
# ## Compute the Nash-Sutcliffe exceedance curves

# %%
first = True

for cc in calibrations:
    if first == True:
        ns_crv = get_exceedance_curve(df_dict[cc].loc[:, 'ns'], cc)
        nslog_crv = get_exceedance_curve(df_dict[cc].loc[:, 'nslog'], cc)
        mon_ns_crv = get_exceedance_curve(df_dict[cc].loc[:, 'mon_ns'], cc)
        first = False
    else:
        ns_crv = pd.merge(ns_crv, get_exceedance_curve(df_dict[cc].loc[:, 'ns'], cc), how='left',
                          left_index=True, right_index=True)
        nslog_crv = pd.merge(nslog_crv, get_exceedance_curve(df_dict[cc].loc[:, 'nslog'], cc), how='left',
                             left_index=True, right_index=True)
        mon_ns_crv = pd.merge(mon_ns_crv, get_exceedance_curve(df_dict[cc].loc[:, 'mon_ns'], cc), how='left',
                              left_index=True, right_index=True)

# %% [markdown]
# ## Merge the curves into a single dataframe

# %%
df = ns_crv
df = pd.merge(df, nslog_crv, how='left', left_index=True, right_index=True)
df = pd.merge(df, mon_ns_crv, how='left', left_index=True, right_index=True)

# Remove rows containing -inf
# df = df[np.isfinite(df).all(1)]

# %% [markdown]
# ## Write the results to a file

# %%
df.to_csv(out_filename, index=True)

# %%
