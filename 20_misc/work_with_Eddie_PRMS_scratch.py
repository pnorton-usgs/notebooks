# ---
# jupyter:
#   jupytext:
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

# %%
import pandas as pd
import numpy as np
import datetime
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.colors as colors
# from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile

import pydot
import networkx as nx

# %% [markdown]
# ## Create accumulations of seginc_* variables to each segment.

# %%

# %%
cvar = 'seginc_gwflow'

base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/src/tests_prms5/tests'
# base_dir = '/Users/pnorton/Projects/National_Hydrology_Model/FY2022_projects/extraction_requests'
filename = f'{cvar}.csv'

plot_filename = f'{base_dir}/seg49294_network.pdf'
param_filename = f'{base_dir}/20220422_NHMv11_byHRU_seg49294/myparam.param'

byHRU_file = f'{base_dir}/20220422_NHMv11_byHRU_seg49294/output/{filename}'
byHW_file = f'{base_dir}/20220422_NHMv11_byHW_seg49294/output/{filename}'
byHW_obs_file = f'{base_dir}/20220422_NHMv11_byHWobs_seg49294/output/{filename}'

# %%
df1 = pd.read_csv(byHRU_file, sep=',', header=0, index_col=0, parse_dates=True)
df1.columns = df1.columns.astype(int)

# %%
df2 = pd.read_csv(byHW_file, sep=',', header=0, index_col=0, parse_dates=True)
df2.columns = df2.columns.astype(int)

# %%
df3 = pd.read_csv(byHW_obs_file, sep=',', header=0, index_col=0, parse_dates=True)
df3.columns = df3.columns.astype(int)

# %%
pfile = ParameterFile(param_filename, verbose=True, verify=True)
params = pfile.parameters

# %%
hru_segment_nhm = params['hru_segment_nhm'].tolist()
hru_segment = params['hru_segment'].tolist()

nhm_id = params['nhm_id'].tolist()
nhm_seg = params['nhm_seg'].tolist()
tosegment_nhm = params['tosegment_nhm'].tolist()


poi = {}
if pfile.parameters.exists('poi_gage_segment'):
    poi_gage_id = pfile.parameters['poi_gage_id'].tolist()
    poi_gage_segment = pfile.parameters['poi_gage_segment'].tolist()
    
    for gg, ss in zip(poi_gage_id, poi_gage_segment):
        poi[ss] = gg

# %% tags=[]
num_outlets = 0
dag_streamnet = nx.DiGraph()

for ii, vv in enumerate(tosegment_nhm):
    if vv == 0 or vv not in nhm_seg:
#     if vv == 0 or tosegment[ii] == 0:
        dag_streamnet.add_node(nhm_seg[ii], style='filled', fontcolor='white', fillcolor='blue')
        dag_streamnet.add_node(vv, style='filled', fontcolor='white', fillcolor='grey')
        dag_streamnet.add_edge(nhm_seg[ii], vv)

        num_outlets += 1
    else:
        dag_streamnet.add_edge(nhm_seg[ii], tosegment_nhm[ii])
        if nhm_seg[ii] > 56460:
            # NOTE: This will only work correctly prior to renumbering the NHM segments
            dag_streamnet.nodes[nhm_seg[ii]]['fillcolor'] = 'deeppink'
            dag_streamnet.nodes[nhm_seg[ii]]['style'] = 'filled'

    if ii + 1 in poi:
        dag_streamnet.nodes[nhm_seg[ii]]['shape'] = 'box'
        dag_streamnet.nodes[nhm_seg[ii]]['style'] = 'filled'
        dag_streamnet.nodes[nhm_seg[ii]]['fillcolor'] = 'gold'
        dag_streamnet.nodes[nhm_seg[ii]]['label'] = f'{nhm_seg[ii]}\n POI: {poi[ii + 1]}'
        
# Output any cycles/loops
# Also update attributes for nodes which are part of a cycle
for xx in nx.simple_cycles(dag_streamnet):
    for yy in xx:
        dag_streamnet.nodes[yy]['style'] = 'filled'
        dag_streamnet.nodes[yy]['fillcolor'] = 'orange'
    print(xx)
    
print("Number of subset nodes: {}".format(dag_streamnet.number_of_nodes()))
print("Number of subset edges: {}".format(dag_streamnet.number_of_edges()))
print("Number of outlet nodes: {}".format(num_outlets))
print("Is a DAG: {}".format(nx.is_directed_acyclic_graph(dag_streamnet)))    

print("Isolate nodes:")
print(list(nx.isolates(dag_streamnet)))

# %%
# Create plot of the graph using pydot
pdot = nx.nx_pydot.to_pydot(dag_streamnet)
pdot.write_pdf(plot_filename)

# %%

# %%
# %%time
df1_cum = df1.copy()
df2_cum = df2.copy()
df3_cum = df3.copy()


# Create upstream graph
dag_us = dag_streamnet.reverse()

# dsmost_segs = list(dag_streamnet.nodes)
dsmost_segs = nhm_seg
# dsmost_segs = [49698]
print(f'Number of nodes: {len(dsmost_segs)}')

for dsmost_seg in dsmost_segs:
    uniq_seg_us = set()
    
    try:
        pred = nx.dfs_predecessors(dag_us, dsmost_seg)
        uniq_seg_us = uniq_seg_us.union(set(pred.keys()).union(set(pred.values())))
    except KeyError:
        print('\nKeyError: Segment {} does not exist in stream network'.format(dsmost_seg))

    dag_ds_subset = dag_streamnet.subgraph(uniq_seg_us).copy()

    node_outlets = [ee[0] for ee in dag_ds_subset.edges()]
#     true_outlets = []
    true_outlets = set([dsmost_seg]).difference(set(node_outlets))

    # Add the downstream segments that exit the subgraph
    for xx in true_outlets:
        nhm_outlet = list(dag_streamnet.neighbors(xx))[0]
        dag_ds_subset.add_node(nhm_outlet, style='filled', fontcolor='white', fillcolor='grey')
        dag_ds_subset.add_edge(xx, nhm_outlet)

    toseg_idx = list(set(xx[0] for xx in dag_ds_subset.edges))
    toseg_idx.insert(0, dsmost_seg)
    
    # csum = 0.0
    # for xx in toseg_idx:
    #     csum += area_by_seg.get(xx, 0.0)
    # fhdl.write(f'{dsmost_seg},{csum}\n')
    
    df1_cum.loc[:, dsmost_seg] = df1.loc[:, toseg_idx].sum(axis=1)
    df2_cum.loc[:, dsmost_seg] = df2.loc[:, toseg_idx].sum(axis=1)
    df3_cum.loc[:, dsmost_seg] = df3.loc[:, toseg_idx].sum(axis=1)

# %%
# Check the segments for an individual outlet

dsmost_seg = 49294
pred = nx.dfs_predecessors(dag_us, dsmost_seg)
uniq_seg_us = uniq_seg_us.union(set(pred.keys()).union(set(pred.values())))

dag_ds_subset = dag_streamnet.subgraph(uniq_seg_us).copy()
node_outlets = [ee[0] for ee in dag_ds_subset.edges()]
true_outlets = set([dsmost_seg]).difference(set(node_outlets))

toseg_idx = list(set(xx[0] for xx in dag_ds_subset.edges))
toseg_idx.insert(0, dsmost_seg)
toseg_idx

# %%
# df1.loc[st:en, toseg_idx]  # .sum(axis=1)
df3.loc[st:en, 49294]

# %%

# %%
# %%time
st = datetime.datetime(2000, 1, 1)
en = datetime.datetime(2007, 12, 31)

df1_cum.head()

# %% [markdown]
# ## Write to csv files

# %%
df1_cum.loc[st:en, :].to_csv(f'/Users/pnorton/tmp/acc_{cvar}_byHRU_NHMv11.csv', sep=',', float_format='%0.4f')
df2_cum.loc[st:en, :].to_csv(f'/Users/pnorton/tmp/acc_{cvar}_byHW_NHMv11.csv', sep=',', float_format='%0.4f')
df3_cum.loc[st:en, :].to_csv(f'/Users/pnorton/tmp/acc_{cvar}_byHWobs_NHMv11.csv', sep=',', float_format='%0.4f')


# %%

# %%
# 49696
seg_chk = 49294

df1_cum.loc[:, seg_chk].plot(color='grey')

# byHW
df2_cum.loc[:, seg_chk].plot(color='goldenrod')

# byHW_obs
df3_cum.loc[:, seg_chk].plot(color='lightgreen')

# Compare the overall average for a segment

# byHRU
df1a = df1_cum.loc[st:en, :].resample('A').mean().mean()

# byHW
df2a = df2_cum.loc[st:en, :].resample('A').mean().mean()

# byHW_obs
df3a = df3_cum.loc[st:en, :].resample('A').mean().mean()

print(f'Overall average for {filename}, segment={seg_chk}')
print(f'   byHRU: {df1a.loc[seg_chk]}')
print(f'    byHW: {df2a.loc[seg_chk]}')
print(f' byHWobs: {df3a.loc[seg_chk]}')

# %%

# %%
# Compare the overall average for a segment

# byHRU
df1a = df1_cum.loc[st:en, :].resample('A').sum().sum()

# byHW
df2a = df2_cum.loc[st:en, :].resample('A').sum().sum()

# byHW_obs
df3a = df3_cum.loc[st:en, :].resample('A').sum().sum()

print(f'Overall total for {filename}, segment={seg_chk}')
print(f'   byHRU: {df1a.loc[seg_chk]}')
print(f'    byHW: {df2a.loc[seg_chk]}')
print(f' byHWobs: {df3a.loc[seg_chk]}')

# %%
# byHRU
df1.loc[st:en, seg_chk].plot(color='grey')

# byHW
df2.loc[st:en, seg_chk].plot(color='cyan')

# byHW_obs
df3.loc[st:en, seg_chk].plot(color='orange')


# %%

# %% [markdown]
# ## Annual mean

# %%

st = datetime.datetime(2000, 1, 1)
en = datetime.datetime(2007, 12, 31)

# byHRU
df1_cum.loc[st:en, seg_chk].resample('A').mean().plot(color='grey')

# byHW
df2_cum.loc[st:en, seg_chk].resample('A').mean().plot(color='goldenrod')

# byHW_obs
df3_cum.loc[st:en, seg_chk].resample('A').mean().plot(color='lightgreen')

# %%
st = datetime.datetime(2000, 1, 1)
en = datetime.datetime(2007, 12, 31)
df1_cum.loc[st:en, seg_chk].resample('A').mean().mean()

# %%
# byHRU
df1a = df1_cum.loc[st:en, :].resample('A').mean().mean()

# byHW
df2a = df2_cum.loc[st:en, :].resample('A').mean().mean()

# byHW_obs
df3a = df3_cum.loc[st:en, :].resample('A').mean().mean()

# %%
# %%time
poi = {}
if pfile.parameters.exists('poi_gage_segment'):
    poi_gage_id = pfile.parameters['poi_gage_id'].tolist()
    poi_gage_segment = pfile.parameters['poi_gage_segment'].tolist()
    
    for gg, ss in zip(poi_gage_id, poi_gage_segment):
        poi[nhm_seg[ss]] = gg

df1a = df1_cum.loc[st:en, :].resample('A').mean().mean()
df1a = df1a.loc[poi.keys()].to_frame()
df1a.reset_index(inplace=True)
df1a.columns = ['seg', 'byHRU']

# byHW
df2a = df2_cum.loc[st:en, :].resample('A').mean().mean()
df2a = df2a.loc[poi.keys()].to_frame()
df2a.reset_index(inplace=True)
df2a.columns = ['seg', 'byHW']

# byHW_obs
df3a = df3_cum.loc[st:en, :].resample('A').mean().mean()
df3a = df3a.loc[poi.keys()].to_frame()
df3a.reset_index(inplace=True)
df3a.columns = ['seg', 'byHWobs']

df_mrg = pd.merge(df1a, df2a, how='left', left_on='seg', right_on='seg')
df_mrg = pd.merge(df_mrg, df3a, how='left', left_on='seg', right_on='seg')
df_mrg.set_index('seg', inplace=True)

poi_id = pd.Series(poi)
poi_id.name = 'poi_id'

df_mrg = pd.merge(df_mrg, poi_id.to_frame(), how='left', left_index=True, right_index=True)
df_mrg.set_index('poi_id', inplace=True)

# %%
df_mrg

# %%

# %%

# %%

# %%

# %%

# %% jupyter={"outputs_hidden": true} tags=[]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(30, 20))
msize = 20
df1a.plot(ax=ax, kind='scatter', x='seg', y='gwflow', s=msize, c='grey', alpha=0.8)
df2a.plot(ax=ax, kind='scatter', x='seg', y='gwflow', s=msize, c='orange', alpha=0.8)
df3a.plot(ax=ax, kind='scatter', x='seg', y='gwflow', s=msize, c='lightgreen', alpha=0.8)

# %% jupyter={"outputs_hidden": true} tags=[]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(30, 20))
df1a.plot(ax=ax, kind='bar', x='seg', y='gwflow', color='grey', alpha=0.5)
df2a.plot(ax=ax, kind='bar', x='seg', y='gwflow', color='orange', alpha=0.5)
df3a.plot(ax=ax, kind='bar', x='seg', y='gwflow', color='lightgreen', alpha=0.5)

# %% jupyter={"outputs_hidden": true} tags=[]
df_mrg.plot(kind='bar', figsize=(30,12))

# %%

# %%
import xarray as xr
from dask.distributed import Client

client = Client()
client

# %%
client.close()

# %%
xdf = xr.open_mfdataset(f'{base_dir}/20220418_NHMv11_byHWobs_columbia_river/*.nc', parallel=True)
xdf

# %%
poi_gage_id = pfile.parameters['poi_gage_id'].tolist()

xdf_ss = xdf.loc[dict(poi_id=poi_gage_id)]

# %%
aa = xdf_ss.to_dataframe()

# %%
aa

# %%
xdf.loc[]

# %%
chk = xr.open_mfdataset(f'/Volumes/USGS_NHM2/NHM/NHM_v11/releases/gm_byHWobs_5.2.1/*.nc', chunks={})
chk = chk.assign_coords(nsegment=chk.nhm_seg)
chk

# %%
cc = chk['seginc_gwflow'].loc[st:en, toseg_idx]  # .to_pandas()

# %%
# cc.loc[dict(nsegment=42499)]
cc['nsegment']

# %%
# chk
chk['nhm_seg'].loc[toseg_idx]

# %%
cc.loc[dict(nsegment=49294)].plot()

# %%
cc.loc[dict(time=slice('2000-01-01', '2000-01-31'), nsegment=49294)].values

# %%

# %%
# self.__dataset = xr.open_dataset(self.__filename, decode_coords=True, chunks={'nsegment': 1000})
# self.__dataset = self.__dataset.assign_coords(nsegment=(self.__dataset.nhm_seg))
src_df = xr.open_dataset(f'/Volumes/USGS_NHM2/NHM/NHM_v11/releases/gm_byHWobs_5.2.1/seginc_gwflow.nc', decode_coords=True, chunks={'segment': 1000})
src_df = src_df.assign_coords(nsegment=src_df.nhm_seg)
src_df


# %%
# %%time

# data = self.__dataset[varname].loc[self.__stdate:self.__endate, self.__nhm_segs].to_pandas()
seg_list = [49294, 49295, 49296, 49297, 49298, 49299, 49300, 49301, 49302, 49303, 49304, 49305, 49306]

dst_df = src_df['seginc_gwflow'].loc[st:en, seg_list].to_pandas()
dst_df

# %%
# data.to_csv(f'{pathname}/{self.__varname}.csv', columns=self.__nhm_segs, sep=',', index=True, header=True, chunksize=50)
dst_df.to_csv('/Users/pnorton/tmp/acc_crap.csv', columns=seg_list, sep=',', index=True, header=True, chunksize=50)

# %%

# %%

# %%
