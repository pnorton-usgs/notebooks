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
#     display_name: Python [conda env:holoview]
#     language: python
#     name: conda-env-holoview-py
# ---

# %%

# %%
import holoviews as hv
import hvplot.pandas
import numpy as np
import pandas as pd

from bokeh.models.formatters import NumeralTickFormatter
from holoviews.element.tiles import OSM
from holoviews.util.transform import lon_lat_to_easting_northing

# %%
wateryears = True

workdir = '/Users/pnorton/Program Development/FY2021_ReFED/High_Flows_SIR/data/NWIS/datapull_CONUS_20210205/trends_wy1960-2019'

filename = f'{workdir}/conus_annual_wy1960-2019_kendall.csv'

# %% [markdown]
# Streamflow records for water years 1960–2019 for 1,853 USGS streamgages with continuous record of annual streamflow were obtained from the USGS National Water Information System (NWIS) database. The streamgages were examined for significant changes in streamflow using the Kendall Tau statistical test (p-value less than or equal to 0.10). A total of 573 sites had a statistically significant upward trend and 182 sites had a statistically significant downward trend.
#
# This file contains streamgage information, the results of the Kendall's Tau test, the mean annual streamflow for first and last ten years of record, and the percentage change in streamflow comparing the first and last ten year means. 

# %%
# site_no	station_nm	dec_lat_va	dec_long_va	drain_area_va	contrib_drain_area_va	pval	tau	trend	first_ten_yr	last_ten_yr	pct_chg
stn_col_names = ['site_no', 'station_nm', 'dec_lat_va', 'dec_long_va',
                 'drain_area_va', 'contrib_drain_area_va',
                 'pval', 'tau', 'trend', 'first_ten_yr', 'last_ten_yr', 'pct_chg']
stn_col_types = [np.str_, np.str_, float, float, float, float,
                 float, float, int, float, float, float]
stn_cols = dict(zip(stn_col_names, stn_col_types))

stations = pd.read_csv(filename, sep='\t', usecols=stn_col_names,
                       dtype=stn_cols)

stations.set_index('site_no', inplace=True)

# Have to force numeric conversion after the fact when there are
# null values in a column
for dd in stations.columns:
    if stn_cols[dd] == np.float_:
        stations[dd] = pd.to_numeric(stations[dd], errors='coerce')
        
stations.loc[:, 'long_merc'], stations.loc[:, 'lat_merc'] = lon_lat_to_easting_northing(stations.dec_long_va,stations.dec_lat_va)

stations.pct_chg *= 100.0

# stations.head()

# %%

# %% [markdown]
# Restrict the dataset to only those sites with a statistically significant (p-value <= 0.10) trend in annual streamflow (trend = ±2 are non-significant trends, trend = ±1 are statistically significant trends.

# %%
df_trends = stations[stations['trend'] > -2]
df_trends = df_trends[df_trends['trend'] < 2]

# df_trends.head()

# %%

# %%
df_trends = stations[stations['trend'] > -2]
df_trends = df_trends[df_trends['trend'] < 2]

da_bins = [-np.inf, 50, 100, 500, 1000, 5000, 10000, 50000, np.inf]
da_names = ['Less than 50', '50 to 100', '100 to 500', '500 to 1,000', '1,000 to 5,000', '5,000 to 10,000', 
            '10,000 to 50,000', 'Greater than 50,000']

da_class = pd.cut(x=df_trends['drain_area_va'], bins=da_bins, labels=da_names)
df_trends.insert(1, 'da_class', da_class)
# df_trends['da_class_counts'] = df_trends.da_class.value_counts()

formatter = NumeralTickFormatter(format="0,0")

# crime.hvplot.bar(x='Year', y=['Violent crime total', 'Property crime total'],
#                  stacked=True, rot=90, width=800, legend='top_left')
sig_da = df_trends.hvplot(y='drain_area_va', by='da_class', kind='hist', alpha=0.5, responsive=True, min_height=300)
sig_da.opts(xformatter=formatter)

sig_pchg = df_trends.hvplot(y='pct_chg', kind='hist', bin_range=(-100, 500), bins=12, responsive=True, min_height=200)
# sig_pchg = df_trends.hvplot(y='pct_chg', kind='hist', responsive=True, min_height=200)
# sig_pchg

link_plots = hv.link_selections.instance()

link_plots(sig_da + sig_pchg)

# %%
da_by_class = df_trends.da_class.value_counts().reset_index().to_numpy()
hv.Bars(da_by_class, hv.Dimension('Drainage area class'), 'Count').opts(xrotation=60)

# %%
# data = [('one',8),('two', 10), ('three', 16), ('four', 8), ('five', 4), ('six', 1)]

# bars = hv.Bars(data, hv.Dimension('Car occupants'), 'Count')

# bars

# %%
df_trends.da_class.value_counts().reset_index().to_numpy()

# %%

# %%

# %%
plot_hv = df_trends.hvplot(x='long_merc', y='lat_merc', c='pct_chg', kind='points', cmap='bwr_r', responsive=True, height=350)

plot_hv.opts(clim=(-500, 500))

pchg_map = OSM() * plot_hv

linked = hv.link_selections.instance()

(linked(pchg_map + sig_pchg + sig_da)).cols(1)

# %%

# %%

# %%
