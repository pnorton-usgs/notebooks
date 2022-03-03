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
#     display_name: Python 2
#     language: python
#     name: python2
# ---

# %%

# %matplotlib inline
import matplotlib as mpl
#from matplotlib.backends.backend_pdf import PdfPages
#mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import pylab
from mpl_toolkits.basemap import Basemap

# %%
stn_info = pd.read_csv('gageNS_xy', r"\s*", header=None, index_col=['site_no'],
                       names=['site_no', 'stat1', 'stat2', 'lat', 'lon'], engine='python')
stn_info.drop(['stat1', 'stat2'], axis=1, inplace=True)

stn_info.head()

# %%
stn_ns = pd.read_csv('luca_ns.csv', index_col=['site_no'])
stn_ns.head()

# %%
stn = pd.merge(stn_info, stn_ns, left_index=True, right_index=True, how='outer')
print stn.shape
stn.dropna(axis=0, how='any', thresh=None, subset=['ns'], inplace=True)
print stn.shape
stn['color'] = ''
stn['color'] = np.where(stn['ns'] > 0.8, '#cc0000', stn['color'])
stn['color'] = np.where(np.logical_and(stn['ns'] > 0.6, stn['ns'] <= 0.8), '#cc33cc', stn['color'])
stn['color'] = np.where(np.logical_and(stn['ns'] > 0.4, stn['ns'] <= 0.6), '#6666ff', stn['color'])
stn['color'] = np.where(np.logical_and(stn['ns'] > 0.2, stn['ns'] <= 0.4), '#339933', stn['color'])
stn['color'] = np.where(np.logical_and(stn['ns'] >= 0.0, stn['ns'] <= 0.2), '#cccc33', stn['color'])
stn['color'] = np.where(stn['ns'] < 0.0, '#333333', stn['color'])
print stn.head()

# %%
# Setup output to a pdf file
#outpdf = PdfPages('NSE_map.pdf')

extent=(50, 22, -64, -119) #Extension for USA

fig, ax = pylab.subplots(nrows=1, ncols=1, figsize=(17,11))
ax = pylab.gca()

north, south, east, west = extent
m = Basemap(llcrnrlon=west, llcrnrlat=south, urcrnrlon=east, urcrnrlat=north,
            projection='aea', resolution='l', lat_0=(south+north)/2, lon_0=(east+west)/2)

# draw parallels.
#parallels = np.arange(0.,90,10.)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=20)

# draw meridians
#meridians = np.arange(180.,360.,10.)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=20)
#m.drawmapboundary()

m.drawcoastlines()
m.drawcountries()
m.drawstates(color='#666666')
m.drawmapboundary(fill_color='#eeeeee')
m.fillcontinents(color='#cc9966',lake_color='#0066cc', zorder=0)

#lats = stn['lat'].tolist()
#lons = stn['lon'].tolist()
#x, y = m(lons,lats)
cols = ['#cc0000', '#cc33cc', '#6666ff', '#339933', '#cccc33', '#333333']
labels = ["> 0.8", "0.6 to 0.8", "0.4 to 0.6", "0.2 to 0.4", "0.0 to 0.2", "< 0.0"]
pl = []

for pp in cols:
    ss = stn[stn['color'] == pp]
    
    lats = ss['lat'].tolist()
    lons = ss['lon'].tolist()
    x, y = m(lons,lats)
    
    pl.append(m.scatter(x, y, 12, marker='o', color=ss['color']))
#plt.title('Locations of %s ARGO floats active between %s and %s' %\
#        (len(lats),date1,date2),fontsize=12)

leg = plt.legend(pl, labels, ncol=1, frameon=True, fontsize=12,
                 handlelength=None, loc = 'lower left', borderpad = 1.8,
                 handletextpad=None, title='NSE', scatterpoints = None)
frame = leg.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')

fig.savefig('NSE_map.png', dpi=250, bbox_inches='tight')
#outpdf.savefig()
#outpdf.close()
#plt.show()



# %%
