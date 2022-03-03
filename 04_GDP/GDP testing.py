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
#     display_name: Python [conda env:idp_gdp]
#     language: python
#     name: conda-env-idp_gdp-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import pyGDP
print('pyGDP version = ' + pyGDP.__version__)

# %%
import requests
s = requests.Session()
s.verify = '/Users/pnorton/files/certs/DOIRootCA2.pem'

print requests.certs.where()

# %%

sbURI = 'https://www.sciencebase.gov/catalogMaps/mapping/ows/54131d77e4b0239f1986bba7'
objGDP = pyGDP.pyGDPwebProcessing(wfs_url=sbURI)

# %%
pyGDP._execute_request.dodsReplace('https://internal.cida.usgs.gov/mows/thredds/dodsC/data/recharge_rates/recharge_2000-2013.nc')

# %%
output = objGDP.submitFeatureWeightedGridStatistics(
                        'sb:nhru_01', 
                        'https://internal.cida.usgs.gov/mows/thredds/dodsC/data/recharge_rates/recharge_2000-2013.nc', 'total_recharge', 
                        '2000-01-01T00:00:00.000Z', '2013-01-01T00:00:00.000Z',
                        attribute='hru_id_reg', value=None, gmlIDs=None, verbose=True,
                        coverage=False, delim='COMMA', weighted=True, stat='MEAN',
                        grpby='FEATURE_ATTRIBUTE')

# %%
with open('test.csv', 'w') as ff:
    ff.write(output)

# %%
print(output)

# %%
dataSetURIs = objGDP.getDataSetURI()

# %%
dataSetURIs

# %%
print(dataSetURIs)

# %%
for ii,jj in enumerate(dataSetURIs):
   if ii > 0:
       print jj[0]
       for kk, mm in enumerate(jj[2]):
           print '\t', kk, mm

# %%
dataSetURIs[1][2]

# %%
objGDP.getDataType('dods://cida.usgs.gov/thredds/dodsC/cmip5_bcca/historical')
