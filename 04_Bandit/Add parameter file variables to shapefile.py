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

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
from future.utils import iteritems

import pyPRMS.ParameterFile as pf
reload(pf)

from osgeo import ogr
import gdal


# %%
workdir='/Users/pnorton/Projects/National_Hydrology_Model/datasets/bandit/jobs/20180919_florida/02450000'
filename = '{}/myparam.param'.format(workdir)

pfile = pf.ParameterFile(filename)

# %%
myparams = pfile.parameters
mydims = pfile.dimensions

print(mydims)

# %%
myparams.check()

# %%
# Loop through the parameter objects in the ParameterSet and return parameters which have the nhru dimension
for vv in myparams.values():
    if 'nhru' in vv.dimensions.keys() and vv.dimensions.ndims == 1:
        print('{}: '.format(vv.name), vv.dimensions.keys())
    
#     for dd in vv.dimensions.values():
#         print('\t{} = {}'.format(dd.name, dd.size))

# %%
# Add parameters to existing shapefile
filename_gis = '{}/GIS/HRU_subset.shp'.format(workdir)

# Load the file - this assumes a file geodatabase
# driver = ogr.GetDriverByName('OpenFileGDB')
# gdb = driver.Open(filename_gis)

# open the shapefile
driver = ogr.GetDriverByName('ESRI Shapefile')
dataSource = driver.Open(filename_gis, 1) # open for rw

if dataSource is None:
    print("ERROR: could not open {} as shapefile!".format(filename_gis))
#     sys.exit(1)

layer = dataSource.GetLayer()
# layer.CreateField(ogr.FieldDefn("area",     ogr.OFTReal))

for vv in myparams.values():
    if 'nhru' in vv.dimensions.keys() and vv.dimensions.ndims == 1 and vv.name != 'nhm_id':
#         print('{}: '.format(vv.name), vv.dimensions.keys())
        layer.CreateField(ogr.FieldDefn(vv.name, ogr.OFTReal))
        
for feature in layer:
    cid = feature.GetField('hru_id_nat')
    print(cid)
    
    for vv in myparams.values():
        if 'nhru' in vv.dimensions.keys() and vv.dimensions.ndims == 1 and vv.name != 'nhm_id':
#             print(vv.name)
            # Set the value for the feature
            cvals = myparams.get_dataframe(vv.name).loc[cid].values[0]
            feature.SetField(vv.name, cvals)


            layer.SetFeature(feature)
    
dataSource = None





# # Another example
# shapefile = ogr.Open(shapefile_path, 1)

# layer = shapefile.GetLayer()
# layer_defn = layer.GetLayerDefn()

# new_field_defn = ogr.FieldDefn("New_field", ogr.OFTReal)
# new_field_defn.SetWidth(50)
# new_field_defn.SetPrecision(11)
# layer.CreateField(new_field_defn)

# # Walk through shapefile, setting new field for each feature
# for feature in layer :
#     geometry = feature.GetGeometryRef()             

#     band_value = 100
#     feature.SetField("New_field",np.double(band_value))

#     # This is the line that crashes the program
#     layer.SetFeature(feature)

# shapefile.Destroy()

# %%
myparams.get_dataframe('va_open_exp').loc[cid].values

# %%
myparams.get_dataframe('nhm_id')

# %%
aa = [8996,8997,8998,9003,9017,9020,9022,9023,9025,9051,9052,9056,9070,9076,9087,9088,9090,9102,9115,9472,9473,9474,9475,9476]

bb = myparams.get_dataframe('va_open_exp').index.values.tolist()


# %%
set(aa).difference(set(bb))

# %%
out_driver = ogr.GetDriverByName('GPKG')
out_dataSource = out_driver.CreateDataSource('tmp.gpkg')

# %%
print(out_dataSource)

# %%
