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
#     display_name: Python [conda env:bandit_py3]
#     language: python
#     name: conda-env-bandit_py3-py
# ---

# %%
from future.utils import iteritems

# %%
import xml.dom.minidom as minidom
import xml.etree.ElementTree as xmlET

# %%
output_vars = ['NHM-PRMS_climate', 'NHM-PRMS_flow', 'NHM-PRMS_snow', 'albedo', 'cap_waterin', 'contrib_fraction', 
               'dprst_area_open', 'dprst_evap_hru', 'dprst_insroff_hru', 'dprst_seep_hru', 'dprst_sroff_hru', 
               'dprst_stor_hru', 'dprst_vol_open', 'dprst_vol_open_frac', 'dunnian_flow', 'freeh2o', 'gw_in_soil', 
               'gw_in_ssr', 'gwres_flow', 'gwres_in', 'gwres_stor', 'hortonian_flow', 'hru_actet', 'hru_impervevap', 
               'hru_impervstor', 'hru_intcpevap', 'hru_intcpstor', 'hru_lateral_flow', 'hru_outflow', 'hru_ppt', 
               'hru_rain', 'hru_snow', 'hru_sroffi', 'hru_sroffp', 'hru_storage', 'hru_streamflow_out', 'infil', 
               'intcp_on', 'net_ppt', 'net_rain', 'net_snow', 'newsnow', 'perv_actet', 'pk_depth', 'pk_ice', 
               'pk_precip', 'pk_temp', 'pkwater_equiv', 'potet', 'pref_flow', 'pref_flow_in', 'pref_flow_infil', 
               'pref_flow_stor', 'prmx', 'recharge', 'slow_flow', 'slow_stor', 'snowcov_area', 'snow_evap', 
               'snow_free', 'snowmelt', 'soil_lower', 'soil_moist', 'soil_moist_tot', 'soil_rechr', 'soil_to_gw', 
               'soil_to_ssr', 'sroff', 'ssres_flow', 'ssres_in', 'ssres_stor', 'ssr_to_gw', 'swrad', 'tavgf', 
               'tmaxf', 'tminf', 'transp_on', 'unused_potet', 'seg_gwflow', 'seginc_gwflow', 'seginc_potet', 
               'seginc_sroff', 'seginc_ssflow', 'seginc_swrad', 'seg_inflow', 'seg_lateral_inflow', 
               'segment_delta_flow', 'seg_outflow', 'seg_sroff', 'seg_ssflow', 'seg_upstream_inflow', 'ssflow_save']


# %%
workfile = '/Users/pnorton/Projects/National_Hydrology_Model/src/fortran/prms6/src/xml/variables.xml'

xml_tree = xmlET.parse(workfile)
xml_root = xml_tree.getroot()

print(xml_root.tag)

# %%
var_dict = {}

# Iterate over child nodes of root
for elem in xml_root.findall('variable'):
    # print(elem.attrib.get('name'))
    name = elem.attrib.get('name')
    dtype = elem.find('type').text
    desc = elem.find('desc').text
    units = elem.find('units').text
    
    # Add dimensions for current parameter
    dims_list = []
    for cdim in elem.findall('./dimensions/dimension'):
        dims_list.append(cdim.attrib.get('name'))
    
    dims = ','.join(dims_list)
    # Print the elements of the node
    # print(elem.find('desc').text)
    # print(elem.find('size').text)
    # dim_size = int(elem.find('size').text)
    var_dict[name] = [dtype, desc, units, dims]
    
variables = list(var_dict.keys())
variables.sort()

# %%

# %%
missing = set(output_vars).difference(set(variables))
missing

# %%

# %%
datatype_map = {'I': 4, 'F': 5, 'D': 6, 'S': 2}
top_comment = '! Datatypes: NC_FLOAT=5, NC_DOUBLE=6, NC_INT=4, NC_CHAR=2'

# %%
var_dict

# %%
for vv, mt in iteritems(var_dict):
    print('Call: {}: {}, {}, {}, {}'.format(vv, mt[0], mt[1], mt[2], mt[3]))

# %%
