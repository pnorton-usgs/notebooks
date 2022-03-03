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

# %%
import xml.dom.minidom as minidom
import xml.etree.ElementTree as xmlET
import numpy as np
from collections import namedtuple

# %%
workdir = "/Users/pnorton/PycharmProjects/xml_scratch"

# %%
dims_tree = xmlET.parse('{}/dimensions.xml'.format(workdir))
dims_root = dims_tree.getroot()

dims = []
for dimension in dims_root.findall('dimension'):
    dims.append(dimension.get('name'))
    
print(dims)

# %%
dims_xml = xmlET.Element('dimensions')

# for kk, vv in iteritems(self.dimensions):
#     dim_sub = xmlET.SubElement(dims_xml, 'dimension')
#     dim_sub.set('name', kk)
#     xmlET.SubElement(dim_sub, 'position').text = str(self.get_position(kk)+1)
#     xmlET.SubElement(dim_sub, 'size').text = str(vv.size)

# %%
dims_root.tag

# %% [markdown]
# ### Output variables

# %%
outvars_tree = xmlET.parse('{}/variables.xml'.format(workdir))
outvars_root = outvars_tree.getroot()

print(outvars_root.tag)

# %%
var_dict = {}

# Iterate over child nodes of root
for elem in outvars_root.findall('variable'):
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
var_dict


# %%
def get_variable(xml_root, var_name):
    OutputVariable = namedtuple('OutputVariable', 'name, dimension, type, units, desc')
    
    for elem in xml_root.findall('variable'):
        # print(elem.attrib.get('name'))
        name = elem.attrib.get('name')
        dtype = elem.find('type').text

        if var_name == name:
            for cdim in elem.findall('.dimensions/dimension'):
                dim = cdim.attrib.get('name')

            desc = elem.find('desc').text
            units = elem.find('units').text
            dimname = dim[1:]

            return OutputVariable(name, dimname, dtype, units, desc)


# %%
get_variable(outvars_root, 'albedo')

# %%

# %%
xmlstr = minidom.parseString(xmlET.tostring(outvars_root)).toprettyxml(indent='    ')

# %%
print(xmlstr)

# %%
