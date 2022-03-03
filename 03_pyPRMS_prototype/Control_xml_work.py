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

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
# from future.utils import iteritems

import xml.dom.minidom as minidom

try:
    import xml.etree.cElementTree as xmlET
except ImportError:
    import xml.etree.ElementTree as xmlET


# %%
# workfile = '/Users/pnorton/Projects/National_Hydrology_Model/src/fortran/prms6/src/xml/control.xml'
workfile = '/Users/pnorton/PycharmProjects/pyPRMS/pyPRMS/xml/control.xml'

xml_tree = xmlET.parse(workfile)
xml_root = xml_tree.getroot()

# print(xml_root.tag)

control = {}
# Iterate over child nodes of root
for elem in xml_root.findall('control_param'):
    # print(elem.attrib.get('name'))
    name = elem.attrib.get('name')
    var_version_major = version_info(elem.attrib.get('version'))[0]
    
    # Print the elements of the node
    # print(elem.find('desc').text)
    # print(elem.find('size').text)
    ctl_type = int(elem.find('type').text)

    
#                     # Add dimensions for current parameter
#                 for cdim in elem.findall('./dimensions/dimension'):
#                     self.parameters.get(name).dimensions.add(cdim.attrib.get('name'))

#                 for cmod in elem.findall('./modules/module'):
#                     self.parameters.get(name).modules = cmod.text
    
    # Lookup what the control variable values mean
    for cvals in elem.findall('./values'):
#         print(cval.attrib.get('name'))
        if var_version_major == 6:
            val_type = cvals.attrib.get('type')
            print(name, var_version_major)

            print('type: {}'.format(val_type))

            for cv in cvals.findall('./value'):
                outvals = []
                for xx in cv.text.split(','):
                    outvals.append(xx)

                print('{}: {}'.format(cv.attrib.get('name'), outvals))

    control[name] = ctl_type
    
ctl_names = list(control)
ctl_names.sort()

# %%
control.keys()


# %%
def version_info(version_str, delim='.'):
    
    if version_str is None:
        return [0, 0, 0]
    
    flds = [int(kk) for kk in version_str.split(delim)]
    return flds
    
    

# %%
version_info('5.1.0')

# %%
