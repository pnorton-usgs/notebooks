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

# %%
import xml.etree.ElementTree as xmlET

# %%
param_root = xmlET.Element('parameter')
param_root.set('name', 'pname')
param_root.set('version', 'ver')

dims_sub = xmlET.Element('dimensions')
dim_sub = xmlET.SubElement(dims_sub, 'dimension')
dim_sub.set('name', 'dimname')
dim_sub.set('position', 'dimpos')
dim_sub.set('size', 'dimsize')

param_root.append(dims_sub)
xmlET.dump(param_root)
tree = xmlET.ElementTree(param_root)
tree.write('crap.xml', xml_declaration=True)

# %%
import xml.dom.minidom as minidom

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = xmlET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="\t")


# %%
prettify(param_root)

# %%
xmlstr = minidom.parseString(xmlET.tostring(param_root)).toprettyxml(indent='    ')
with open('crap.xml', 'w') as f:
    f.write(xmlstr.encode('utf-8'))

# %%
