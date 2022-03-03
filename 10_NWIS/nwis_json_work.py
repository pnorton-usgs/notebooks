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
from __future__ import (absolute_import, division, print_function)
from future.utils import iteritems
# from tests.utils import resource_file
from owslib.waterml.wml11 import WaterML_1_1 as wml

# %%
import re
from urllib2 import urlopen, Request, HTTPError, URLError

# %%
BASE_NWIS_URL = 'http://waterservices.usgs.gov/nwis'

# Regex's for stripping unneeded clutter from the rdb file
t1 = re.compile('^#.*$\n?', re.MULTILINE)  # remove comment lines
t2 = re.compile('^5s.*$\n?', re.MULTILINE)  # remove field length lines

# %%
# https://waterservices.usgs.gov/nwis/dv/?format=waterml,2.0&sites=01465500&statCd=00003&parameterCd=00060&startDT=2014-9-25&endDT=2014-10-4&siteStatus=all

req_html = '{}/dv/?format=json&sites=01465500&statCd=00003&parameterCd=00060&startDT=2014-9-25&endDT=2014-10-4&siteStatus=all'.format(BASE_NWIS_URL)
print(req_html)

# %%
# Open the URL
streamgage_obs_page = urlopen(req_html)

# %%
streamgage_observations = streamgage_obs_page.read()

# %%
print(streamgage_observations)

# %%
import json
import pprint


# %%
def traverse(obj, path=None, callback=None):
    """
    Traverse an arbitrary Python object structure (limited to JSON data
    types), calling a callback function for every element in the structure,
    and inserting the return value of the callback as the new value.
    """
    if path is None:
        path = []
     
    if isinstance(obj, dict):
        value = {k: traverse(v, path + [k], callback) for k, v in obj.items()}
    elif isinstance(obj, list):
        value = [traverse(elem, path + [[]], callback) for elem in obj]
    else:
        value = obj
    
    if value is None:
        return '-----------------'
    else:
        return value

    #if callback is None:
    #    return value
    #else:
    #    return callback(path, value) 

def removeNulls2(obj):
    # Recursively traverse the JSON object and remove all dictionary entries 
    # where the value is None or contains an empty list.
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, dict):
                removeNulls2(v)
            elif isinstance(v, list):
                if len(v) == 0:
                    # Empty list
                    del obj[k]
                else:
                    removeNulls2(v)
            else:
                if v is None:
                    del obj[k]
    elif isinstance(obj, list):
        for elem in obj:
            if isinstance(elem, dict):
                removeNulls2(elem)

def convertDateTime(obj):
    # Convert dateTime entries from YYYY-MM-DDTHH:mm:ss-TZ to YYYY-MM-DD
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, dict):
                convertDateTime(v)
            elif isinstance(v, list):
                convertDateTime(v)
            else:
                if k == u'dateTime':
                    if len(v) > 10:
                        # Chop off everything but the date
                        obj[k] = v[:10]
    elif isinstance(obj, list):
        for elem in obj:
            if isinstance(elem, dict):
                convertDateTime(elem)
    

def changeKeys(obj, origKey, newKey):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, dict):
                changeKeys(v, origKey, newKey)
            elif isinstance(v, list):
                changeKeys(v, origKey, newKey)
            #else:
            if k == origKey:
                obj[newKey] = obj[k]
                del obj[k]
    elif isinstance(obj, list):
        for elem in obj:
            if isinstance(elem, dict):
                changeKeys(elem, origKey, newKey)



def deleteKeys(obj, delKey):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, dict):
                deleteKeys(v, delKey)
            elif isinstance(v, list):
                deleteKeys(v, delKey)
            #else:
            if k == delKey:
                #obj[newKey] = obj[k]
                del obj[k]
    elif isinstance(obj, list):
        for elem in obj:
            if isinstance(elem, dict):
                deleteKeys(elem, delKey)


def delMissing(obj, missingVal):
    # Recursively search the JSON data for missingVal. If the missing value
    # occurs within a list of dictionary items then remove the parent list
    # element (which removes the all the contained dictionary items).
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, dict):
                delMissing(v, missingVal)
            elif isinstance(v, list):
                delMissing(v, missingVal)
            elif v == missingVal:
                # Flag when the value is equal to the missingVal
                return True
                #break
        return False
    elif isinstance(obj, list):
        for elem in obj[:]:
            if isinstance(elem, dict):
                if delMissing(elem, missingVal):
                    # Remove the list element if there is a missing value 
                    # for at least one of the contained dictionary keys.
                    obj.remove(elem)


# %%
tmp = json.loads(streamgage_observations)

# %%
type(tmp)

# %%
for xx, yy in iteritems(tmp):
    print(xx)
    
    if isinstance(yy, dict):
        for zz in yy:
            print('  ', zz)


# %%
delMissing(tmp, u'-999999')
removeNulls2(tmp)
convertDateTime(tmp)
pprint.pprint(tmp)

# %%
for xx, yy in enumerate(tmp['value']['timeSeries'][0]['values']):
    for tt, zz in enumerate(yy['value']):
        print(zz['dateTime'], zz['value'])

# for xx, yy in enumerate(tmp['value']['timeSeries'][0]['values'][0]['value']):
#     print(xx, yy)

# %%
tmp['value']['timeSeries'][0]['values'][0]['value']

# %%
