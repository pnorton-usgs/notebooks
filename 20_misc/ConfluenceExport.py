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
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% code_folding=[0]
# #!/usr/bin/python
import re
import sys
import os
import os.path
import getpass
from urllib.request import urlretrieve
from urllib.parse import unquote
import requests
import errno
import json
import zipfile
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

session = requests.Session()
retry = Retry(connect=10, backoff_factor=1)
adapter = HTTPAdapter(max_retries=retry)
session.mount('http://', adapter)
session.mount('https://', adapter)

pageid=sys.argv[1]

# %% code_folding=[0]
# Need to build from user input
user=input('AD Username (FULL EMAIL ADDRESS): ')
passwd = getpass.getpass(prompt='Password: ')

# %% code_folding=[0]
# YOU MUST EDIT THIS CELL TO SPECIFY WHICH PARENT PAGE TO START EXPORTING
#Right click 'edit' on the top-level Confluence page you'd like to edit, and paste in the page ID URL below. 

# EXAMPLE: 
# pageid = '148963437' # Monthly mtgs 2011

# pageid ='460324995' # REPLACE NUMBER AND MAKE COMMENT about what you're exporting
pageid = '84443304'  # MOWS

# pageid = '155025420'  # near Kersey

# pageid = '211615823'  # near Estes

# pageid = '544049182' # snow depletion curve images

# pageid = '86278233'
# pageid = '570080806'   # Subset visualizations

# %% code_folding=[0]
base_url = 'https://my.usgs.gov/confluence/confluence'

def fix_names(tofix):    
    # 'MoWS - Applications - NHDPlus Region 10 Lower - Highest Res HRUs above Big Thompson River near Estes Park (gage 06735500)'
    # Characters that OneDrive and SharePoint don't like
    # /\:*?"<>|
    u_tofix = unquote(tofix)
    u_tofix = re.sub('[ \\/\*\?\"\<\>\-,()|]', '_', u_tofix)
    u_tofix = re.sub('[:]', '-', u_tofix)
    
    return u_tofix.replace('___', '_').replace('__', '_').rstrip('_')
    # return ''.join([x if x.isalnum() else '_' for x in tofix]).replace('___', '_').replace('__', '_').rstrip('_')
    
def download(url, filename):
    response = session.get(url, auth=(user, passwd), stream=True)
    if response.status_code == 200:
        if response.headers.get('Content-Disposition'):
            #print(response.headers.get('Content-Length'))
            #check content-length != 118 (empty zip)
            if response.headers.get('Content-Length') != '118':
                open(filename, 'wb').write(response.content)
            else:
                print('Empty zipfile encountered: ' + filename)

def buildPath(pageid):
    # fullPath="./"
    fullPath = ""
    # https://my.usgs.gov/confluence/rest/api/content/544051479?expand=ancestors
    # ancestorUrl = 'https://my.usgs.gov/confluence/rest/api/content/' + pageid + '?expand=ancestors'
    ancestorUrl = f'{base_url}/rest/api/content/{pageid}?expand=ancestors'
    response = session.get(ancestorUrl, auth=(user, passwd))
    jsonResp = response.json()
    
    for item in jsonResp["ancestors"]:
        # Use OS path separators
        # fullPath+="".join([x if x.isalnum() else "_" for x in item["title"]])+os.path.sep
        fullPath += fix_names(item["title"]) + os.path.sep
    if not fullPath:
        # fullPath=fullPath="".join([x if x.isalnum() else "_" for x in jsonResp["title"]])+os.path.sep
        fullPath = fix_names(jsonResp["title"]) + os.path.sep
                    
    return fullPath

def processPage(pageid, pageName):
    cleanPageName = fix_names(pageName).replace('/', '_')
    print(f'*** cleanPageName: {cleanPageName}')
    pagePath = buildPath(pageid).replace('.', '')

    if not os.path.exists(pagePath):
        try:
            os.makedirs(pagePath)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    # Export a page as a Word doc:
    # wordUrl = 'https://my.usgs.gov/confluence/exportword?pageId=' + pageid
    wordUrl = f'{base_url}/exportword?pageId={pageid}'
    output_filename = f'{pagePath}{cleanPageName}.doc'
    download(wordUrl, output_filename)

    # Download page attachments:
#     listAttachmentsUrl = 'https://my.usgs.gov/confluence/rest/api/content/' + pageid + '/child/attachment?limit=1000'
    listAttachmentsUrl = f'{base_url}/rest/api/content/{pageid}/child/attachment?limit=1000'
    response = session.get(listAttachmentsUrl, auth=(user, passwd))
    jsonResp = response.json()
    
    for item in jsonResp['results']:
        # downloadUrl = 'https://my.usgs.gov/confluence' + item['_links']['download']
        downloadUrl = f'{base_url}{item["_links"]["download"]}'
        # output_path = pagePath + '/' + cleanPageName + '_attachments/'
        output_path = f'{pagePath}/{cleanPageName}_attachments'
                    
        if not os.path.exists(output_path):
          try:
            os.makedirs(output_path)
          except OSError as e:
            if e.errno != errno.EEXIST:
              raise
                    
        print('-'*40)
        print(f'** Original name: {item["title"]};  path: {output_path}')
        print(f'Attempting to download file: {downloadUrl}')
        # download(downloadUrl, output_path + fix_names(item['title'])
        download(downloadUrl, f'{output_path}/{fix_names(item["title"])}')
    
    # This will return all children, need to call recursivly:
    # childrenUrl = 'https://my.usgs.gov/confluence/rest/api/content/search?limit=750&cql=parent=' + pageid
    childrenUrl = f'{base_url}/rest/api/content/search?limit=750&cql=parent={pageid}'
    response = session.get(childrenUrl, auth=(user, passwd))
    jsonResp = response.json()
    
    for item in jsonResp["results"]:
      print(f'Processing item id: {item["id"]}')
      buildPath(pageid)
      processPage(item["id"], item["title"])
    #https://my.usgs.gov/confluence/rest/api/content/489357375/child/page

initPageUrl = f'{base_url}/rest/api/content/{pageid}'
# response = requests.get(initPageUrl, auth=(user, passwd))
# jsonResp=response.json()
# processPage(pageid,"".join([x if x.isalnum() else "_" for x in jsonResp["title"]]))


# %%
response = requests.get(initPageUrl, auth=(user, passwd))
jsonResp = response.json()

# %%
processPage(pageid,"".join([x if x.isalnum() else "_" for x in jsonResp["title"]]))

# %%
fix_names(jsonResp["title"])

# %%
"".join([x if x.isalnum() else "_" for x in jsonResp["title"]]).replace('___', '_').replace('__', "_") + os.path.sep

# %%
aa = 'MoWS_Home/MoWs_Current_project_management_information/MoWS_Model_application_plans/MoWS_Applications_NHDPlus_Region_10_Lower_South_Platte_River_near_Kersey_CO._subset_of_10c//MoWS_Applications_NHDPlus_Region_10_Lower_Highest_Res_HRUs_above_Big_Thompson_River_near_Estes_Park_gage_06735500_attachments'

# %%
os.path.basename(aa)

# %%
fix_names('image2016-9-21%2015%3A12%3A12.png')

# %%

# %%
