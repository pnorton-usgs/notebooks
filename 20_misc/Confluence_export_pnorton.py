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
#     display_name: Python [conda env:.conda-bandit_nhgf]
#     language: python
#     name: conda-env-.conda-bandit_nhgf-py
# ---

# %%
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from urllib.parse import unquote

import argparse
import errno
import getpass
import os
import os.path
import re
import requests

# %%
# Set the starting pageid
pageid = '84443304'  # MOWS
pageid = '211615823'

# %%
base_url = 'https://my.usgs.gov/confluence'

session = requests.Session()
retry = Retry(connect=10, backoff_factor=1)
adapter = HTTPAdapter(max_retries=retry)
session.mount('http://', adapter)
session.mount('https://', adapter)

def fix_names(tofix):
    # Characters that OneDrive and SharePoint don't like
    # /\:*?"<>|
    u_tofix = unquote(tofix)
    u_tofix = re.sub('[ \\/\*\?\"\<\>\-,()|]', '_', u_tofix)
    u_tofix = re.sub('[:]', '-', u_tofix)

    return u_tofix.replace('___', '_').replace('__', '_').rstrip('_')


def download(url, filename, user, passwd):
    response = session.get(url, auth=(user, passwd), stream=True)
    if response.status_code == 200:
        if response.headers.get('Content-Disposition'):
            if response.headers.get('Content-Length') != '118':
                open(filename, 'wb').write(response.content)
            else:
                print('Empty zipfile encountered: ' + filename)


def build_path(pageid, user, passwd):
    full_path = ""
    ancestor_url = f'{base_url}/rest/api/content/{pageid}?expand=ancestors'
    response = session.get(ancestor_url, auth=(user, passwd))
    json_resp = response.json()

    for item in json_resp["ancestors"]:
        # Use OS path separators
        full_path += fix_names(item["title"]) + os.path.sep
    if not full_path:
        full_path = fix_names(json_resp["title"]) + os.path.sep

    return full_path


def process_page(pageid, page_name, user, passwd):
    clean_page_name = fix_names(page_name).replace('/', '_')
    page_path = build_path(pageid, user, passwd).replace('.', '')
    print(f'*** {clean_page_name=}')

    if not os.path.exists(page_path):
        try:
            os.makedirs(page_path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    # Export a page as a Word doc:
    word_url = f'{base_url}/exportword?pageId={pageid}'
    output_filename = f'{page_path}{clean_page_name}.doc'
    download(word_url, output_filename, user, passwd)

    # Download page attachments:
    list_attachments_url = f'{base_url}/rest/api/content/{pageid}/child/attachment?limit=1000'
    response = session.get(list_attachments_url, auth=(user, passwd))
    json_resp = response.json()

    for item in json_resp['results']:
        download_url = f'{base_url}{item["_links"]["download"]}'
        output_path = f'{page_path}/{clean_page_name}_attachments'

        if not os.path.exists(output_path):
            try:
                os.makedirs(output_path)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise

        # print('-' * 40)
        # print(f'** {item["title"]=}')
        # print(f'{output_path=}')
        print(f'Downloading: {download_url}')
        download(download_url, f'{output_path}/{fix_names(item["title"])}', user, passwd)

    # This will return all children, need to call recursivly:
    child_url = f'{base_url}/rest/api/content/search?limit=750&cql=parent={pageid}'
    response = session.get(child_url, auth=(user, passwd))
    json_resp = response.json()

    for item in json_resp["results"]:
        # print(f'Processing item id: {item["id"]}')
        build_path(pageid, user, passwd)
        process_page(item["id"], item["title"], user, passwd)


# %%
# Need to build from user input
user = input('AD Username (FULL EMAIL ADDRESS): ')
passwd = getpass.getpass(prompt='Password: ')

# %%
init_page_url = f'{base_url}/rest/api/content/{pageid}'
response = requests.get(init_page_url, auth=(user, passwd))

json_resp = response.json()

process_page(pageid, json_resp['title'], user, passwd)

# %%
