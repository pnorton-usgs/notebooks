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

# %%

# %%
import os

# %%
Tcratch notebook to figure out why th

# %%
work_dir = '/Users/pnorton/OneDrive - DOI/MoWS_Home'

# %%
fname = []

for root, dir_names, file_names in os.walk(work_dir):
    for cfile in file_names:
        fname.append(os.path.join(root, cfile))

# %%
largest = 0

for ff in fname:
    largest = max(largest, len(ff))
    if len(ff) >= 400:
        print(ff)

# %%
# Illegal characters
# " * : < > ? / \ |
bad_chars = '"*:<>?|\\'
new_filenames = []

def containsAny(str, set):
    """ Check whether sequence str contains ANY of the items in set. """
    return 1 in [c in str for c in set]


for ff in fname:
    if containsAny(ff, bad_chars):
        new_filenames.append([ff, ff.replace(':', '-')])
        print(ff)

# %%
for fpair in new_filenames:
    os.rename(fpair[0], fpair[1])

# %%

# %%

# %%
