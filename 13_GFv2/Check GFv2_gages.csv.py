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

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %% [markdown]
# Various checks for GFv2_gages.csv file from Andy Bock. 
#
# Version locations:
#     .../FY2021_projects/modeling_fabric/GFv2.0/20210812_POIs

# %%
import datetime
import geopandas as gpd
import pandas as pd
import re

# %%
work_dir = '/Users/pnorton/Projects/National_Hydrology_Model/FY2021_projects/modeling_fabric/GFv2.0/20210812_POIs'

file_prefix = 'GFv2_gages'
file_suffix = 'csv'
filename = f'{work_dir}/{file_prefix}.{file_suffix}'

out_filename = f'{work_dir}/{file_prefix}_cleaned.{file_suffix}'

# %%
df = pd.read_csv(filename, sep=',', quotechar='"', index_col=False)

# %%
df.info()

# %%

# %%

# %% [markdown]
# ## Read the data file

# %%
fhdl = open(filename, 'r')

rawdata = fhdl.read().splitlines()
fhdl.close()
it = iter(rawdata)

# next(it)  # Skip header
lines = []
first = True

for rec in it:
    if first:
        hdr_str = rec
        first = False
        continue
    lines.append(rec)

# %%
lines[0]

# %%
# Integer or "Integer".
# Real or "Real".
# String or "String".
# Date (format "YYYY-MM-DD"), Time (format "HH:MM:SS+nn") and DateTime (format "YYYY-MM-DD HH:MM:SS+nn"), whereas nn is the timezone.
# "WKT" (preferred over Point(X/Y)). All WKT geometry types are allowed: Point, LineString, Polygon, Multipoint, MultiLinestring, MultiPolygon, GeometryCollection, Arcs, ... (see OGC WKT).
# "CoordX","CoordY" (preferred) or "Point(X)","Point(Y)". Two separate colums in either order and not necessary neighboring of type Integer or Float: one containing the easting coordinate, and another containing northing coordinate separated by a comma.
# All values of that WKT column MAY contain the same geometry (sub)type.

# %% [markdown]
# ## Cleanup data and guess the datatypes for the fields
# Write new CSV file with cleaned data and also a CSVT file containing the datatypes

# %%
# known_fields = {'TotDASqKM_orig': 'Real'}

hdr_flds = hdr_str.replace('"', '').split(',')

# Test line needs a little cleanup for the geom field
tt = lines[0]
tt = re.sub(r'\([^()]+\)', lambda x: x.group().replace(',', ''), tt)

# Use the first line of data to guess the datatypes
a_str = tt.split(',')

fld_dtypes = []

for xx, yy in zip(a_str, hdr_flds):
    c_dtype = ''
    
    if yy == 'CLASS':
        c_dtype = 'String'
    elif yy in ['index', 'idx', '']:
        c_dtype = 'Integer'
    elif xx.startswith('"'):
        c_dtype = 'String'
    elif xx.startswith('c('):
        c_dtype = 'WKT'
    elif xx.endswith(')'):
        continue
    elif xx.isdecimal():
        c_dtype = 'Integer'
    else:
        try:
            float(xx)
            c_dtype = 'Real'
        except ValueError:
            try:
                d_lst = re.split('[-:+ ]', xx)
                d_int = [int(x) for x in d_lst]
                
                if len(d_lst) == 7:
                    datetime.datetime(*d_int)
                    c_dtype = 'DateTime'
                elif len(d_lst) == 4:
                    c_dtype = 'Time'
                elif len(d_lst) == 3:
                    if ':' in xx:
                        datetime.time(*d_int)
                        c_dtype = 'Time'
                    else:
                        datetime.date(*d_int)
                        c_dtype = 'Date'
                else:
                    c_dtype = 'ERROR'
            except ValueError:
                c_dtype = 'ERROR'
                
#     fld_dtypes.append(f'"{c_dtype}"')
    fld_dtypes.append(c_dtype)
    print(f'{yy}, {c_dtype}: {xx}')

# %%
fld_dtypes

# %%

# %% [markdown]
# ## Cleanup the data
# Clean the data and write to a new file

# %%
csvt_hdl = open(f'{out_filename}t', 'w')

csvt_hdl.write(f'{",".join(fld_dtypes)}\n')

csvt_hdl.close()

# %%
ohdl = open(out_filename, 'w')
hdr_str = ','.join(hdr_flds)

ohdl.write(f'{hdr_str}\n')

# lines_cln = []

for xx in lines:
    # Remove the comma that is between the parentheses
    # See: https://stackoverflow.com/questions/5099405
    xx = re.sub(r'\([^()]+\)', lambda x: x.group().replace(',', ''), xx)
    
    # Change the 'geom' field from c(.. ..) to POINT(.. ..)
    xx = xx.replace('c(', 'POINT(')
    
    # Remove the double-quotes
    xx = xx.replace('"', '')
    
    ohdl.write(f'{xx}\n')
    # lines_cln.append(xx)
    
ohdl.close()

# %%

# %%
df = pd.read_csv(out_filename, sep=',', index_col=False)

# %%
df.info()

# %%
df.head()

# %%
df['da_ratio'] = df['contrib_dr'] / df['drain_area']

# %%
df[df['da_ratio'] > 1.0]

# %%

# %%

# %%
