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
#     display_name: Python [conda env:nco]
#     language: python
#     name: conda-env-nco-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
import datetime
# import fsspec
import os
# import tempfile
import xarray as xr

from nco import Nco

# %%
work_dir = '/Volumes/USGS_NHM2/datasets/wrf-conus404'
wrf_dir = '/Volumes/USGS_NHM2/datasets/wrf-conus404/originals'

# %%
# Create NCO object
nco = Nco(debug=True)

# %%
st_date = datetime.datetime(2019, 10, 1)
en_date = datetime.datetime(2019, 10, 2)

num_days = 1

delta = datetime.timedelta(days=num_days)   # Number of days

# %%
# %%time
c_start = st_date
fi = 1

# options for ncrcat
opts = ['-h',
        '-v Times,XTIME,T2',
        f'-p {wrf_dir}']  

while c_start < en_date:
    # print('-'*40)
    job_files = []
    
    for dd in range(num_days):
        cdate = c_start + datetime.timedelta(days=dd)
        
        wy_dir = f'WY{cdate.year}'
        if cdate >= datetime.datetime(cdate.year, 10, 1):
            wy_dir = f'WY{cdate.year+1}'

        for hh in range(24):
            fdate = cdate + datetime.timedelta(hours=hh)
            
            file_pat = f'{wy_dir}/wrf2d_d01_{fdate.strftime("%Y-%m-%d_%H:%M:%S")}'
            job_files.append(file_pat)
            
    # Concatenate files
    i_flist = ' '.join(job_files)
    nco.ncrcat(input=i_flist, output=f'{work_dir}/wrf2d_d01_1day.nc', options=opts)
    
    # Rechunk
    time_cnk = 1
    lat_cnk = 175
    lon_cnk = 175

    opts = ['-O', '-4', '-L 2', 
            '--cnk_map=dmn', 
            f'--cnk_dmn Time,{time_cnk}',
            f'--cnk_dmn south_north,{lat_cnk}', 
            f'--cnk_dmn west_east,{lon_cnk}']

    nco.ncks(input=f'{work_dir}/wrf2d_d01_1day.nc', output=f'{work_dir}/wrf2d_d01_1day_cnk.nc', options=opts)    
    
    c_start += delta

# %%

# %%
os.environ['PATH']

# %%

# %%

# %%

# %%
