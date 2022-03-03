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
# import xarray as xr
import fsspec
import os
import tempfile
import xarray as xr

from nco import Nco
from nco.custom import Atted

# %%
# Sample gridmet filename: 2020_gm_tmin_2021_03_31.nc
#     year ----------------^^^^
#     variable --------------------^^^^
#     date retrieved -------------------^^^^^^^^^^

var = 'tmin'
src_dir = '/Volumes/USGS_NHM2/datasets/gridmet/gridmet_raw'
output_dir = '/Volumes/USGS_NHM2/datasets/gridmet/test_output'
output_filename = f'gridmet_{var}_1979-2020'    # .nc is added later

# %%
# Open a filesystem
fs = fsspec.filesystem('file')

# Create NCO object
nco = Nco(debug=False)

# %%
flist = sorted(fs.glob(f'{src_dir}/*_gm_{var}_*.nc'))
print(flist[0])
print(flist[-1])
len(flist)

# %%
# ncatted -a history,global,a,c,’Data version 2.0\n’ in.nc

# %% [markdown]
# The original retrieved Gridmet data (on denali) was saved in netCDF Classic format but still included the _ChunkSizes attribute which confuses some tools. Removing the attribute can be done in-place if you have read-write access to the files; otherwise you have to make a copy of the original files without the _ChunkSizes attribute.

# %%
# Original NCO command
# ncatted -a _ChunkSizes,,d,, ${ff}

# opts = ['-a _ChunkSizes,,d,,']
# nco.ncatted(input=somefile.nc)

opts = ['-h', Atted(mode='delete', att_name='_ChunkSizes')]

for ff in flist[0:2]:
    nco.ncatted(input=ff, options=opts)

# %% [markdown]
# The retrieved Gridmet data has a fixed time dimension named day. In order to concatenate individual files (1 year, 1 variable per file) into a single file per variable for the period of record we need to make the time dimension into a record dimension (e.g. unlimited).

# %%
# %%time
# ncks -O --mk_rec_dmn day ${ff} -o merged/${ff}

# Adding -h prevents adding entries into the global history


record_dim = 'day'
opts = [f'-O -h --mk_rec_dmn {record_dim}']

for ff in flist[0:2]:
    tmp_file = f'tmp_{os.path.basename(ff)}'
    print(f'Processing {tmp_file}')
    
    nco.ncks(input=ff, output=f'{output_dir}/{tmp_file}', options=opts)

# %% [markdown]
# The individual files are then concatenated together.

# %%
# %%time
# ncrcat ${var}.nc -o gridmet_${var}_1979-2020.nc

o_flist = sorted(fs.glob(f'{output_dir}/*.nc'))
opts = ['-h']

nco.ncrcat(input=o_flist, output=f'{output_dir}/{var}.nc', options=opts)



# %%
# Construct a simplified/sanitized history entry
input_files = ' '.join([os.path.basename(ff) for ff in o_flist])
history_text = f'ncrcat {input_files} -o {var}.nc'

# ncatted -a history,global,a,c,’Data version 2.0\n’ in.nc
opts = ['-h', Atted(mode='append', att_name='history', var_name='global', value=history_text, stype='c')]
# opts = [f'-h -a history,global,a,c,"{history_text}"']
# xx = nco.ncatted(input=f'{output_dir}/{var}.nc', output=f'{output_dir}/{var}_crap.nc', options=opts)
xx = nco.ncatted(input=f'{output_dir}/{var}.nc', options=opts)

# %% [markdown]
# ### Now rechunk the concatenated file.

# %%
# %%time
# ncks -O -4 -L 2 --cnk_map=dmn --cnk_dmn day,122 --cnk_dmn lat,98 --cnk_dmn lon,231 ${ff} -o ../${ff}
time_cnk = 122
lat_cnk = 98
lon_cnk = 231

opts = ['-O', '-4', '-L 2', 
        '--cnk_map=dmn', 
        f'--cnk_dmn {record_dim},{time_cnk}',
        f'--cnk_dmn lat,{lat_cnk}', 
        f'--cnk_dmn lon,{lon_cnk}']

nco.ncks(input=f'{output_dir}/{var}.nc', output=f'{output_dir}/{var}_nc4.nc', options=opts)


# %%

# %%
ds = xr.open_dataset(f'{output_dir}/{output_filename}.nc', chunks='auto')
ds

# %%
ds.daily_minimum_temperature

# %%
ds.attrs.keys()

# %%

# %%

# %%

# %% [markdown]
# ## Another way of doing this

# %% [markdown]
# Here we create a temporary directory to contain the intermediate files. The directory and its contents are automatically removed once the code block has completed. The final output file is written to the output directory.

# %%
# %%time

clean_history = True

with tempfile.TemporaryDirectory() as tmp_dir:
    print(f'Working in {tmp_dir}')
    
    # Make the time dimension a record dimension
    print('    Make record dimension')
    record_dim = 'day'
    opts = [f'-O --mk_rec_dmn {record_dim}']

    for ff in flist[0:2]:
        tmp_filename = f'tmp_{os.path.basename(ff)}'
        print(f'Processing {tmp_filename}')

        nco.ncks(input=ff, output=f'{tmp_dir}/{tmp_filename}', options=opts)
        
    # Concatenate the individual files
    print('    Concatenate files')
    o_flist = sorted(fs.glob(f'{tmp_dir}/*.nc'))

    nco.ncrcat(input=o_flist, output=f'{tmp_dir}/{var}_concat.nc')   # , options=opts)        
        
    # Rechunk the data
    print('    Rechunk data')
    time_cnk = 122
    lat_cnk = 98
    lon_cnk = 231

    opts = ['-O', '-4', '-L 2', 
            '--cnk_map=dmn', 
            f'--cnk_dmn {record_dim},{time_cnk}',
            f'--cnk_dmn lat,{lat_cnk}', 
            f'--cnk_dmn lon,{lon_cnk}']

    nco.ncks(input=f'{tmp_dir}/{var}_concat.nc', output=f'{output_dir}/{output_filename}.nc', options=opts)
    

    if clean_history:
        print('    Clean history')
        # Read the final file and create a modified history
        ds = xr.open_dataset(f'{output_dir}/{output_filename}.nc', chunks='auto')

        history = ds.attrs['History']
        aa = history.split('\n')

        new_hist = []

        for io, dd in enumerate(aa):
            bb = dd.split()
            keep = True

            for ii, cc in enumerate(bb):
                if cc in ['ncatted']:
                    keep = False
                    continue
                elif '--output' in cc:
                    bb[ii] = '--output=' + os.path.basename(cc.split('=')[1])
                else:
                    if os.path.isfile(cc):
                        bb[ii] = os.path.basename(cc)
            if keep:
                new_hist.append(' '.join(bb))

        ds.close()
        
        # Remove History global attribute
        opts = ['-h', Atted(mode='delete', att_name='History', var_name='global')]
        nco.ncatted(input=f'{output_dir}/{output_filename}.nc', options=opts)        
        
        # Add new history global attribute (all lowercase)
        # NOTE: the double backslash for the value argument is needed so that \n is passed to ncatted correctly
        opts = ['-h', Atted(mode='create', att_name='history', var_name='global', value='\\n'.join(new_hist), stype='c')]
        nco.ncatted(input=f'{output_dir}/{output_filename}.nc', options=opts)        


# %%
