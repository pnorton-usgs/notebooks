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
# original script at: https://bitbucket.csiro.au/projects/CMAR_RS/repos/netcdf-tools/browse/chunking/chunk_shape_3D.py?at=97c17c836371b07d9d08c9ffa05b3a8db623e0f1

# chunk_shape_3D: Tries to work out best chunking shape
# Author: Russ Rew, Unidata

#-this.is.my.updated.version.(August.06.2018).for.the.script.to.run.in.Python3x


import math
import operator
from functools import reduce


def binlist(n, width=0):
    """
    Return list of bits that represent a non-negative integer.
    n      -- non-negative integer
    width  -- number of bits in returned zero-filled list (default 0)
    """
    return list(map(int, list(bin(n)[2:].zfill(width))))

def numVals(shape):
    """
    Return number of values in chunk of specified shape, given by a list of dimension lengths.
    shape -- list of variable dimension sizes"""
    if(len(shape) == 0):
        return 1
    return reduce(operator.mul, shape)

def perturbShape(shape, onbits):
    """
    Return shape perturbed by adding 1 to elements corresponding to 1 bits in onbits
    shape  -- list of variable dimension sizes
    onbits -- non-negative integer less than 2**len(shape)
    """
    return list(map(sum, zip(shape, binlist(onbits, len(shape)))))

def chunk_shape_3D(varShape, valSize=4, chunkSize=4096):
    """
    Return a 'good shape' for a 3D variable, assuming balanced 1D, 2D access
    varShape  -- length 3 list of variable dimension sizes
    chunkSize -- maximum chunksize desired, in bytes (default 4096)
    valSize   -- size of each data value, in bytes (default 4)
    Returns integer chunk lengths of a chunk shape that provides balanced access of 1D subsets
    and 2D subsets of a netCDF or HDF5 variable var with shape (T, X, Y), where the 1D subsets
    are of the form var[:,x,y] and the 2D slices are of the form var[t,:,:]. typically 1D time
    series and 2D spatial slices.
    'Good shape' for chunks means that the number of chunks accessed to read either kind of 1D
    or 2D subset is approximately equal, and the size of each chunk (uncompressed) is no more
    than chunkSize, which isoften a disk block size.
    """
    rank = 3  # this is a special case of n-dimensional function chunk_shape
    chunkVals = chunkSize / float(valSize) # ideal number of values in a chunk
    numChunks = varShape[0] * varShape[1] * varShape[2] / chunkVals # ideal number of chunks
    axisChunks = numChunks**0.25 # ideal number of chunks along each 2D axis
    cFloor = [] # will be first estimate of good chunk shape
    
    
    print(f'chunkVals = {chunkVals:0.1f}')
    print(f'numChunks = {numChunks:0.1f}')
    print(f'axisChunks = {axisChunks:0.1f}')
    
    #cFloor = [varShape[0]//axisChunks**2,varShape[1]//axisChunks,varShape[2]//axisChunks]
    #except that each chunk shape dimension must be at least 1
    #chunkDim = max(1.0,varShape[0]//axisChunks**2)
    if varShape[0] / axisChunks**2 < 1.0:
        chunkDim = 1.0
        axisChunks = axisChunks / math.sqrt(varShape[0] / axisChunks**2)
    else:
        chunkDim = varShape[0] // axisChunks**2
        
    print(f'axisChunks**2 = {axisChunks**2:0.3f}')
        
    print(f'axisChunks = {axisChunks:0.3f}')
    print(f'chunkDim = {chunkDim:0.3f}')
    
    cFloor.append(chunkDim)
    prod = 1.0  # factor to increase other dims if some must be increased to 1.0
    for i in range(1, rank):
        if varShape[i] / axisChunks < 1.0:
            prod *= axisChunks / varShape[i]
            
    for i in range(1, rank):
        if varShape[i] / axisChunks < 1.0:
            chunkDim = 1.0
        else:
            chunkDim = (prod * varShape[i]) // axisChunks
        cFloor.append(chunkDim)
        
    """
    cFloor is typically too small, (numVals(cFloor) < chunkSize) Adding 1 to each shape dim
    results in chunks that are too large, (numVals(cCeil) > chunkSize).  Want to just add 1
    to some of the axes to get as close as possible to chunkSize without exceeding it. Here 
    we use brute force, compute numVals(cCand) for all 2**rank candidates and return the one
    closest to chunkSize without exceeding it.
    """
    bestChunkSize = 0
    cBest = cFloor
    print('--- Candidates ---')
    for i in range(8):
        # cCand = map(sum,zip(cFloor,binlist(i,rank)))
        cCand = perturbShape(cFloor, i)
#         print(cCand)
        
        thisChunkSize = int(valSize * numVals(cCand))
        print(f'{list(map(int, cCand))}: Total size per chunk: {thisChunkSize} (ratio: {numVals(cCand) / chunkVals:0.3f})')
        
        if bestChunkSize < thisChunkSize <= chunkSize:
            bestChunkSize = thisChunkSize
            cBest = list(cCand) # make a copy of best candidate so far
    return list(map(int, cBest))

# %%

# sz_dims = [8759, 224, 464]
# sz_dims = [672, 1015, 1367]    # Total size of each dimension
# sz_dims = [1, 1015, 1367]
# sz_dims = [1512, 596, 1385]    # USGS Monthly Water Balance Model
sz_dims = [12, 8075, 7814]
# 1048576 = 1M
# 4194304 = 4M
# 8388608 = 8M
# 10485760 = 10M
# 104857600 = 100M

block_size = 1048576 # Bytes 
data_size = 4    # Bytes  (float32 = 4, float_64 = 8)
D = block_size / data_size    # max number of values in a chunk
print(f'Number of values per block: {D}')

print('='*40)
sel_chunks = chunk_shape_3D(sz_dims, valSize=data_size, chunkSize=block_size)
print('='*40)

sel_cnk_num_vals = sel_chunks[0] * sel_chunks[1] * sel_chunks[2]
sel_cnk_sz = sel_cnk_num_vals * data_size
ratio = sel_cnk_num_vals / D

print(f'Selected chunk sizes: {sel_chunks}')
print(f'Number of values in chunk: {sel_cnk_num_vals}, size: {sel_cnk_sz}')
print(f'Efficiency ratio: {ratio:0.3f}')

# Size of each chunk
cnk_sz_MB = sel_cnk_sz / 1024 / 1024

# Time-series
num_cnk_per_access_ts = math.ceil(sz_dims[0] / sel_chunks[0])

# Spatial slice
num_cnk_per_access_spatial = math.ceil(sz_dims[1] / sel_chunks[1]) * math.ceil(sz_dims[2] / sel_chunks[2])

print(f'Time-series: {num_cnk_per_access_ts}; read: {num_cnk_per_access_ts * cnk_sz_MB:0.1f} MB')
print(f'Spatial slice: {num_cnk_per_access_spatial}; read: {num_cnk_per_access_spatial * cnk_sz_MB:0.1f} MB')

# %%
# 4, 80, 3
# 4, 3, 80


# %% [markdown]
# Different chunk shapes are optimal for different patterns of access. For example, if you know that the data will primarily be accessed by fixed slices of the first dimension, optimal chunks would be in the shape 1 x 37 x 256 x 512 or a similar shape that fits within one physical disk block, such as 1 x 5 x 32 x 64, which for 4-byte values fits within a 64K disk block.
#
# If there is no most common access pattern, let D be the number of values you want in a chunk (preferably <= what one disk access can read), then let c = (D/(25256 * 37 * 256 * 512)) ^ (1/4), and use a chunk shape resulting from multiplying each dimension size by c and truncating to an integer.

# %%
# sz_dim1 = 168*40*52
# sz_dim2 = 1015
# sz_dim3 = 1367

sz_dim1 = 672
sz_dim2 = 1015
sz_dim3 = 1367

# sz_dim1 = 1512
# sz_dim2 = 596
# sz_dim3 = 1385


# 1048576 = 1M
# 4194304 = 4M
# 104857600 = 100M

block_size = 10485760
data_size = 4    # bytes
D = block_size / data_size    # number of values in a block
# c = D / (sz_dim1 * sz_dim2 * sz_dim3)**0.25

# %%
print(f'Number of values in a block: {D}')

# %%
# cnk_dim1 = 168
# cnk_dim2 = 350
# cnk_dim3 = 350

cnk_dim1 = 72
cnk_dim2 = 175
cnk_dim3 = 175

# cnk_dim1 = 120
# cnk_dim2 = 300
# cnk_dim3 = 700

dim1_cnk_sz = cnk_dim1 * data_size
dim2_cnk_sz = cnk_dim2 * data_size
dim3_cnk_sz = cnk_dim3 * data_size

this_chunk_num_vals = cnk_dim1 * cnk_dim2 * cnk_dim3
this_chunk_sz = this_chunk_num_vals * data_size
# ratio = this_chunk_sz / block_size
ratio = this_chunk_num_vals / D

print(f'Number of values in chunk: {this_chunk_num_vals}, size: {this_chunk_sz}')
print(f'Efficiency ratio: {ratio}')


print(f'{dim1_cnk_sz=}')
print(f'{dim2_cnk_sz=}')
print(f'{dim3_cnk_sz=}')

# %%
print(sz_dim2 / 3)
print(sz_dim3 / 3)

print(cnk_dim2 / sz_dim2)
print(cnk_dim3 / sz_dim3)

# %%
sz_dim2 % cnk_dim2

# %%
sz_dim2 / cnk_dim2

# %%

# %%
# Number of chunks per access
cnk_sz_MB = this_chunk_sz / 1024 / 1024

# Time-series
num_cnk_per_access_ts = math.ceil(sz_dim1 / cnk_dim1)

# Spatial slice
num_cnk_per_access_spatial = math.ceil(sz_dim2 / cnk_dim2) * math.ceil(sz_dim3 / cnk_dim3)

print(f'Time-series: {num_cnk_per_access_ts}; read: {num_cnk_per_access_ts * cnk_sz_MB} MB')
print(f'Spatial slice: {num_cnk_per_access_spatial}; read: {num_cnk_per_access_spatial * cnk_sz_MB} MB')

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
