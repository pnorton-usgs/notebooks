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
#     display_name: Python [conda env:bandit_38]
#     language: python
#     name: conda-env-bandit_38-py
# ---

# %%
import math

# %%
# chunk_shape_3D
#
# Tries to work out best chunking shape
#
# Ref: http://www.unidata.ucar.edu/blogs/developer/en/entry/chunking_data_choosing_shapes
#
# Source: http://www.unidata.ucar.edu/staff/russ/public/chunk_shape_3D.py
#
# Author: Russ Rew, Unidata

import math
import operator

def binlist(n, width=0):
    """Return list of bits that represent a non-negative integer.

    n      -- non-negative integer
    width  -- number of bits in returned zero-filled list (default 0)
    """
    return map(int, list(bin(n)[2:].zfill(width)))

def numVals(shape):
    """Return number of values in chunk of specified shape, given by a list of dimension lengths.

    shape -- list of variable dimension sizes"""
    if(len(shape) == 0):
        return 1
    return reduce(operator.mul, shape)

def perturbShape(shape, onbits):
    """Return shape perturbed by adding 1 to elements corresponding to 1 bits in onbits

    shape  -- list of variable dimension sizes
    onbits -- non-negative integer less than 2**len(shape)
    """
    return map(sum, zip(shape, binlist(onbits, len(shape))))

def chunk_shape_3D(varShape, valSize=4, chunkSize=4096):
    """
    Return a 'good shape' for a 3D variable, assuming balanced 1D, 2D access

    varShape  -- length 3 list of variable dimension sizes
    chunkSize -- maximum chunksize desired, in bytes (default 4096)
    valSize   -- size of each data value, in bytes (default 4)

    Returns integer chunk lengths of a chunk shape that provides
    balanced access of 1D subsets and 2D subsets of a netCDF or HDF5
    variable var with shape (T, X, Y), where the 1D subsets are of the
    form var[:,x,y] and the 2D slices are of the form var[t,:,:],
    typically 1D time series and 2D spatial slices.  'Good shape' for
    chunks means that the number of chunks accessed to read either
    kind of 1D or 2D subset is approximately equal, and the size of
    each chunk (uncompressed) is no more than chunkSize, which is
    often a disk block size.
    """

    rank = 2  # this is a special case of n-dimensional function chunk_shape
    chunkVals = chunkSize / float(valSize) # ideal number of values in a chunk
    numChunks  = varShape[0] * varShape[1] / chunkVals # ideal number of chunks

    axisChunks = numChunks ** 0.5 # ideal number of chunks along each 2D axis
    cFloor = [] # will be first estimate of good chunk shape
    # cFloor  = [varShape[0] // axisChunks**2, varShape[1] // axisChunks, varShape[2] // axisChunks]
    
    # Except that each chunk shape dimension must be at least 1
    # chunkDim = max(1.0, varShape[0] // axisChunks**2)
    print('chunkVals = {}'.format(chunkVals))
    print('numChunks = {}'.format(numChunks))
#     print('axisChunks = {}'.format(axisChunks))
    
    if varShape[0] / axisChunks**2 < 1.0:
        chunkDim = 1.0
        axisChunks = axisChunks / math.sqrt(varShape[0] / axisChunks**2)
    else:
        chunkDim = varShape[0] // axisChunks**2
    
    print('axisChunks**2 = {}'.format(axisChunks**2))
    print('axisChunks = {}'.format(axisChunks))
    print('chunkDim = {}'.format(chunkDim))
    
    cFloor.append(chunkDim)
    prod = 2.0  # factor to increase other dims if some must be increased to 1.0
    
    for i in range(1, rank):
        if varShape[i] / axisChunks < 1.0:
            prod *= axisChunks / varShape[i]
            print('Changed prod: {}'.format(prod))
    for i in range(1, rank):
        if varShape[i] / axisChunks < 1.0:
            chunkDim = 1.0
        else:
            chunkDim = (prod * varShape[i]) // axisChunks
        cFloor.append(chunkDim)

    # cFloor is typically too small, (numVals(cFloor) < chunkSize)
    # Adding 1 to each shape dim results in chunks that are too large,
    # (numVals(cCeil) > chunkSize).  Want to just add 1 to some of the
    # axes to get as close as possible to chunkSize without exceeding
    # it.  Here we use brute force, compute numVals(cCand) for all
    # 2**rank candidates and return the one closest to chunkSize
    # without exceeding it.
    bestChunkSize = 0
    cBest = cFloor
    print('chunkSize = {}'.format(chunkSize))
    print('cBest = {}'.format(cBest))
    
    for i in range(8):
        # cCand = map(sum,zip(cFloor, binlist(i, rank)))
        cCand = perturbShape(cFloor, i)
        thisChunkSize = valSize * numVals(cCand)
        print('thisChunkSize = {}'.format(thisChunkSize))
        print(cCand)
        
        if bestChunkSize < thisChunkSize <= chunkSize:
            bestChunkSize = thisChunkSize
            cBest = list(cCand) # make a copy of best candidate so far
    return map(int, cBest)




# %%
def chunk_shape_2D(varShape, valSize=4, chunkSize=4096):
    num_vals_in_chunk = float(chunkSize / valSize)
    
    num_vals_total = varShape[0] * varShape[1]
    
    chunk_percent = (num_vals_in_chunk / num_vals_total)**.5
    
    if chunk_percent > 1.0:
        print('ERROR: Total size of array is too small for reliable chunks at this chunkSize')
    # print('num_vals_in_chunk = {}'.format(num_vals_in_chunk))
    # print('num_vals_total = {}'.format(num_vals_total))
    # print('chunk_percent = {}'.format(chunk_percent))
    
    starting_chunk = [int(varShape[0] * chunk_percent), int(varShape[1] * chunk_percent)]
    # starting_chunk = [int(math.ceil(varShape[0] * chunk_percent)), int(math.ceil(varShape[1] * chunk_percent))]
    starting_size = starting_chunk[0] * starting_chunk[1]
    # print('starting_chunk = {}'.format(starting_chunk))
    # print('starting_size = {}'.format(starting_size))
    
    curr_best = starting_chunk
    for x1 in range(0, 2):
        for x2 in range(0, 2):
            # print(starting_chunk[0] + x1, starting_chunk[1] + x2)
            cnk_size = (starting_chunk[0] + x1) * (starting_chunk[1] + x2)
            # print('  size = {}'.format(cnk_size))
            
            if starting_size < cnk_size <= num_vals_in_chunk:
                curr_best = [starting_chunk[0] + x1, starting_chunk[1] + x2]
                
    # print('curr_best = {}'.format(curr_best))
    return curr_best

def compute_chunk_cache(chunk_sizes, dim_size, var_type):
    chunks_per_domain = math.ceil(float(dim_size) / float(chunk_sizes[0]))

    # Compute initial cache size
    # var_cache_size = 2**exponent(real(cnk_per_domain * product(cnk_sizes) * 4))
    var_cache_size = None
    if var_type == 'f4':
        var_cache_size = 2**math.frexp(float(chunks_per_domain * chunk_sizes[0] * chunk_sizes[1] * 4))[1]
    elif var_type == 'f8':
        var_cache_size = 2**math.frexp(float(chunks_per_domain * chunk_sizes[0] * chunk_sizes[1] * 8))[1]
    return var_cache_size


# %%
# NCO chunking discussion: https://github.com/nco/nco/issues/70

# chunk_shape_3D(varShape, valSize=4, chunkSize=4096)
# chunk_shape_2D([14,3], chunkSize=192)
# chunk_shape_2D([109951, 13241], chunkSize=32768)
# chunk_shape_2D([3065, 7305], chunkSize=32768)
# chunk_shape_2D([13515, 109951], chunkSize=16384)

# dimsizes [13241, 109951]
# 32768 [31, 262]
# 65536 [44, 369]
# 262144 [88, 738]
# 524288 [125, 1044]
# 1048576 [177, ]
# 2097152 [251, 2087]

chunk_shape_2D([7814, 8075], valSize=2, chunkSize=32768)

# %%

# %%
compute_chunk_cache([31, 262], 13241, 'f4')

# %%
nhru = 366
ntime = 114958
dst_cnk_size = 40 * 1024 * 1024   # in bytes

cnk_sizes = chunk_shape_2D([nhru, ntime], valSize=4, chunkSize=dst_cnk_size)
cnk_size_mb = (cnk_sizes[0] * cnk_sizes[1] * 4) / 1024 / 1024
cnk_per_domain = math.ceil(float(nhru) / float(cnk_sizes[1]))
var_cache_size = compute_chunk_cache(cnk_sizes, nhru, 'f4')
var_cache_size_mb = var_cache_size / 1024 / 1024
var_cache_slots = cnk_per_domain * 10

print(f'{cnk_sizes=}')
print(f'{cnk_size_mb=} MB')
print(f'{cnk_per_domain=}')
print(f'{var_cache_size=} bytes')
print(f'{var_cache_size_mb=} MB')
print(f'{var_cache_slots=}')

# %%
cnk_sizes

# %%
9%4

# %%
nhru = 114958
for xx in range(1, nhru):
    if nhru % xx == 0:
        print(f'{xx}')

# %%
math.sqrt(nhru)


# %%
def prime3(a):
    if a < 2: return False
    if a == 2 or a == 3: return True # manually test 2 and 3   
    if a % 2 == 0 or a % 3 == 0: return False # exclude multiples of 2 and 3
 
    maxDivisor = a**0.5
    d, i = 5, 2
    while d <= maxDivisor:
        if a % d == 0: return False
        d += i 
        i = 6 - i # this modifies 2 into 4 and viceversa
 
    return True


# %%
prime3(2131)

# %%
import datetime

st_date = datetime.datetime(1980,10,1)
en_date = datetime.datetime(1989,9,30)

# %%
(en_date - st_date).days + 1

# %%
import math

ex = math.frexp(78375600.)

2**ex[1]

# %%
