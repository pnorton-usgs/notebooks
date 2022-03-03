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

# %%

# %%
def rawcount(filename):
    f = open(filename, 'rb', buffering=0)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read
    # read_f = f.raw.read

    buf = read_f(buf_size)
    while buf:
        lines += buf.count(b'\n')
        buf = read_f(buf_size)

    return lines


# %%
work_dir = '/Volumes/parker_rocks/NHM_output'
filename = f'{work_dir}/nhru_transp_on.csv'

# %%
# 2.92 s ± 103 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
# # %timeit rawcount(filename)

# wc -l
# 0m0.581s

rawcount(filename)

# %%
filename = f'{work_dir}/nhru_swrad.csv'

# %%
# (after caching by calling wc -l)
# 13.3 s ± 206 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
# %timeit rawcount(filename)

# wc -l 
# 2m49.721s
# after caching
# 0m2.495s

# %%
rawcount(filename)

# %%
