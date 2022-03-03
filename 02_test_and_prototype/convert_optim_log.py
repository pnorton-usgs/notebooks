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
#     display_name: Python 2
#     language: python
#     name: python2
# ---

# %%
# %matplotlib inline
import matplotlib.pyplot as plt

import prms_lib as prms
import pandas as pd
import numpy as np
import math as mth
import datetime
import re

#print len('Determining starting parameters...')
print len('Current generation for generation ')
print len('Results for multi-objective global optimization:')

# %% [markdown]
# ### Read optim_log.log and convert to a CSV format

# %%
basedir = '/Users/pnorton/Projects/National_Hydrology_Model/PRMS_testing/mocom_t4'
basinid = '4.39_03410500'
runid = '2015-04-29_1424'

workdir = '%s/%s/runs/%s' % (basedir, basinid, runid)
opt_file = '%s/optim_log.log' % workdir
opt_out = '%s/optim_fixed.log' % workdir

infile = open(opt_file, 'r')
rawdata = infile.read().splitlines()
infile.close()
it = iter(rawdata)

outfile = open(opt_out, 'w')

bad_chars = '():='
rgx = re.compile('[%s]' % bad_chars)

for line in it:
    if line[0:34] == 'Determining starting parameters...':
        # This is the group of random starting sets
        next(it)
        gennum = 0
        
        tmp_hdr = next(it).split()
        tmp_hdr.insert(0,'setnum')
        hdr_flag = True
        
        print "header length: %d" % len(tmp_hdr)
        #tmp = 'setnum ' + next(it) + ' test0 test1 rank soln_num gennum'
        #outfile.write('%s\n' % ','.join(tmp.split()))
        #print tmp.split()
        
        while True:
            x = next(it)
            if x[0:1] == '':
                break

            # Strip out the junk characters ():=
            x = rgx.sub('', x) + ' ' + str(gennum)
            x = x.split()
            
            if hdr_flag:
                # Header info from starting population is incomplete, fill it it out
                # with information inferred from the first line of data
                print "First date line length: %d" % len(x)
                for pp in range(0,(len(x) - len(tmp_hdr) - 3)):
                    tmp_hdr.append('test%d' % pp)
                tmp_hdr.append('rank')
                tmp_hdr.append('soln_num')
                tmp_hdr.append('gennum')
                
                # Write the header out to the file
                outfile.write('%s\n' % ','.join(tmp_hdr))
                #print tmp_hdr

                hdr_flag = False
            
            outfile.write('%s\n' % ','.join(x))
            #print x.split()
            
    if line[0:34] == 'Current generation for generation ':
        gennum = int(line.split(' ')[-1].rstrip(':'))+1
        #print 'gennum:',gennum
        next(it)    # skip one line
        next(it)    # skip one line
        next(it)    # skip one line
        
        while True:
            x = next(it)
            if x[0:1] == '':
                break
            
            # Strip out the junk characters ():=
            x = rgx.sub('', x) + ' ' + str(gennum)
            
            outfile.write('%s\n' % ','.join(x.split()))
            #print x.split()
    elif line[0:48] == 'Results for multi-objective global optimization:':
        gennum = int(next(it).split()[1])+1
        print 'gennum:',gennum
        next(it)    # skip one line
        next(it)    # skip one line
        next(it)    # skip one line
        
        while True:
            x = next(it)
            if x[0:1] == '':
                break
            
            # Strip out the junk characters ():=
            x = rgx.sub('', x) + ' ' + str(gennum)
            
            outfile.write('%s\n' % ','.join(x.split()))
outfile.close()
print 'Total generations:', gennum - 1

# %%

# %%
