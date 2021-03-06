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
#     display_name: Python [conda env:py27]
#     language: python
#     name: conda-env-py27-py
# ---

# %%
from future.utils import iteritems

import os
import subprocess32
import sys
from collections import OrderedDict

# from Queue import *
# from threading import Thread

import bandit_cfg as bc


# %%
import os, time
import threading, Queue
import subprocess

class WorkerThread(threading.Thread):
    """ A worker thread that takes directory names from a queue, finds all
        files in them recursively and reports the result.

        Input is done by placing directory names (as strings) into the
        Queue passed in dir_q.

        Output is done by placing tuples into the Queue passed in result_q.
        Each tuple is (thread name, dirname, [list of files]).

        Ask the thread to stop by calling its join() method.
    """
    def __init__(self, input_q, result_q):
        super(WorkerThread, self).__init__()
        self.input_q = input_q
        self.result_q = result_q
        self.stoprequest = threading.Event()

    def run(self):
        # As long as we weren't asked to stop, try to take new tasks from the
        # queue. The tasks are taken with a blocking 'get', so no CPU
        # cycles are wasted while waiting.
        # Also, 'get' is given a timeout, so stoprequest is always checked,
        # even if there's nothing in the queue.
        while not self.stoprequest.isSet():
            try:
                thecmd = self.input_q.get(True, 0.05)
                retcode = self.run_cmd(thecmd)
                self.result_q.put((self.name, thecmd, retcode))
            except Queue.Empty:
                continue

    def join(self, timeout=None):
        self.stoprequest.set()
        super(WorkerThread, self).join(timeout)

    def run_cmd(self,thecmd):
        try:
            retcode = subprocess.call(thecmd, shell=True)
        except:
            # Error running command
            print "Error: shutting down"
            retcode = -1
            self.stoprequest.set()
        return retcode


# %%
hw_src = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/conus_headwaters/hwSegsVll.csv'

hw_file = open(hw_src, 'r')
hw_file.next()

# First column is hwAreaId
# Second and following columns are seg_id_nat
head_waters = OrderedDict()

for line in hw_file:
    cols = line.strip().replace(" ", "").split(',')
    cols = [int(xx) for xx in cols]
    head_waters[cols[0]] = cols[1:]

# %%

# %%
num_threads = 4

# ****************************************************************************
# Initialize the threads
cmd_q = Queue.Queue()
result_q = Queue.Queue()

# Create pool of threads
pool = [WorkerThread(input_q=cmd_q, result_q=result_q) for i in range(num_threads)]

# Start the threads
for thread in pool:
    try:
        thread.start()
    except (KeyboardInterrupt, SystemExit):
        # Shutdown the threads when the program is terminated
        thread.join()
        sys.exit(1)

# %% [markdown]
# ### For each head_waters
# - create directory hw_# (where # is the hwAreaId)
# - copy default bandit.cfg into directory
# - run bandit on the directory
#

# %%
default_config_file = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/conus_headwaters/bandit.cfg'

jobdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/conus_headwaters/hw_job'
hwdir = 'hw_{:04d}'

cmd_bandit = '/Users/pnorton/PycharmProjects/Bandit/bandit.py'

if not os.path.exists(jobdir):
    try:
        os.makedirs(jobdir)
    except OSError, err:
        print "\tError creating directory: %s" % err
        exit(1)
        
st_dir = os.getcwd()
os.chdir(jobdir)

# Read the default configuration file
config = bc.Cfg(default_config_file)

work_count = 0

for kk, vv in iteritems(head_waters):
    cdir = hwdir.format(kk)

    # Create the headwater directory if needed
    if not os.path.exists(cdir):
        try: 
            os.makedirs(cdir)
        except OSError, err:
            print "\tError creating directory: %s" % err
            exit(1)
            
    # Update the outlets in the basin.cfg file and write into the headwater directory
    config.update_value('outlets', vv)
    config.update_value('output_dir', '{}/{}'.format(jobdir, cdir))
    config.write('{}/bandit.cfg'.format(cdir))
    
    # Run bandit
    # Add the command to queue for processing
    work_count += 1
#     cmd_q.put(" ".join(cmd))
    
    os.chdir(cdir)
    cmd_q.put(cmd_bandit)
#     subprocess.call(cmd_bandit, shell=True)
    os.chdir(jobdir)
    
    if work_count==6:
        break
        
print("work_count = {:d}".format(work_count))

# Output results
while work_count > 0:
    result = result_q.get()

    sys.stdout.write("\r" + "work_count: {:4d}".format(work_count))
    sys.stdout.flush()
    
#     print "Thread %s return code = %d" % (result[0], result[2])
    work_count -= 1
    
    if result[2] != 0:
        # An error occurred running the command
        work_count = 0

# Ask for the threads to die and wait for them to do it
for thread in pool:
    thread.join()


# %%

# %%

# %%

# %%
