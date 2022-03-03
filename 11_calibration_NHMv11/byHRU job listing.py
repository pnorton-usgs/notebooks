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

# %%
import datetime

# %%

# %%
workdir = '/Users/pnorton/Projects/National_Hydrology_Model/calibrations/NHMv11'

notedir = f'{workdir}/notes'
slurmdir = f'{workdir}/slurms'

# %%
# Read file containing directory listing of slurm job outputs
#   created with: find . -mindepth 1 -maxdepth 1 -name "*-byHRUall.out" -printf '%f %Tx\n' > ~/slurms_dir_list.txt
fhdl = open(f'{notedir}/slurms_dir_list.txt')

jobdates = {}

for row in fhdl:
    flds = row.strip().split(' ')
    
    jobid = int(flds[0].split('-')[0])
    jdate = datetime.datetime.strptime(flds[1], '%m/%d/%Y')
    
    jobdates[jobid] = jdate

# %%

# %%

# %%
fhdl = open(f'{notedir}/byHRU_jobs.txt', 'r')

job_to_hrus = {}
byHRU_job = {}
byHRU_dir = {}

for row in fhdl:
    # Sample row: /scratch/5288060/DIR24_4790/HRU114938
    flds = row.strip().split('/')
    
    jobid = int(flds[2])
    jobdir = flds[3]
    hru = int(flds[4][3:])
    
    if hru not in byHRU_job:
        byHRU_job[hru] = set()
    byHRU_job[hru].add(jobid)
    # byHRU_job[hru].append(jobid)
    
    if hru not in byHRU_dir:
        byHRU_dir[hru] = set()
    byHRU_dir[hru].add(jobdir)
    # byHRU_dir[hru].append(jobdir)
    
    if jobid not in job_to_hrus:
        job_to_hrus[jobid] = []
    job_to_hrus[jobid].append(hru)

# %%
job_to_hrus

# %%
# Open file of missing HRU ids
fhdl = open(f'{workdir}/v11_byHRU_missing.csv', 'r')

# skip first line
fhdl.readline()

job_to_hrus = {}
miss_hru_jobs = {}

for row in fhdl:
    flds = row.strip().split(',')
    
    last_job = None
    last_date = datetime.datetime(2000,10,1)

    mhru = int(flds[0])
    
    try:
        for jj in byHRU_job[mhru]:
            if jobdates[jj] > last_date:
                last_job = jj
                last_date = jobdates[jj]

        miss_hru_jobs[mhru] = last_job
        
        if last_job not in job_to_hrus:
            job_to_hrus[last_job] = []
        job_to_hrus[last_job].append(mhru)
    except KeyError:
        print(f'Missing HRU, {mhru}, not found in the job list')

fhdl.close()

# %%
hru = 112644

last_job = None
last_date = datetime.datetime(2000,10,1)

for jj in byHRU_job[hru]:
    if jobdates[jj] > last_date:
        last_job = jj
        last_date = jobdates[jj]
        
print(f'Calibration jobs: {byHRU_job[hru]}')
print(f' Most-recent job: {last_job}')
print(f'   Job directory: {byHRU_dir[hru]}')


# %%
# Check the most recent job for any errors
job_to_hrus

# %%
miss_hru_jobs

# %%
ohdl = open(f'{workdir}/job_check_list.txt', 'w')
ohdl.write('nhm_id,jobid\n')

for kk, vv in miss_hru_jobs.items():
    ohdl.write(f'{kk},{vv}\n')
ohdl.close()

# %%
ohdl = open(f'{workdir}/job_to_hrus.txt', 'w')
ohdl.write('jobid,nhm_ids\n')

for kk, vv in job_to_hrus.items():
    ohdl.write(f'{kk}\t{vv}\n')
ohdl.close()

# %%
