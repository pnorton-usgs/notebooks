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

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
from collections import OrderedDict
from pyPRMS.ParamDb import ParamDb
from pyPRMS.ParameterFile import ParameterFile
from Bandit.bandit_helpers import parse_gages, set_date, subset_stream_network
import networkx as nx

# %%
basedir = '/Users/pnorton/Projects/National_Hydrology_Model'
paramdb = f'{basedir}/datasets/paramdb_v11/paramdb_v11_gridmet_CONUS'

workdir = f'{basedir}/FY2021_projects/water_use/models/20201113_gm_09380000'
filename = f'{workdir}/myparam.param'

# %%
# Load the parameter database
pdb = ParamDb(paramdb_dir=paramdb, verbose=True, verify=True)
nhm_params = pdb.parameters

# Load parameter file
pfile = ParameterFile(filename, verbose=True, verify=True)
params = pfile.parameters

# %%
poi_ids = params['poi_gage_id'].tolist()
poi_segments = params['poi_gage_segment'].tolist()


hru_segment = nhm_params.get('hru_segment_nhm').tolist()
nhm_id = nhm_params['nhm_id'].tolist()
nhm_id_to_idx = nhm_params.get('nhm_id').index_map
nhm_seg = nhm_params['nhm_seg'].tolist()
tosegment_nhm = nhm_params['tosegment_nhm'].tolist()
poi_id_to_seg = nhm_params.poi_to_seg

print(f'Number of NHM hru_segment entries: {len(hru_segment)}')

# %% [markdown]
# ## Build the stream network

# %%
# Build the stream network
dag_ds = pdb.parameters.stream_network(tosegment='tosegment_nhm', seg_id='nhm_seg')

# %%
# Read the POIs from a given model and produce a CSV file containing each POIs
# and the associated number of segments and HRUs. 
outfile = open(f'{workdir}/poi_hrus.csv', 'w')
outfile.write('POI,segments,HRUs\n')

uscutoff_seg = list()

for pidx, poi_id in enumerate(poi_ids):
    print(pidx)
    
    dsmost_seg = [poi_id_to_seg[poi_id]]
    
    dag_ds_subset = subset_stream_network(dag_ds, uscutoff_seg, dsmost_seg)
    
    # Create list of toseg ids for the model subset
    toseg_idx = list(set(xx[0] for xx in dag_ds_subset.edges))

    # Use the mapping to create subsets of nhm_seg, tosegment_nhm, and tosegment
    # NOTE: toseg_idx and new_nhm_seg are the same thing
    new_nhm_seg = [ee[0] for ee in dag_ds_subset.edges]

    # Using a dictionary mapping nhm_seg to 1-based index for speed
    new_nhm_seg_to_idx1 = OrderedDict((ss, ii+1) for ii, ss in enumerate(new_nhm_seg))

    # Generate the renumbered local tosegments (1-based with zero being an outlet)
    new_tosegment = [new_nhm_seg_to_idx1[ee[1]] if ee[1] in new_nhm_seg_to_idx1
                     else 0 for ee in dag_ds_subset.edges]

    # Create a dictionary mapping hru_segment segments to hru_segment 1-based indices filtered
    # by new_nhm_seg and hru_noroute.
    seg_to_hru = OrderedDict()
    hru_to_seg = OrderedDict()

    for ii, vv in enumerate(hru_segment):
        # Contains both new_nhm_seg values and non-routed HRU values
        # keys are 1-based, values in arrays are 1-based
        if vv in new_nhm_seg:
            hid = nhm_id[ii]
            seg_to_hru.setdefault(vv, []).append(hid)
            hru_to_seg[hid] = vv

    # HRU-related parameters can either be output with the legacy, segment-oriented order
    # or can be output maintaining their original HRU-relative order from the parameter database.
    hru_order_subset = [kk for kk in hru_to_seg.keys()]

    new_hru_segment = [new_nhm_seg_to_idx1[kk] if kk in new_nhm_seg else 0 if kk == 0 else -1 for kk in
                       hru_to_seg.values()]

    hru_order_subset0 = [nhm_id_to_idx[xx] for xx in hru_order_subset]
    outfile.write(f'{poi_id},{len(toseg_idx)},{len(hru_order_subset)}\n')

outfile.close()

# %%

# %%
