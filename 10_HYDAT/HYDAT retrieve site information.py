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
import netCDF4 as nc
import numpy as np
import pandas as pd
import sqlite3

# %%
src_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/HYDAT/20220502_update'
# src_dir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/HYDAT'
src_db = f'{src_dir}/Hydat.sqlite3'

# %%
connection = sqlite3.connect(src_db)

# %%
cursor = connection.cursor()

# %%
my_query = 'SELECT STATION_NUMBER, STATION_NAME, LATITUDE, LONGITUDE, DRAINAGE_AREA_GROSS, DRAINAGE_AREA_EFFECT, PROV_TERR_STATE_LOC FROM STATIONS'
cursor.execute(my_query)

# %%
# all_rows = cursor.fetchall()

# for row in all_rows:
#     print(row)

# %%
df = pd.read_sql_query(my_query, connection)

# %%
df.head()

# %%
df.info()

# %%
field_map = {'STATION_NUMBER': 'poi_id',
             'STATION_NAME': 'poi_name',
             'LATITUDE': 'latitude',
             'LONGITUDE': 'longitude',
             'DRAINAGE_AREA_GROSS': 'drainage_area',
             'DRAINAGE_AREA_EFFECT': 'drainage_area_contrib'}

# %%
df.rename(columns=field_map, inplace=True)
df.set_index('poi_id', verify_integrity=True, inplace=True)
df.sort_index(inplace=True)

# %%

# %%
# Add missing columns
df['poi_agency'] = 'EC'

# %%
ec_gages = ['01AG003', '01AD004', '01AD003', '01AF009', '01AF003', '01AJ004', '01AF002', '01AK001', '01AJ010', 
            '01AK004', '01AF010', '01AF007', '01AG002', '01AJ001', '01AK005', '01AD001', '01AK007', '01AJ009', 
            '01AK008', '01AJ003', '01AJ011', '01AQ008', '01AR004', '01AR005', '01AQ002', '01AQ001', '01AR011', 
            '01AR006', '02OE018', '05NB001', '05NB039', '05NF014', '05NB021', '05NB017', '05NB027', '05ND001', 
            '05NB036', '05NB035', '05NB041', '05NB031', '05NF016', '05ND007', '05NA003', '05NA005', '05NA004', 
            '05NF002', '05NF013', '05NC001', '05NF015', '05NF007', '05ND010', '05NB040', '05NB011', '05NF006', 
            '05NB018', '05NB014', '05NB034', '05ND004', '05NB009', '05NB033', '05NF001', '05ND011', '05NB030', 
            '05NF010', '05NF008', '05OB019', '05OC001', '05OD004', '05OB021', '05OD001', '05OA006', '05OA008', 
            '05OA007', '05OB010', '05OB001', '05OB016', '05OC015', '05OA010', '05OA015', '05OB006', '05OB025', 
            '05OC022', '05OB007', '05OD028', '05OB023', '05OC019', '05PD026', '05PB009', '05PC018', '05PB014', 
            '05PC022', '05PC011', '05PB021', '05PB015', '05PC019', '05PB018', '05AD007', '05AE904', '05AD010', 
            '05AE006', '05AD005', '05AE005', '05AE002', '05AD901', '05AE043', '05AD035', '05AD028', '05AD003', 
            '05AD041', '05AE041', '11AB118', '11AB103', '11AB083', '11AC023', '11AB117', '11AC001', '11AC025', 
            '11AC017', '11AB101', '11AC051', '11AB075', '11AA039', '11AA001', '11AB008', '11AB082', '11AB009', 
            '11AC062', '11AB096', '11AA029', '11AA025', '11AA005', '11AB108', '11AA026', '11AF005', '11AE010', 
            '11AE014', '11AE015', '11AE008', '08NK019', '08NP003', '08NK029', '08NP004', '08NK018', '08NG053', 
            '08NK002', '08NJ026', '08NF006', '08NH132', '08NG076', '08NK020', '08NH131', '08NH119', '08NH032', 
            '08NH031', '08NG077', '08NG086', '08NJ009', '08NE010', '08NK021', '08NK005', '08NH005', '08NJ061', 
            '08NH101', '08NH084', '08NJ160', '08NH118', '08NH016', '08NJ013', '08NK022', '08NF001', '08NH021', 
            '08NF005', '08NG012', '08NE074', '08NJ168', '08NH007', '08NH130', '08NP002', '08NK028', '08NG002', 
            '08NK027', '08NG046', '08NP001', '08NH120', '08NG085', '08NG078', '08NF002', '08NG065', '08NH004', 
            '08NH006', '08NJ027', '08NE114', '08NJ167', '08NK016', '08NJ158', '08NN019', '08NL039', '08NB013', 
            '08NE077', '08NL007', '08NN003', '08NM050', '08NB012', '08NE049', '08NE087', '08NL070', '08NM142', 
            '08NN022', '08NM116', '08NL038', '08NL071', '08NN014', '08NM041', '08ND013', '08NB015', '08NM229', 
            '08NE058', '08NL069', '08NM165', '08NM085', '08NM138', '08NN026', '08NL022', '08NL050', '08NE001', 
            '08NM053', '08NC004', '08NM174', '08NM160', '08ND019', '08NM168', '08NN013', '08NL045', '08NA006', 
            '08NM037', '08NA002', '08NM022', '08NM020', '08ND018', '08NM232', '08NE117', '08NN023', '08NN015', 
            '08NE006', '08NE110', '08NM127', '08NL004', '08NA045', '08NB014', '08NM065', '08NM171', '08NM002', 
            '08NM173', '08NB005', '08ND012', '08NN002', '08NN012', '08NM145', '08NB019', '08NA011', '08NE039', 
            '08NM134', '08NE008', '08NN028', '07B055', '10C070', '08MH016', '08PA004', '08MH156', '08MH152', 
            '08MH153', '08MH155', '08HA033', '08MH056', '08HA026', '08MH001', '08MH029', '08HA010', '08HA070', 
            '08MH103', '08C110', '07A090', '05B070', '09A190', '07D050', '10A110', '07D130', '01G070', '04C070']
ec_gages.sort()

# %%
# Subset to just those stations used in the NHM
df2 = df[df.index.isin(ec_gages)]

# %%
df2

# %%
df2[df2['PROV_TERR_STATE_LOC'] == 'QC']

# %%
poiname_list = df2['poi_name'].tolist()

# %%
nc.stringtochar(np.array(poiname_list).astype('S'))

# %%
type(poiname_list)

# %%
poiname_list[0]

# %%
poiname_list[0].encode('latin-1')

# %%
ascii(poiname_list[0])

# %%
import unicodedata

# %%
# Produce ASCII output without multibyte unicode characters
poiname_list_norm = [unicodedata.normalize('NFKD', aa).encode('ascii', 'ignore') for aa in poiname_list]

# unicodedata.normalize('NFKD', poiname_list[0]).encode('ascii', 'ignore')

# %%
nc.stringtochar(np.array(poiname_list_norm).astype('S'))
# poiname_list_norm

# %%

# %%

# %%
df2.to_csv(f'{src_dir}/hydat_sites.tab', sep='\t')

# %%
# df[1] = df[1].apply(lambda x: x + 1)
df2.loc[:, ('drainage_area')] = df2.loc[:, ('drainage_area')].apply(lambda x: x * 0.386102)
df2.loc[:, ('drainage_area_contrib')] = df2.loc[:, ('drainage_area_contrib')].apply(lambda x: x * 0.386102)

# %%
df2.to_csv(f'{src_dir}/hydat_sites.tab', sep='\t')

# %%
