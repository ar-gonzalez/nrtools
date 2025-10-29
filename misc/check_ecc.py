import json
import os
from nrtools.initialdata.initialdata import *

import re

def extract_last_ecc_from_file(file_path):
    last_ecc = None  # Initialize to store the last ecc value
    with open(file_path, 'r') as file:
        # Read all lines
        lines = file.readlines()
        
        # Iterate through all lines
        for line in lines:
            # Search for the line containing 'ecc' and extract the value
            match = re.search(r'#\s*ecc\s*=\s*([0-9e\.\-]+)', line)
            if match:
                # Update last_ecc each time we find a match
                last_ecc = float(match.group(1))
                
    return last_ecc  # Return the last ecc value found, or None if not found

basedir = '/home/mo63poy/BHNS_INITIAL_DATA'
id_exe = '/home/mo63poy/Elliptica/Exe/elliptica'

file1 = '../dat/simnames_dbkeys.json'  # dic with dbkeys
file2 = '../dat/eob_data_2023.json'  # json with all data

# Read the first JSON file
with open(file1, 'r') as f:
    simnames = json.load(f)
    
# Read the second JSON file
with open(file2, 'r') as f:
    worknames = json.load(f)

    
# Create a lookup dictionary for fast access to database_key by simulation_name
key_dict = {entry['simulation_name']: entry['database_key'] for entry in simnames}



mylist = []
runerror = []
# Update list2 with the database_key if simulation_name matches
for entry in worknames:
    sim_name = entry.get('simulation_name')
    if sim_name in key_dict:
        working_name = entry.get('working_name')
        name = working_name.split('/')[0]
        simpath = os.path.join(basedir,name)
        evopath = os.path.join(basedir,working_name)

        bhns_id = Initial_Data(path=simpath,params=None, id_exe=id_exe)
        id_output = bhns_id.ou
        id_output.read_md()
        dic = id_output.id_dic

        omega = dic['BHNS_angular_velocity']
        mass = dic['BHNS_ADM_mass']
        movpunc = os.path.join(evopath,'bam_evo/moving_puncture_distance.lxyz6')        
        comm = 'python EccRed_noForceBal.py --Mass '+mass+' --Omega '+omega+' --tskip 50 --dmin 40 '+movpunc+' > fit.data'
        try:
            os.system(comm)
            ecc_value = extract_last_ecc_from_file('fit.data')
            #print(working_name,ecc_value)
            '''
            if ecc_value>0.02 or ecc_value==None:
                continue
            else:
                mylist.append((working_name,ecc_value))
            '''
            mylist.append((working_name,ecc_value))
        except TypeError:
            runerror.append(working_name)
print(runerror)

for entry in mylist:
    print(entry)
