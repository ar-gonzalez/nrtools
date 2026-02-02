from watpy.coredb.coredb import *
from watpy.utils.coreh5 import *
from watpy.wave.wave import *
import glob

db_path = '/data/numrel/DATABASE/CoRe_DB_clone' 


to_correct = ['SLy_BH_m4.8_s0.6--NS_m1.6_s0--d50/bam_128_8_LLF','MS1b_BH_m2.8_s-0.3--NS_m1.6_s0--d40/bam_96_8_LLF','MS1b_BH_m2.8_s-0.6--NS_m1.6_s0--d40/bam_96_8_LLF',
'SLy_BH_m3.2_s0.3--NS_m1.4_s0--d45/bam_96_8_LLF','SLy_BH_m3.2_s0.6--NS_m1.4_s0--d45/bam_96_8_LLF','ALF2_BH_m3.2_s-0.3--NS_m1.6_s0--d40/bam_96_8_LLF',
'SLy_BH_m3.2_s-0.3--NS_m1.4_s0--d45/bam_96_8_LLF','SLy_BH_m3.2_s-0.6--NS_m1.4_s0--d45/bam_96_8_LLF','SLy_BH_m4.8_s0.6--NS_m1.4_s0--d50/bam_128_8_LLF',
'SLy_BH_m4.8_s0.7--NS_m1.4_s0--d50/bam_128_8_LLF','ALF2_BH_m2.7_s0.3--NS_m1.35_s0--d40/bam_96_8_LLF','ALF2_BH_m2.7_s-0.7--NS_m1.35_s0--d40/bam_96_8_LLF',
'ALF2_BH_m2.7_s-0.6--NS_m1.35_s0--d40/bam_96_8_LLF','ALF2_BH_m2.7_s-0.5--NS_m1.35_s0--d40/bam_96_8_LLF','ALF2_BH_m2.7_s-0.3--NS_m1.35_s0--d40/bam_96_8_LLF',
'MS1b_BH_m2.7_s-0.7--NS_m1.35_s0--d40/bam_96_8_LLF','MS1b_BH_m2.7_s-0.6--NS_m1.35_s0--d40/bam_96_8_LLF','MS1b_BH_m2.7_s-0.5--NS_m1.35_s0--d40/bam_96_8_LLF']


def get_dbkey(working_name):
    # input can be string or list of strings
    json_file = '../../../bhns_eobnr/dat/eob_data_may2025.json'
    json_dbkeys = '../../../bhns_eobnr/dat/simnames_dbkeys_all.json'
    with open(json_file, 'r') as jf:
        data = json.load(jf)
        for entry in data:
            try:
                if entry['working_name']==working_name:
                    simname = entry['simulation_name']
                    with open(json_dbkeys, 'r') as jf:
                        data2 = json.load(jf)
                        for entry2 in data2:
                            if entry2['simulation_name']==simname:
                                dbkey = entry2['database_key']
                                break
            except:
                if entry['working_name'] in working_name:
                    simname = entry['simulation_name']
                    with open(json_dbkeys, 'r') as jf:
                        data2 = json.load(jf)
                        for entry2 in data2:
                            if entry2['simulation_name']==simname:
                                dbkey = entry2['database_key']
    return dbkey

name = to_correct[0]

CoM_dir = './'+name.split('/')[0]+'_CoM'
dbk = get_dbkey(name).replace('_',':')
cdb = CoRe_db(db_path)
sim = cdb.sim
data = sim[dbk].run['R01'].data # h5 object
datamd = sim[dbk].run['R01'].md # md object

# get dictionary and dset
datadic = datamd.data
dset = data.read_dset()

# move the CoM corrected ones to the sim path
os.system('cp '+CoM_dir+'/*.txt '+data.path)

new_files = [os.path.split(x)[1] for x in glob.glob('{}/{}'.format(data.path,'Rh_l?_m?_rInf.txt'))]
wm = mwaves(path = data.path, code = 'core', filenames = new_files, 
            mass = datadic['id_mass'], f0 = datadic['id_gw_frequency_Momega22'],
            ignore_negative_m=False)

data.write_EJ_to_txt()
# now add new data to the h5 dset
for (l,m) in wm.modes:
    data.write_strain_to_txt(lm=[(l,m)])
    data.write_psi4_to_txt(lm=[(l,m)])
    dset['rh_{}{}'.format(l,m)] = [os.path.split(x)[1] for x in glob.glob('{}/Rh_l{}_m{}_r*.txt'.format(data.path,l,m))] 

data.create_dset(dset)
data.dump()


# clean up
sim[dbk].run['R01'].clean_txt()
print(data.path)
