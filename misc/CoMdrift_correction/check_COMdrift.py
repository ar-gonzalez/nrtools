import numpy as np
from nrtools.initialdata.initialdata import *
from nrtools.evolution.evolution import *
from nrtools.utils.utils import windowing
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import EOBRun_module
from watpy.wave.gwutils import *
from watpy.wave.wave import *
from watpy.utils.units import MSun_meter
import matplotlib
from bajes.obs.gw.noise import get_design_sensitivity
matplotlib.rcParams['text.usetex']= True
matplotlib.rcParams['font.serif']= 'Palatino' 
matplotlib.rcParams['font.size']= 15 #28
import argparse
import os
import json

def modes_to_k(modes):
    """
    Map multipolar (l,m) -> linear index k
    """
    return [int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes]

def search_json(obj, search_string):
    if isinstance(obj, dict):
        for key, value in obj.items():
            if search_string in str(key) or search_string in str(value):
                return True
            if search_json(value, search_string):
                return True
    elif isinstance(obj, list):
        for item in obj:
            if search_json(item, search_string):
                return True
    elif search_string in str(obj):
        return True
    return False


not_done = ['MS1b_BH_m2.7_s0.6--NS_m1.35_s0--d35/bam_96_8_LLF','MS1b_BH_m2.7_s0.3--NS_m1.35_s0--d35/bam_96_8_LLF','ALF2_BH_m2.8_s-0.6--NS_m1.4_s0--d35/bam_96_8_LLF','ALF2_BH_m2.8_s-0.6--NS_m1.4_s0--d35/bam_96_8_LLF','ALF2_BH_m2.8_s-0.5--NS_m1.4_s0--d35/bam_96_8_LLF','ALF2_BH_m3.2_s0.6--NS_m1.6_s0--d40/bam_96_8_LLF','ALF2_BH_m2.7_s0.7--NS_m1.35_s0--d35/bam_96_8_LLF']


if __name__ == '__main__':

    basedir = '/home/mo63poy/BHNS_INITIAL_DATA'
    id_exe = '/home/mo63poy/Elliptica/Exe/elliptica'
    ev_path = '/home/mo63poy/BAM'
    Msuns  = 4.925491025543575903411922162094833998e-6
    Mpc_m  = 3.085677581491367278913937957796471611e22
    Msun_m = 1.476625061404649406193430731479084713e3
    PC_SI  = 3.085677581491367e+16 # m
    MPC_SI = 1e6 * PC_SI

    
    data_file = "../dat/eob_data_latest.json"
    core_release = "../dat/simnames_dbkeys.json"

    all_dics = json.load(open(data_file))
    core_dics = json.load(open(core_release))

    dics = [dic for dic in all_dics]# if dic['working_name'] not in not_done]

    q_tp = []
    c_tp = []
    umm = []
    for dic in dics:
        nname = dic['working_name']
        print(nname)
        if search_json(core_dics,dic['simulation_name']):
            print('>> IN CURRENT RELEASE')
        else:
            print('>> NOT IN THE RELEASE, CHECK!!!')
        simname = nname.split('/')[0]
        evo = nname.split('/')[1]
        resolution, lmax, flux = ev_folder_parser(evo)
     
        simpath = os.path.join(basedir,simname)
        bhns_id = Initial_Data(path=simpath,params=None, id_exe=id_exe)
        bhns_ev = Evolution(path=os.path.join(simpath,evo), ev_path=ev_path, initial_data=bhns_id, resolution=int(resolution), lmax=int(lmax), lmax2=6,flux=flux)
        ev_output = bhns_ev.ou
        ev_output.plot_moving_puncture()


        