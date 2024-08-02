import os
import numpy as np
from ..utils.utils import get_id_gw_frequency_Hz, get_id_gw_frequency_Momega22, get_kappa2t

########################################
# ID Output class
########################################

class Output():
    """
    Reads the produced `BHNS_properties.txt` files and
    produces metadata compatible with WATPy
    
    ------------------
    Initialization:
    ------------------
    simname : name of the simulation
    path    : where the initial data should be produced
    status  : status of the initial data run
    id_outdir  : path to the output directory
    """
    def __init__(self, simname, path, status, id_outdir, verbose=False):
        if status=='Not started':
            print('===> Error: ID has not started yet')
        else:
            self.path = path
            self.simname = simname
            self.id_outdir = id_outdir
            self.verbose = verbose
            if status=='Done':
                highr = 0
                for ir,reso in enumerate(sorted(os.listdir(self.id_outdir))):
                    rr = float(reso.split('_')[7].split('x')[0])
                    if float(rr)>highr:
                        highr = rr
                        ihighr = ir
                    elif rr==highr:
                        if reso.endswith('_01'):
                            highr = rr
                            ihighr = ir
                self.hr_path = sorted(os.listdir(self.id_outdir))[ihighr] # Path to highest resolution ID
                self.sr_path = sorted(os.listdir(self.id_outdir))[ihighr-1] # Path to second highest resolution ID
                self.outpath = os.path.join(self.path,self.simname+'_00')
                self.txt_path = os.path.join(os.path.join(self.outpath,self.hr_path),'BHNS_properties.txt')
                self.read_md()

    def get_resolutions(self):
        '''
        Returns arrays of: 
        - resolution folder names (string)
        - the available resolutions (integer)
        '''
        resolutions = []
        for res_folder in sorted(os.listdir(self.outpath)):
            dims = res_folder.split('_')[-2]
            resolutions.append(int(dims.split('x')[0]))
        return sorted(os.listdir(self.outpath)), resolutions

    def read_md(self):
        '''
        Reads metadata into a dictionary of the highest resolution
        '''
        dic = {}
        with open(self.txt_path) as file:
            dic['simname'] = self.simname
            for i,line in enumerate(file):
                if i<3 or len(line.strip()) == 0:
                    continue
                else:
                    k, v = line.strip().split('=')
                    dic[k.strip()] = v.strip()
        file.close()
        self.id_dic = dic
        if self.verbose:
            print('==> Metadata: ',dic)

    def read_md_sr(self):
        '''
        Reads metadata into a dictionary from the second highest resolution
        '''
        dic = {}
        txt_sr_path = os.path.join(os.path.join(self.outpath,self.sr_path),'BHNS_properties.txt')
        with open(txt_sr_path) as file:
            dic['simname'] = self.simname
            for i,line in enumerate(file):
                if i<3 or len(line.strip()) == 0:
                    continue
                else:
                    k, v = line.strip().split('=')
                    dic[k.strip()] = v.strip()
        file.close()
        self.id_dic2 = dic
        
    def get_value_from_resolutions(self,value):
        '''
        Input: Parameter value (see e.g. parfile.py)
        Returns: Array of that value recovered in available resolutions
        '''
        res_dirs = os.listdir(self.outpath)
        res_dirs.sort()
        value_array = []
        for res in res_dirs:
            res_file = os.path.join(os.path.join(self.outpath,res),'BHNS_properties.txt')
            with open(res_file) as file:
                for i,line in enumerate(file):
                    if i<3 or len(line.strip()) == 0:
                        continue
                    else:
                        k, v = line.strip().split('=')
                        if k.strip() == value:
                            value_array.append(float(v.strip()))
            file.close()
        return value_array
    
    def get_gw_freqs(self):
        '''
        Returns ID GW frequencies for CoRe:
        id_gw_frequency_Hz: Initial data: initial GW frequency (Hz)
        id_gw_frequency_Momega22: Initial data: Mass-rescaled initial GW frequency (c=G=Msun=1 units)
        '''
        _, _, mtot = self.get_msun_masses()
        omega = float(self.id_dic['BHNS_angular_velocity'])
        id_gw_frequency_Hz = get_id_gw_frequency_Hz(omega)
        id_gw_frequency_Momega22 = get_id_gw_frequency_Momega22(omega, mtot)
        return id_gw_frequency_Hz, id_gw_frequency_Momega22
    
    def get_msun_masses(self):
        '''
        Outputs the grav mass of
        M_BH: BH's mass
        M_NS: NS' mass
        Mtot: total binary mass
        For a NS with M_b = 1.6
        '''
        '''
        eos = self.id_dic['NS_EoS_description']
        mb = float(self.id_dic['NS_baryonic_mass_current'])
        if eos=='SLy':
            if mb>1.5: # for mb=1.6
                M_NS = 1.4199062240028333
            else: # for mb=1.4
                M_NS = 1.259943000559463
        elif eos=='MS1b':
            if mb>1.5: # for mb=1.6
                M_NS = 1.4611674103103354
            else: # for mb=1.4
                M_NS = 0
        elif eos=='ALF2':
            if mb>1.5: # for mb=1.6
                M_NS = 1.441691281566144
            else: # for mb=1.4
                M_NS = 0
        else:
            M_NS = 0
            print("===> Error: EoS not recognized, please add to /nrtools/initialdata/output.py")
        '''
        M_BH = float(self.id_dic['BH_Christodoulou_mass_current'])
        M_NS = float(self.id_dic['NS_TOV_ADM_mass'])
        return M_BH, M_NS, M_BH + M_NS

    def get_tidal_params(self):
        eos = self.id_dic['NS_EoS_description']
        mb = float(self.id_dic['NS_baryonic_mass_current'])
        mbh, mg, _ = self.get_msun_masses()
        if eos=='SLy':
            if mb>1.5: # for mb=1.6
                lam = 276.48036295542073
                kappa = get_kappa2t(lam,mbh,mg)
            else: # for mb=1.4
                lam = 599.2794774456522
                kappa = get_kappa2t(lam,mbh,mg)
        elif eos=='MS1b':
            if mb>1.5: # for mb=1.6
                lam = 1006.2885125921738
                kappa = get_kappa2t(lam,mbh,mg)
            elif mb>1.39: # for mb=1.4
                lam = 1931.501612858557
                kappa = get_kappa2t(lam,mbh,mg)
            elif mb<1.36: # for mb=1.35
                lam = 2312.470950148368
                kappa = get_kappa2t(lam,mbh,mg)
        elif eos=='ALF2':
            if mb>1.5: # for mb=1.6
                lam = 493.7973627693258
                kappa = get_kappa2t(lam,mbh,mg)
            elif mb>1.39: # for mb=1.4
                lam = 1006.3969784624323
                kappa = 0
            elif mb<1.36: # for mb=1.35
                lam = 1206.8906420652897
                kappa = get_kappa2t(lam,mbh,mg)
        else:
            lam = 0
            kappa = 0
            print("===> Error: EoS not recognized, please add to /nrtools/initialdata/output.py")
        
        return lam, kappa

    def get_ADM_qtys(self):
        '''
        Outputs
        M_ADM: binary ADM mass
        J_ADM: binary ADM angular momentum
        '''
        M_ADM = float(self.id_dic['BHNS_ADM_mass'])
        jx = float(self.id_dic['BHNS_Jx_ADM'])
        jy = float(self.id_dic['BHNS_Jy_ADM'])
        jz = float(self.id_dic['BHNS_Jz_ADM'])
        J_ADM = np.sqrt(jx*jx + jy*jy + jz*jz)
        return M_ADM, J_ADM
    
