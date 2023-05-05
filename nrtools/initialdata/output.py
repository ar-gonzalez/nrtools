import os
import numpy as np

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
            self.hr_path = sorted(os.listdir(self.id_outdir))[-1] # Path to highest resolution ID
            self.sr_path = sorted(os.listdir(self.id_outdir))[-2] # Path to second highest resolution ID
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
    
    def get_mtot_msun(self):
        '''
        Outputs the total grav mass of the binary
        '''
        return self.id_dic['BH_irreducible_mass_current'] + self.id_dic['NS_TOV_ADM_mass']
    
    def get_gw_freq(self):
        '''
        Gets the GW initial frequency in mass rescaled units
        Tip: Use watpy/utils/units.py to convert frequency from
        geometric units to SI (Herz)
        '''
        mtot = self.M
        omega = self.id_dic['BHNS_angular_velocity']
        return mtot*omega / (2*np.pi)