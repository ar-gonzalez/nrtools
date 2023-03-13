import os

########################################
# ID Metadata class
########################################

class Metadata():
    """
    Reads the produced `BHNS_properties.txt` file and
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
        if status!='Done':
            print('===> Error: Class available only for finished ID')
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

    def read_md(self):
        # Reads md into a dictionary of the highest resolution
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
        # Reads md into a dictionary from the second highest resolution
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
        
