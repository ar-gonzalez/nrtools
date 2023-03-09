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
    path    : where the initial data should be produced
    params  : dictionary with the basic parameters of the binary
    """
    def __init__(self, simname, path='.', status='Done', id_outdir):
        if status!='Done':
            print('===> Error: Class available only for finished ID')
        else:
            self.path = path
            self.simname = simname
            self.id_outdir = id_outdir
            self.hr_path = sort(os.listdir(self.id_outdir))[-1] # Path to highest resolution ID
            self.sr_path = sort(os.listdir(self.id_outdir))[-2] # Path to second highest resolution ID
            self.txt_path = os.path.join(self.hr_path,'BHNS_properties.txt')
            self.read_md()

    def read_md(self):
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
        
