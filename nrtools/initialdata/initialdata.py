import numpy as np
import os
from .parfile import *
from .output import *
from ..evolution.evolution import *
import matplotlib.pyplot as plt

########################################
# Main class for Initial Data
########################################

class Initial_Data():
    """
    ------------------
    Initialization:
    ------------------
    path    : where the initial data should be produced
    params  : dictionary with the basic parameters of the binary
    id_exe  : path/to/initialdata/executable
    """
    def __init__(self, path, params, id_exe):
        self.path = path
        self.user_params = params
        self.parfile = Parameter_File(self.path, self.user_params)

        if params==None and id_exe==None:
            print("==> I will assume ID exists already in ", self.path)
            self.simname = self.path.split('/')[-1]
            self.check_status()  
            self.ou = Output(self.simname, self.path, self.status, self.id_outdir)
        else:
            print("==> Create new ID")
            self.id_exe = id_exe
            self.make_parfile()  
            self.write_bashfile()

    def check_status(self):
        self.id_outdir = os.path.join(self.path,self.simname+'_00')
        try:
            num_res = os.listdir(self.id_outdir)
            self.status = 'Ongoing'
        except FileNotFoundError:
            self.status = 'Not started'
        
        if self.status=='Ongoing':
            log_file = [i for i in os.listdir(self.path) if i.endswith('.log')][0]
            with open(os.path.join(self.path,log_file),'r') as hrf:
                cont = hrf.read()
                fertig = "} construct_initial_data :))"
                if fertig in cont:
                    self.status = 'Done'
        print("==> Initial Data Status: ",self.status)

    def make_parfile(self):
        pardic = self.parfile.pardic
        mm = float(pardic['BH_irreducible_mass'])
        cc = pardic['BH_chi_z']
        mbNS = pardic['NS_baryonic_mass']
        sep = pardic['BHNS_separation']

        self.simname = pardic['NS_EoS_description'] + '_BH_m' + str(round(mm,1)) + '_s' + str(cc) + '--NS_m' + str(mbNS) + '_s0--d' + str(sep)
        self.simpath = os.path.join(self.path,self.simname)

        # Create directory
        try:
            os.mkdir(self.simpath)
        except FileExistsError:
            print('Directory exists already: ',self.simpath)

        # Write par file
        with open(os.path.join(self.simpath,self.simname+'.par'), 'w') as f:
            for key, value in pardic.items():
                f.write('%s =   %s\n' % (key, value))

    def write_bashfile(self,bashname= 'run_elliptica.sh'):
        bss = open(os.path.join(self.simpath,bashname), 'a')
        self.bashname = bashname
        bss.write('#!/bin/bash \n')
        bss.write('#SBATCH --partition s_standard \n')
        bss.write('#SBATCH -J '+self.simname+'\n')
        bss.write('#SBATCH -o '+os.path.join(self.simpath,'out.log')+' \n')
        bss.write('#SBATCH -N 1 \n')
        bss.write('#SBATCH -n 1 \n')
        bss.write('#SBATCH -t 8-8:00:00 \n')
        bss.write('#SBATCH --mail-user=alejandra.gonzalez@uni-jena.de \n')
        bss.write('#SBATCH --mail-type=begin \n')
        bss.write('#SBATCH --mail-type=end \n')
        bss.write('#SBATCH --cpus-per-task=32 \n')
        bss.write('#SBATCH  --mem-per-cpu=MaxMemPerCPU \n\n')
        bss.write('export OUTDIR='+self.simpath+' \n')
        bss.write('export PAR='+self.simname+' \n')
        bss.write('export ELLIPTICA='+self.id_exe+' \n')
        bss.write('export OMP_NUM_THREADS=16 \n\n')
        bss.write('module purge \n')
        bss.write('module load compiler/gcc/10.2.0 \n')
        bss.write('module load compiler/intel/2020-Update2 \n\n')
        bss.write('time srun $ELLIPTICA $OUTDIR/$PAR.par > $OUTDIR/job.log \n')
        bss.close()

    def run_job(self):
        submitjob = 'sbatch ' + os.path.join(self.simpath,self.bashname)
        os.system(submitjob)

    def check_accuracy(self):
        # Expected values read from parfile:
        pardic = self.parfile.pardic
        bhmass_expected = float(pardic['BH_irreducible_mass'])
        bhchiz_expected = float(pardic['BH_chi_z'])
        nsmass_expected = float(pardic['NS_baryonic_mass'])
        # Current obtained values:
        dic = self.ou.id_dic
        bhmass_current = float(dic['BH_irreducible_mass_current'])
        bhchiz_current = float(dic['BH_chi_z_current'])
        nsmass_current = float(dic['NS_baryonic_mass_current'])
        # Get Errors:
        bhmass_error = np.abs(bhmass_current-bhmass_expected)
        bhchiz_error = np.abs(bhchiz_current-bhchiz_expected)
        nsmass_error = np.abs(nsmass_current-nsmass_expected)
        # Print status
        print("==> Checking Accuracy: \n")
        print("~ BH_irreducible_mass: expected: ", bhmass_expected, "| current: ", bhmass_current,"| Error= ",bhmass_error,"(",bhmass_error*100/bhmass_expected,"%) \n")
        print("~ BH_chi_z: expected: ", bhchiz_expected, "| current: ", bhmass_current,"| Error= ",bhchiz_error,"(",bhchiz_error*100/bhchiz_expected,"%) \n")
        print("~ NS_baryonic_mass: expected: ", nsmass_expected, "| current: ", bhmass_current,"| Error= ",nsmass_error,"(",nsmass_error*100/nsmass_expected,"%)\n")
        return bhmass_error, bhchiz_error, nsmass_error

    def check_convergence(self, patch = 'right_BH_around_front'):
        res_folders, resolutions = self.ou.get_resolutions()

        # Get ADM mass
        madm = self.ou.get_value_from_resolutions('BHNS_ADM_mass')
        plt.title(self.simname)
        plt.plot(resolutions,madm,c='#1b9e77',marker='*')
        plt.grid()
        plt.xlabel(r'Resolution N')
        plt.ylabel(r'$M_{\rm ADM}$')
        plt.show()

        # Get Hamiltonian and Momentum constraints
        plt.title(self.simname)
        hami = []
        momx = []
        momy = []
        momz = []
        for i, res in enumerate(res_folders):
            diagnostics_dir = os.path.join(os.path.join(self.ou.outpath,res),'Diagnostics_00')
            file_path = os.path.join(diagnostics_dir,patch + '_X0Y0Z0_0d.txt')
            ham, mx, my, mz = np.loadtxt(file_path,comments='#',usecols=(3,9,6,12),unpack=True)
            hami.append(ham[-1])
            momx.append(mx[-1])
            momy.append(my[-1])
            momz.append(mz[-1])
        
        plt.plot(resolutions,hami,label='H',c='#1b9e77',marker='*')
        plt.plot(resolutions,momx,label=r'$M^x$',c='#d95f02',marker='*')
        plt.plot(resolutions,momy,label=r'$M^y$',c='#7570b3',marker='*')
        plt.plot(resolutions,momz,label=r'$M^z$',c='#e7298a',marker='*')

        plt.grid()
        plt.yscale("log")
        plt.xlabel(r'Resolution N')
        plt.ylabel(r'$L^2$ norm')
        plt.legend()
        plt.show()

    def evolve(self, ev_path, resolution=64, lmax=10, lmax2=6, flux='LLF'):
        '''
        Input: path/to/evolution/code, resolution, refinement levels for BH,
                refinement levels for NS, flux reconstruction scheme
        Returns: initializes Evolution object
        '''
        if self.status=='Done':
            ev_name = "bam_"+str(resolution)+"_"+str(lmax)+'_'+flux
            path = os.path.join(self.path,ev_name)
            try:
                os.mkdir(path)
            except FileExistsError:
                print('Directory exists: ',path)
            evolution = Evolution(path, ev_path, self.ou, resolution, lmax, lmax2, flux)
        else:
            print("Error: Initial data is not finished!")
            evolution = None
        return evolution
