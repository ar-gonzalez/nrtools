import os
from .parfile import *

########################################
# Main classes for Initial Data
########################################

class Initial_Data():
    """
    ------------------
    Initialization:
    ------------------
    path    : where the initial data should be produced
    params  : dictionary with the basic parameters of the binary
    """
    def __init__(self, path='.', params, id_exe='$HOME/BHNS_Elliptica/Elliptica/Exe/elliptica'):
        self.path = path

        if params==None:
            print("Need specific parameters of binary to simulate")
        
        self.id_exe = id_exe
        self.make_parfile(params)

        # Create directory
        simpath = os.path.join(self.path,simname)
        try:
            os.mkdir(simpath)
        except FileExistsError:
            print('Directory exists already: ',simname)

        self.simpath = simpath
        self.write_bashfile()
        self.run_job()

    def make_parfile(self, params):
        parfile = Parameter_File(self.path, params)
        pardic = parfile.pardic

        mm = float(pardic['BH_irreducible_mass'])
        cc = pardic['BH_chi_z']
        mbNS = pardic['NS_baryonic_mass']
        sep = pardic['BHNS_separation']

        self.simname = pardic['NS_EoS_description'] + '_BH_m' + str(round(mm,1)) + '_s' + str(cc) + '--NS_m' + str(mbNS) + '_s0--d' + str(sep)

        # Write par file
        with open(self.path+self.simname+'/'+self.simname+'.par', 'w') as f:
            for key, value in pardic.items():
                f.write('%s =   %s\n' % (key, value))

    def write_bashfile(self,bashname= 'run_elliptica.sh'):
        bss = open(os.path.join(self.simpath,bashname), 'a')
        self.bashname = bashname
        bss.write('#!/bin/bash \n')
        bss.write('#SBATCH --partition s_standard \n')
        bss.write('#SBATCH -J '+simname+'\n')
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