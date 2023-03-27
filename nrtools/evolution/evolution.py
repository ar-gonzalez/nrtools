import os
from .ev_parfile import *

########################################
# Main class for Initial Data
########################################

class Evolution():
    """
    ------------------
    Initialization:
    ------------------
    path          : where the evolution should be produced
    ev_path       : path/to/evolution/code/
    initial_data  : initial data object
    resolution    : resolution for the simulation
    lmax, lmax2   : refinement levels for obj1 and obj2
    flux          : flux reconstruction scheme (LLF or EFL)
    """
    def __init__(self, path, ev_path, initial_data, resolution, lmax, lmax2, flux):
        self.path = path
        self.evname = self.path.split('/')[-1]
        self.ev_path = ev_path
        self.check_status()
        self.flux = flux
        self.parfile = Ev_Parameter_File(self.path, self.ev_path, initial_data, resolution, lmax, lmax2, self.flux)

    def check_status(self):
        ev_log = [i for i in os.listdir(self.path) if i.endswith('.log')]
        if len(ev_log)==0:
            self.status = 'Not started'
        else:
            self.status = 'Ongoing'
            log_file = [i for i in ev_log if i.startswith('bam')][0]
            with open(os.path.join(self.path,log_file),'r') as hrf:
                    cont = hrf.read()
                    fertig = "Thank you for running b a m."
                    if fertig in cont:
                        self.status = 'Done'
        
        print("==> Evolution Status: ",self.status)

    def write_bashfile(self, bashname = 'run_bam.sh', cluster='ARA'):
        if cluster == 'ARA':
            partition = 'b_standard'
            time = '8-0:00:00'
            memcpu = 'MaxMemPerCPU'
            modules = ['mpi/openmpi/2.1.3-gcc-7.3.0','libs/fftw/3.3.7-gcc-7.3.0','apps/mathematica/12.0','mpi/intel/2019-Update3','compiler/intel/2019-Update3']
        elif cluster == 'DRACO':
            partition = 'long' # compute, standard
            time = '14-0:00:00' # inf, 3-0:00:00
            memcpu = '2G'
            modules = ['icc/latest','mkl/latest','compiler/gcc/10.2.1','mpi/openmpi/4.1.1']
        else:
            print('ERROR: Unknown cluster name. Currently available: ARA, DRACO.')

        bss = open(os.path.join(self.path,bashname), 'a')
        self.bashname = bashname
        bss.write('#!/bin/bash \n')
        bss.write('#SBATCH --partition '+partition+' \n')
        bss.write('#SBATCH -J '+self.evname+'\n')
        bss.write('#SBATCH -o bam_out.log \n')
        bss.write('#SBATCH -e error.err \n')
        bss.write('#SBATCH --nodes 2 \n')
        bss.write('#SBATCH --ntasks-per-node=4 \n')
        bss.write('#SBATCH -t '+time+' \n')
        bss.write('#SBATCH --mail-user=alejandra.gonzalez@uni-jena.de \n')
        bss.write('#SBATCH --mail-type=begin \n')
        bss.write('#SBATCH --mail-type=end \n')
        bss.write('#SBATCH --cpus-per-task=6 \n')
        bss.write('#SBATCH  --exclusive \n')
        bss.write('#SBATCH  --mem-per-cpu='+memcpu+' \n\n')
        bss.write('export OMP_NUM_THREADS=6 \n')
        bss.write('export I_MPI_DEBUG=5 \n')
        bss.write('export KMP_AFFINITY=verbose,granularity=fine,scatter \n\n')
        bss.write('module purge \n')

        for mod in modules:
            if mod == modules[-1]:
                bss.write('module load '+mod+' \n\n')
            else:
                bss.write('module load '+mod+' \n')
                
        if self.flux=='EFL':
            exe_file = os.path.join(self.ev_path,'exe/bam_wEFL')
        else:
            exe_file = os.path.join(self.ev_path,'exe/bam_noEFL')

        bss.write('time mpirun '+exe_file+' -nt 6 bam_evo.par \n')
        bss.close()

    def run_job(self):
        os.chdir(self.path)
        submitjob = 'sbatch ' + os.path.join(self.path,self.bashname)
        os.system(submitjob)
