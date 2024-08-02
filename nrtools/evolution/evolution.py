import os
from .ev_parfile import *
from .ev_output import *
from ..initialdata.output import *
from ..utils.utils import lin_momentum_from_wvf
from watpy.utils.num import diff1
from watpy.wave.wave import write_headstr, rinf_float_to_str

def ev_folder_parser(folder_name):
    pars = folder_name.split('_')
    resolution = pars[1]
    lmax = pars[2]
    flux = pars[3]
    return resolution, lmax, flux

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
        self.idata = initial_data
        self.check_status()
        self.flux = flux
        self.parfile = Ev_Parameter_File(self.path, self.ev_path, self.idata, resolution, lmax, lmax2, self.flux)
        self.ou = Ev_Output(self.path, self.status, lmax)

        core_out = os.path.join(self.path,'CoReDB')
        try:
            os.mkdir(core_out)
        except FileExistsError:
            print('Directory exists: ',core_out)
        self.core_out = core_out

    def check_status(self):
        ev_log = [i for i in os.listdir(self.path) if i.endswith('out.log')]
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

    def get_grid_spacing_min(self):
        grid_file = [i for i in os.listdir(self.path) if i=='grid_setup.log' or i.startswith('res') or i.startswith('grid')][0]
        with open(os.path.join(self.path,grid_file), 'r') as file:
            for line in file:
                if 'dxyz_finest_BH' in line.strip():
                    dxyz = line.strip().split('= ')[-1].split(' # ')[0]
                    break
        return float(dxyz)
    
    def write_bashfile(self, batchsys='slurm', cluster='ARA'):
        if batchsys=='slurm':
            jobsub = 'run_bam.sh'
            self.write_bashfile_slurm(jobsub, cluster)
        elif batchsys=='pbs':
            jobsub = 'run_bam.pbs'
            self.write_bashfile_pbs(jobsub, cluster)
        else:
            print('ERROR: Unknown batch system. Currently available: slurm, pbs.')

    def run_job(self, batchsys='slurm'):
        if batchsys=='slurm':
            self.run_job_slurm()
        elif batchsys=='pbs':
            self.run_job_pbs()

    def write_bashfile_slurm(self, jobsub = 'run_bam.sh', cluster='ARA'):
        if cluster == 'ARA':
            partition = 'b_standard'
            time = '8-0:00:00'
            memcpu = 'MaxMemPerCPU'
            modules = ['intel/oneapi/2024.0.1','mpi/latest','mpi/openmpi/5.0.2/gcc']
        elif cluster == 'DRACO':
            partition = 'standard' # compute, standard
            time = '3-0:00:00' # inf, 3-0:00:00
            memcpu = '2G'
            modules = ['use.intel-oneapi','mpi/latest']
        elif cluster == 'LRZ':
            partition = 'micro' 
            time = '20:00:00'  
            memcpu = '90G'
            modules = None
        else:
            print('ERROR: Unknown cluster name. Currently available: ARA, DRACO.')

        bss = open(os.path.join(self.path,jobsub), 'a')
        self.bashname = jobsub
        whole_name = self.idata.simname+'-'+self.evname
        bss.write('#!/bin/bash \n')
        bss.write('#SBATCH --partition '+partition+' \n')
        if cluster == 'DRACO':
            bss.write('#SBATCH --qos=multi-node \n')
        bss.write('#SBATCH -J '+whole_name+'\n')
        bss.write('#SBATCH -o bam_out.log \n')
        bss.write('#SBATCH -e error.err \n')
        bss.write('#SBATCH --nodes 4 \n')
        bss.write('#SBATCH --ntasks-per-node=4 \n')
        bss.write('#SBATCH -t '+time+' \n')
        bss.write('#SBATCH --mail-user=alejandra.gonzalez@uni-jena.de \n')
        bss.write('#SBATCH --mail-type=begin \n')
        bss.write('#SBATCH --mail-type=end \n')
        bss.write('#SBATCH --cpus-per-task=6 \n')
        bss.write('#SBATCH  --exclusive \n')
        if cluster == 'LRZ':
            bss.write('#SBATCH --no-requeue \n')
            bss.write('#SBATCH --get-user-env \n')
            bss.write('#SBATCH --account=pn39go \n')
        bss.write('##SBATCH  --mem-per-cpu='+memcpu+' \n\n')
        bss.write('export OMP_NUM_THREADS=6 \n')
        bss.write('export I_MPI_DEBUG=5 \n')
        bss.write('export KMP_AFFINITY=verbose,granularity=fine,scatter \n\n')

        if cluster!='LRZ':
            bss.write('module purge \n')
            for mod in modules:
                if mod == modules[-1]:
                    bss.write('module load '+mod+' \n\n')
                else:
                    bss.write('module load '+mod+' \n')

        exe_file = os.path.join(self.ev_path,'exe/bam_latest')

        bss.write('time mpirun -n 16 '+exe_file+' -nt 6 bam_evo.par \n')
        bss.close()

    def run_job_slurm(self):
        os.chdir(self.path)
        submitjob = 'sbatch ' + os.path.join(self.path,self.bashname)
        os.system(submitjob)

    def write_bashfile_pbs(self, jobsub = 'run_bam.pbs', cluster='HAWK'):
        if cluster == 'HAWK':
            node = 'rome'
            time = '24:00:00'
            modules = ['hlrs-software-stack/.9','hlrs-software-stack/current','impi/2021.9.0','intel/19.1.3','mkl/2023.1.0']
        else:
            print('ERROR: Unknown cluster name. Currently available: HAWK')

        bss = open(os.path.join(self.path,jobsub), 'a')
        self.bashname = jobsub
        bss.write('#!/bin/bash \n')
        bss.write('#PBS -N '+self.evname+'\n')
        bss.write('#PBS -o /lustre/hpe/ws10/ws10.3/ws/xujapigo-bhns/'+self.evname+'/bam_out.log \n')
        bss.write('#PBS -e /lustre/hpe/ws10/ws10.3/ws/xujapigo-bhns/'+self.evname+'/error.err \n')
        bss.write('#PBS -l select=4:node_type='+node+':mpiprocs=32:ompthreads=4 \n')
        bss.write('#PBS -l walltime='+time+' \n')
        bss.write('#PBS -M alejandra.gonzalez@uni-jena.de \n')
        bss.write('export OMP_NUM_THREADS=4 \n')
        bss.write('export I_MPI_DEBUG=5 \n')
        bss.write('export KMP_AFFINITY=verbose,granularity=fine,scatter \n\n')
        bss.write('module purge \n')

        for mod in modules:
            if mod == modules[-1]:
                bss.write('module load '+mod+' \n\n')
            else:
                bss.write('module load '+mod+' \n')

        exe_file = os.path.join(self.ev_path,'exe/bam_latest')

        bss.write('time mpirun -n 128 '+exe_file+' -nt 4 /lustre/hpe/ws10/ws10.3/ws/xujapigo-bhns/'+self.evname+'/bam_evo.par \n')
        bss.close()

    def run_job_pbs(self):
        os.chdir(self.path)
        submitjob = 'qsub ' + os.path.join(self.path,self.bashname)
        os.system(submitjob)

    def get_core_wm_object(self):
        ev_output = self.ou
        id_output = self.idata.ou
            
        _, _, mtot = id_output.get_msun_masses()
        _, Momg22 = id_output.get_gw_freqs()
        f0 = Momg22 / (2*np.pi) / mtot 
        dfiles = [os.path.split(x)[1] for x in glob.glob('{}/{}'.format(ev_output.out_inv_dir,'Rpsi4mode??_r*.l0'))]
        wm = mwaves(path = ev_output.out_inv_dir, code = 'bam', filenames = dfiles, 
            mass = mtot, f0 = f0,
            ignore_negative_m=True)
        return wm

    def get_core_data(self):
        id_output = self.idata.ou
            
        mbh, mns, _ = id_output.get_msun_masses()
        wm = self.get_core_wm_object()

        core_out = self.core_out
            
        # Get waveforms
        for r in wm.radii:
            for (l,m) in wm.modes:
                psilm = wm.get(var='Psi4',l=l, m=m, r=r)
                psilm.write_to_txt('Psi4', core_out)
                hlm = wm.get(l=l, m=m)
                hlm.write_to_txt('h', core_out)

        # Get energetics
        madm, jadm = id_output.get_ADM_qtys()
        wm.energetics(mbh, mns, madm, jadm, path_out = core_out)

    def get_lin_momentum(self):
        wm = self.get_core_wm_object()
        h     = {}
        h_dot = {}
        u     = {}

        for rad in wm.radii:
            for lm in wm.modes: 
                w         = wm.get(l=lm[0], m=lm[1], r=rad)
                t         = w.time
                u[lm]     = w.time_ret()
                h[lm]     = w.h
                h_dot[lm] = diff1(t, h[lm])
            
            Px, Py, Pz, P = lin_momentum_from_wvf(h, h_dot, t, u, wm.modes)
            headstr  = write_headstr(rad,wm.mass)
            headstr += "Px:0 Py:1 Pz:2 P:3 t:4 u:5"
            data = np.c_[Px, Py, Pz, P, t, w.time_ret()]
            rad_str = rinf_float_to_str(rad)
            fname = "P_r"+rad_str+".txt"
            np.savetxt('{}/{}'.format(self.core_out,fname), data, header=headstr)


    def output_metadata(self, author='alejandra'):
        '''
        Outputs a dictionary of the simulation's important parameters
        and print them in a file `metadata.txt`
        '''
        simkeys = {}
        simname = "/".join([self.idata.simname,self.evname])
        simkeys['simulation_name'] = simname

        if author=='alejandra':
            simkeys['author'] = "Alejandra Gonzalez <alejandra.gonzalez@uni-jena.de>"
        elif author=='fabio':
            simkeys['author'] = "Fabio Magistrelli <fabio.magistrelli@uni-jena.de>"
        elif author=='francesco':
            simkeys['author'] = "Francesco Brandoli <francesco.brandoli@uni-jena.de>"
        else:
            print("Who is ",author,"?")

        # Initial Data
        simkeys['ID_type'] = "BHNS"
        simkeys['ID_code'] = "Elliptica"
        simkeys['initial_separation'] = self.idata.simname.split('_')[-1][-2:]
        #simkeys['initial_orbital_angular_velocity'] = self.idata.ou.id_dic['BHNS_angular_velocity']
        fhz, fm22 = self.idata.ou.get_gw_freqs()
        simkeys['initial_gw_frequency_Hz'] = "{:.2f}".format(fhz)
        simkeys['initial_gw_frequency_Momega22'] = "{:.5f}".format(fm22)
        simkeys['EoS'] = self.idata.simname.split('_')[0]
        m1, m2, _ = self.idata.ou.get_msun_masses()
        simkeys['initial_mass1'] = "{:.2f}".format(m1) #self.idata.simname.split('_')[2][1:]
        simkeys['initial_mass2'] = "{:.2f}".format(m2) #self.idata.simname.split('_')[4][1:]

        # Evolution Data
        try:
            self.get_core_data()
            h22_file = sorted([i for i in os.listdir(self.core_out) if i.startswith('Rh_l2_m2_r')])[-1]
            time22, amp22 = np.loadtxt(fname=os.path.join(self.core_out,h22_file), comments='#', usecols=(0,4), unpack=True)
            tmrg = time22[np.argmax(amp22)]
        
            _, _, px2, _, tpx2 = self.ou.extract_objects_tracks() # of the BH
            until_merger = np.where((tpx2<(tmrg+1))&(tpx2>(tmrg-1)))[0]
            until_merger = until_merger[0]
            px2_um = px2[:until_merger]
            orbits = list(px2_um).count(px2_um[0]) - 1 # minus the starting point
        except:
            tmrg = 0
            orbits = 0

        simkeys['evolution_code'] = "BAM"
        simkeys['grid_spacing_min'] = self.get_grid_spacing_min()
        simkeys['merger_time'] = "{:.2f}".format(tmrg)
        simkeys['number_of_orbits'] = orbits

        # Remnant and PM
        t, mbh, sx, sy, sz, s = self.ou.final_bh_properties()
        simkeys['final_time'] = "{:.2f}".format(t)
        simkeys['remnant_mass'] = "{:.4f}".format(mbh)
        simkeys['remnant_dimensionless_spin'] = ", ".join([str(sx/(mbh*mbh)), str(sy/(mbh*mbh)), str(sz/(mbh*mbh))])

        # Print
        with open(os.path.join(self.path,'metadata.txt'), 'w') as f:
            f.write("#-----------------------------------------------\n")
            f.write("# General information about the simulation\n")
            f.write("#-----------------------------------------------\n")
            for key, value in simkeys.items():
                if key=='ID_type':
                    f.write("\n#-----------------------------------------------\n")
                    f.write("# Initial Data\n")
                    f.write("#-----------------------------------------------\n")
                    f.write('%s =   %s\n' % (key, value))
                elif key=='evolution_code':
                    f.write("\n#-----------------------------------------------\n")
                    f.write("# Evolution\n")
                    f.write("#-----------------------------------------------\n")
                    f.write('%s =   %s\n' % (key, value))
                elif key=='final_time':
                    f.write("\n#-----------------------------------------------\n")
                    f.write("# Remnant Properties\n")
                    f.write("#-----------------------------------------------\n")
                    f.write('%s =   %s\n' % (key, value))
                else:
                   f.write('%s =   %s\n' % (key, value)) 

        return simkeys



