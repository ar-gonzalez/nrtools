import os
import matplotlib.pyplot as plt
import numpy as np

########################################
# EV Output class
########################################

class Ev_Output():
    """
    Reads the produced output files and
    plots current results
    
    ------------------
    Initialization:
    ------------------
    simname : name of the simulation
    path    : evolution main dir
    status  : status of the ievolution run
    lmax    : refinement levels for obj1
    """
    def __init__(self, path, status, lmax):
        if status=='Not started':
            print('===> Error: Evolution has not started yet')
        else:
            folder_name = [i for i in os.listdir(path) if os.path.isdir(os.path.join(path,i)) and i.startswith('bam')][0]
            self.outpath = os.path.join(path,folder_name) # where the sim output is
            self.plotsdir = os.path.join(path,'plots')
            try:
                os.mkdir(self.plotsdir)
            except FileExistsError:
                print('Directory exists: ',self.plotsdir)

            self.lmax = lmax

            # Output directories:
            self.out_0d_dir = os.path.join(self.outpath, 'output_0d')
            self.out_1d_dir = os.path.join(self.outpath, 'output_1d')
            self.out_2d_dir = os.path.join(self.outpath, 'output_2d')
            self.out_inv_dir = os.path.join(self.outpath, 'Invariants')

    def plot_moving_puncture(self):
        try:
            mp_file = [i for i in os.listdir(self.outpath) if i.startswith('moving_puncture_distance.lxyz')][0]
            self.mp_path = os.path.join(self.outpath,mp_file)

            px_ns, py_ns, px_bh, py_bh = np.loadtxt(fname=self.mp_path, comments='"', usecols=(0,1,3,4), unpack=True)

            plt.scatter(px_ns,py_ns,color='#7fbf7b',label="NS",marker='.')
            plt.scatter(px_bh,py_bh,color='#af8dc3',label="BH",marker='.')
            plt.grid()
            plt.legend()
            plt.savefig(os.path.join(self.plotsdir,'moving_punctures.pdf'))
            plt.show()
        except IndexError:
            print("===> Error: Time integration hasn't started yet")

    def plot_0d_output(self, var='all'):
        '''
        Input: either all variables available or a specific one, save plot or not
        Returns: plot of the variable
        '''
        norm_vars = ['alpha', 'ham', 'momx', 'momy', 'momz', 'betax', 'rpsi4','ipsi4']
        plots_0d = os.path.join(self.plotsdir,'0d_plots')
        try:
            os.mkdir(plots_0d)
        except FileExistsError:
            print('Directory exists: ',plots_0d)

        if len(os.listdir(self.out_0d_dir))==0:
            print('===> Error: No output produced yet')
        else:
            lmax = self.lmax
            if var=='all':
                paras = [i for i in os.listdir(self.outpath) if i.endswith('.par')]
                params_file = os.path.join(self.outpath, paras[0])
                params = {}
                with open(params_file) as file:
                    for i,line in enumerate(file):
                        if line.strip().startswith('#') or len(line.strip()) == 0:
                            continue
                        else:
                            try:
                                k, v = line.strip().split('=')
                            except:
                                k, v = line.strip().split('=',1) # Split on the first occurrence of '='
                                v = v.split('#', 1)[0].strip() # Remove any comments starting with '#'

                            if k.strip()=='0doutput':
                                try:
                                    v = v.split('#', 1)[0].strip()
                                except:
                                    continue
                                params[k.strip()] = v.strip()
                            else:
                                continue
                out0 = params['0doutput'].split() # variables in 0d output
                for var in out0:
                    if var in norm_vars:
                        outname = '_norm.l'
                    else:
                        outname = '_integral.l'
                    #fig = plt.figure()
                    plt.title(var)
                    for lvl in range(lmax+1):
                        try:
                            var0file = os.path.join(self.out_0d_dir, var + outname + str(lvl))
                            t0, v0 = np.loadtxt(fname=var0file, usecols=(0,1), unpack=True)
                        except OSError:
                            var0file = os.path.join(self.out_0d_dir, var + outname + str(lvl) + 'a')
                            t0, v0 = np.loadtxt(fname=var0file, usecols=(0,1), unpack=True)
                        plt.scatter(t0,v0,label=str(lvl), marker='.')
                    plt.legend()
                    plt.savefig(os.path.join(plots_0d,var+'.pdf'))
                    plt.show()

            else:
                if var in norm_vars:
                    outname = '_norm.l'
                else:
                    outname = '_integral.l'
                #fig = plt.figure()
                plt.title(var)
                for lvl in range(lmax+1):
                    var0file = os.path.join(self.out_0d_dir, var + outname + str(lvl))
                    t0, v0 = np.loadtxt(fname=var0file, usecols=(0,1), unpack=True)
                    plt.scatter(t0,v0,label=str(lvl),marker='.')
                plt.legend()
                plt.savefig(os.path.join(plots_0d,var+'.pdf'))
                plt.show() 



