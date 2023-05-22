import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from itertools import chain
import pyvista as pv
import imageio
from watpy.wave.wave import mwaves

########################################
# 2D Data class
########################################
class Data2D():
    '''
    Single field 2d data on multiple ref.levels and boxes
    '''    
    def __init__(self, path = '.', var = 'grhd_rho.xy', 
                 lev = [5,6] # each level can be made of one of more box (a,b,...)
                 ,outdir = 'viz', verbose = True):
        self.path = os.path.abspath(path)
        self.var = var
        self.lev = lev
        self.plane = var.split('.')[-1]
        self.outdir = os.path.abspath(outdir)
        self.verbose = verbose

        os.makedirs("{}/{}".format(self.outdir,self.var),
                    exist_ok = True)

        # collect data at each iteration
        self.timelev = dict.fromkeys(self.lev)
        for l in self.lev:
            self.timelev[l] = {}
            if self.verbose: print('=> level {}'.format(l))
            files = glob.glob("{}/{}{}_vtk/{}{}*_????.vtk".format(self.path,self.var,l,var,l))
            print("{}/{}{}_vtk/{}{}*_????.vtk".format(self.path,self.var,l,var,l))
            #print(self.path+'/'+self.var)
            for f in files:
                grid,iteration = os.path.basename(f).split('.')[-2].split('_') # .xy6a_0000. -> 'xy6a','0000'
                iteration = int(iteration)
                if self.verbose: print(' iteration {} grid {} file {}'.format(iteration,grid,f))
                if iteration not in self.timelev[l]:
                    self.timelev[l][iteration] = {}
                    self.timelev[l][iteration]['files'] = []
                self.timelev[l][iteration]['files'].append(f)
            
    def get_iterations(self):
        iterations = np.array(list(chain.from_iterable(self.timelev[l].keys() for l in self.timelev.keys())))
        return np.unique(iterations)

    def plot_slice(self, iteration,
                   add_warp=1, # Modifies point coordinates by moving points along point normals by the scalar amount times the scale factor.
                   show_edges=True, # Show the edges of all geometries within a mesh
                   show_grid=True, # Show gridlines and axes labels.
                   add_axes=True, # Add an interactive axes widget in the bottom left corner.
                   show=False):
        
        #pv.global_theme.show_edges = show_edges
        #pv.global_theme.edge_color = 'white'
        scalar_bar_args={'title': self.var}
        
        pl = pv.Plotter(off_screen=not show)
        data = []
        for l in self.lev:
            if iteration in self.timelev[l].keys():
                for df in self.timelev[l][iteration]['files']:
                    data = pv.read(df)
                    if add_warp > 1: 
                        data = data.warp_by_scalar(factor=add_warp)
                    pl.add_mesh(data, scalar_bar_args=scalar_bar_args, label='RL{}'.format(l))#, log_scale=True)
                    if show_edges:
                        edges = data.extract_all_edges()
                        act = pl.add_mesh(edges, line_width=0.001, color='w',opacity=0.5)
                
        _ = pl.add_title('iter = {}'.format(iteration),font_size=12)
        if show_grid: _ = pl.show_grid()
        if add_axes: _ = pl.add_axes(line_width=5)
        if show: pl.show(screenshot='{}/{}/{:06d}.png'.format(self.outdir,self.var,iteration))
        else: pl.screenshot('{}/{}/{:06d}.png'.format(self.outdir,self.var,iteration))
        pl.close()
        
    def plot_slices(self,imin=0,imax=np.inf,di=1,
                    add_warp=1,
                    show_edges=True,
                    show_grid=True,
                    add_axes=True,
                    show=False):
        
        if self.verbose: print('Rendering...')
        for it in self.get_iterations():
            if it < imin: continue
            if it > imax: continue
            if it % di: continue
            if self.verbose: print('iter = {}'.format(it))
            self.plot_slice(it,
                            add_warp=add_warp,
                            show_edges=show_edges,
                            show_grid=show_grid,
                            add_axes=add_axes,
                            show=show)

    def make_movie(self,codec="mpeg4"):
        cwd = os.getcwd()
        os.chdir("{}/{}".format(self.outdir,self.var))
        if codec == "mpeg4":
            os.system("ffmpeg -framerate 25 -i %06d.png -codec mpeg4 movie.mp4")
        os.chdir(cwd)

    def make_gif(self):
        image_dir = "{}/{}".format(self.outdir,self.var)
        images = [os.path.join(image_dir, f) for f in os.listdir(image_dir) if f.endswith('.png')]
        images = sorted(images)
        output_file = image_dir+'/'+self.var+'.gif'
        with imageio.get_writer(output_file, mode='I') as writer:
            for image in images:
                writer.append_data(imageio.imread(image))


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

    def get_horizon_area(self):
        '''
        Returns time and the coordinate area of the apparent horizon
        '''
        try:
            hfile = os.path.join(self.outpath,'horizon_0')
            t, ca = np.loadtxt(fname=hfile, comments='#', usecols=(0,9), unpack=True)
        except IndexError:
            t, ca = None
            print("===> Error: Time integration hasn't started yet")
        return t, ca

    def plot_horizon_area(self):
        try:
            hfile = os.path.join(self.outpath,'horizon_0')
            t, ca = np.loadtxt(fname=hfile, comments='#', usecols=(0,9), unpack=True)
            plt.scatter(t,ca, marker='.',label='AH coord. area')
            plt.grid()
            plt.legend()
            plt.savefig(os.path.join(self.plotsdir,'AH_coord_area.pdf'))
            plt.show()
        except:
            print("===> Error: Time integration hasn't started yet")

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

    def get_0d_variable(self, var):
        '''
        Returns variable at the finest level where the puncture is
        '''
        norm_vars = ['alpha', 'ham', 'momx', 'momy', 'momz', 'rpsi4','ipsi4']
        if len(os.listdir(self.out_0d_dir))==0:
            print('===> Error: No output produced yet')
            t0 = v0 = None
        else:
            lmax = self.lmax
            if var in norm_vars:
                outname = '_norm.l'
            else:
                outname = '_integral.l'
            var0file = os.path.join(self.out_0d_dir, var + outname + str(lmax) + 'b')
            t0, v0 = np.loadtxt(fname=var0file, usecols=(0,1), unpack=True)
        return t0, v0

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

    def do_2d_movies(self,var='all'):
        movies_out = os.path.join(self.outpath,'viz')
        try:
            os.mkdir(movies_out)
        except FileExistsError:
            print('Directory exists: ',movies_out)

        if var=='all':
            vars = ['grhd_rho.xy','alpha.xy','bssn_chi.xy','grhd_epsl.xy','grhd_p.xy','grhd_v2.xy','grhd_vx.xy','grhd_D.xy','hydroa_Db.xy','hydroa_Dh.xy','hydroa_Du.xy','hydroa_etot.xy','hydroa_uesc.xy']
            for var in vars:
                data2d = Data2D(path = self.out_2d_dir, var=var, outdir=movies_out)
                data2d.plot_slices() 
                data2d.make_movie()
                data2d.make_gif()
        else:
            data2d = Data2D(path = self.out_2d_dir, var=var, outdir=movies_out)
            data2d.plot_slices() 
            data2d.make_movie()
            data2d.make_gif()

    def get_mp_Rpsi4(self,Mtot,Momg22):
        '''
        Get watpy multipolar wave object
        Mtot: binary gravitational mass (initial_data output)
        Momg22: initial GW frequency in geometric units (id_gw_frequency_Momega22)
        '''
        try:
            fnames = [os.path.split(x)[1] for x in glob.glob('{}/{}'.format(self.out_inv_dir,'Rpsi4mode??_r*.l0'))]
            f0 = Momg22 / 2*np.pi / Mtot
            mwave = mwaves(path = self.out_inv_dir, code = 'bam', filenames = fnames, mass = Mtot, f0 = f0, ignore_negative_m=True)
        except IndexError:
            mwave = None
            print("===> Error: Time integration hasn't started yet")
        return mwave


