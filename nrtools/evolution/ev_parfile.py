import os
import re

ID_HEADER = {
    'physics' : 'adm z4 Gauge matter eos grhd EIDdataReader Invariants AHmod ADM_mass hydroanalysis',
    # initial data
    'EIDdataReader_exe'      : '@@',
    'EIDdataReader_checkpoint'  : '@@',
    'EIDdataReader_save'     : 'yes',
    'EIDdataReader_physics'  : 'BHNS',
    'EIDdataReader_BHfiller' : 'ChebTn_Ylm_perfect_s2',
    # eos
    'eos'                    : 'pwp',
    'eos_tab_file'           : '@@',
    # checkpointing
    'checkpoint'             : 'yes',
    'checkpoint_dt_hours'    : '10',
    'checkpoint_variables'   : 'all',
    'ExitIfNAN'              : 'yes'
}

GRID_SETUP = {
    # basic setup
    'order_centered'              : '4',
    'order_advection'             : '4',
    'advection_lopsided6'         : '2',
    'advection_lopsided'          : '1',
    'order_dissipation'           : '6',
    'dissipation_factor'          : '0.5',
    'dissipation_factor_shells'   : '0.1',
    'bampi_timer_on'              : 'yes',
    'bampi_lowlatency'            : 'yes',
    'bampi_nghosts'               : '6',
    'order_RP'                    : '6',
    'order_RP_shells'             : '6',
    'amr_nbuffer'                 : '6',
    'amr_npunctures'              : '1',
    'amr_lmax'                    : '@@',
    'amr_move_nxyz'               : '@@',
    'nxyz'                        : '@@',
    'amr_move_lcube'              : '@@',
    'dxyz'                        : '@@',
    'amr_bo_lmin'                 : '4',
    'amr'                         : 'bo newfmr move',
    'amr_fmr'                     : 'nestedboxes',
    'grid'                        : 'box bitant'
}

HYDRO = {
    # hydro
    'hrsc_nghosts'                : '4',
    'hrsc_rec'                    : 'WENOZ',
    'hrsc_rec_HO'                 : 'WENOoptimal',
    'hrsc_TVD_limiter'            : 'MC2',
    'hrsc_rec_metric'             : '@@',
    'hrsc_flux'                   : '@@',
    'ExitIfNAN_vars'              : 'grhd_D grhd_Tau grhd_Sx grhd_rho grhd_epsl grhd_vx grhd_p',
    'grhd_C2P'                    : 'p_root',
    'grhd_C2P_NewtonRaphsonTR'    : '1e-9',
    'grhd_C2P_NewtonRaphsonNR'    : '100',
    'grhd_vmax'                   : '0.999',
    'grhd_Wlor_max'               : '1e5',
    'grhd_use_atmosphere'         : 'yes',
    'grhd_use_atmosphere'         : 'ColdStatic',
    'grhd_atm_factor'             : '100',
    'grhd_atm_level'              : '1e-10',
    'grhd_use_atmosphere_mask'    : 'yes',
    'grhd_use_atmosphere_prerhs'  : 'yes',
    'grhd_use_atmosphere_postrhs' : 'no',
    'matter_interpolate_order'    : '4',
    'matter_interpolate_scheme_restriction'  : 'linear',
    'matter_interpolate_scheme_prolongation' : 'linear',
    'conservative_amr'            : 'yes',
    'camr_treshold'               : '1e20',
    'grhd_use_excision'           : 'yes',
    'grhd_excision_rfct'          : '0.9',
    'grhd_excision_modus'         : 'atm'
}
    
ENTROPY_VISCOSITY = {
    # entropy-viscosity
    'grhd_visc_method_cE'          : '1',
    'grhd_visc_method_cmax'        : '1',
    'grhd_visc_method'             : 'phys',
    'grhd_visc_compute'            : 'prestep', 
    'grhd_visc_method_on_t'        : '0',
    'grhd_visc_method_off_t'       : '0'
}

EVOLUTION = {
    # evolution
    'evolution_method'            : 'rk',
    'evolution_method_rk'         : 'rk4g',
    'dtfac'                       : '0.25',
    'finaltime'                   : '20000',
    'z4_normalizedetg'            : 'yes',
    'z4_subtractA'                : 'yes',
    'z4_chi_div_floor'            : '1e-5',
    'z4_initial_lapse'            : 'donothing',
    'z4_initial_shift'            : 'zero',
    'z4_lapse'                    : '1+log withshift',
    'z4_shift'                    : 'withoutB withShiftadv',
    'z4_shiftalphapower'          : '0.0',
    'z4_shiftdriver'              : '@@', 
    'z4_shiftgammacoeff'          : '1.0',
    'z4_kappa1'                   : '0.02', 
    'z4_kappa2'                   : '0.0',
    'punctures_lapse'             : 'psiBL^(-2)'
}

OUTPUT = {
    # output
    '0douttime'                   : '@@',
    '0doutput'                    : 'alpha ham momx grhd_D grhd_rho rpsi4 ipsi4',
    '#1douttime'                   : '@@',
    '#1doutput'                    : 'alpha grhd_rho',
    '#1doutinterpolate'            : 'no',
    '#1doutputall'                 : 'no',
    '#2douttime'                   : '@@',
    '#2doutput'                    : 'alpha bssn_chi grhd_rho ham momx grhd_vx grhd_v2 grhd_epsl grhd_p hydroa_uesc hydroa_Du hydroa_Db hydroa_Dh hydroa_etot',
    '#2dformat'                    : 'vtk binary float',
    '#2doutinterpolate'            : 'no',
    '#2doutputr'                   : 'sphere_data',
    '#3douttime'                   :  '@@',
    '#3doutput'                    :  'bssn_chi grhd_rho hydroa_Du hydroa_Db',
}

BOUNDARY_AND_GAUGE = {
    # boundary
    'boundary'                               : 'background radcentered',
    # gauge
    'Gauge'                                  : 'moving_puncture',
    'moving_puncture_track_method'           : 'extint',
    'moving_puncture_fixz'                   : 'none',
    'compute_moving_puncture_distance'       : 'line',
    'compute_moving_puncture_distance_method': 'AHsimpleNSBH'
}
    
AHMOD = {
    # AHmod
    'AHmod_verbose'               : 'no',
    'AHmod_ntheta'                : '100',
    'AHmod_nphi'                  : '100',
    'AHmod_nhorizons'             : '1',
    'AHmod_searchMTS'             : '1 0.0  30000.0  1',
    'AHmod_uselast'               : 'yes',
    'AHmod_initial_guess_expand'  : '1.0',
    'AHmod_LevelOffset'           : '1',
    'AHmod_UseOptimalLevel'       : 'no',
    'AHmod_flow_iter'             : '5000',
    'AHmod_mass_tol'              : '5.0e-06',
    'AHmod_hmean_tol'             : '100.0',
    'AHmod_time'                  : '@@', 
    'AHmod_shrinking'             : '2.0',
    'AHmod_output'                : 'yes',
    'AHmod_output_xyt'            : 'no',
    'AHmod_output_lm'             : 'no'
}  
    
INVARIANTS = {
    # invariants
    'ntheta'                                  : '47',
    'nphi'                                    : '46',
    'invariants_compute_modes'                : 'yes',
    'invariants_modes_r'                      : '300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 2000 2500', 
    'invariants_modes_lmin'                   : '0',
    'invariants_modes_lmax'                   : '4',
    'invariants_energy_r'                     : '300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 2000 2500', 
    'gauss_codacci_mainardi'                  : 'standard',
    'invariants_order'                        : '4',
    'invariants_compute_modes_general'        : 'no',
    'mode_lmax'                               : '6',
    'invariants_modes_output'                 : 'all',
    'Invariants_output_time'                  : '@@'
}  

HYDROANALYSIS = {
    # hydroanalysis
    'hydroanalysis_ejecta_spheres'            : 'yes',
    'hydroanalysis_ejecta_spheres_radius'     : '200',
    'hydroanalysis_ejecta_nradius'            : '12',
    'hydroanalysis_ejecta_dradius'            : '100',
    'hydroanalysis_Mbar_radius'               : '8',
    'hydroanalysis_Mbar_nradius'              : '10',
    'hydroanalysis_Mbar_dradius'              : '0.5',
    'hydroanalysis_rATM'                      : '1e-25',
    'hydroa_mode_projection'                  : 'yes'
}   
    
ADM_MASS = {
    # ADM Mass
    'ADM_mass_ncircles'         : '101',
    'ADM_mass_npoints'          : '80',
    'ADM_mass_lmin'             : '0',
    'ADM_mass_lmax'             : '7',
    'ADM_mass_r'                : '300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 2000 2500' 
}

########################################
# Parameter File class
########################################

class Ev_Parameter_File():
    """
    ------------------
    Initialization:
    ------------------
    path          : where to put the parameter file
    ev_exe        : path/to/executable
    resolution    : resolution for the simulation
    lmax, lmax2   : refinement levels for obj1 and obj2
    flux          : flux reconstruction scheme (LLF or EFL)

    Writes the parameter file with initialization if non existent.
    """
    def __init__(self, path, ev_path, initial_data, resolution, lmax, lmax2, flux):
        self.path = path
        if len([i for i in os.listdir(self.path) if i.endswith('.par')])>0:
            EV_PARDICR = {}
            parfile = os.path.join(self.path,[i for i in os.listdir(self.path) if i.endswith('.par')][0])
            with open(parfile) as f:
                for line in f:
                    if len(line.strip())==0 or line.strip().startswith('#'):
                        continue
                    else:
                        if '#' in line: # to remove comments
                            keyval = line.strip().split('#')[0]
                            key, value = keyval.split('=')
                        else:
                            key, value = line.strip().split('=')
                        EV_PARDICR[key] = value

            self.pardic = EV_PARDICR
        else: 
            self.id_header = ID_HEADER
            self.grid_setup = GRID_SETUP
            self.hydro = HYDRO
            self.evolution = EVOLUTION
            self.output = OUTPUT
            self.boundary_and_gauge = BOUNDARY_AND_GAUGE
            self.ahmod = AHMOD
            self.invariants = INVARIANTS
            self.hydroanalysis = HYDROANALYSIS
            self.adm_mass = ADM_MASS
            self.entropy_viscosity = ENTROPY_VISCOSITY
            self.pardic = self.parameter_dictionary(ev_path, initial_data, resolution, lmax, lmax2, flux)

    def parameter_dictionary(self, ev_path, initial_data, resolution, lmax, lmax2, flux):
        # get grid specific parameters for the binary:
        grid_params = self.get_grid_params(initial_data, resolution, lmax, lmax2, flux)
        eos_tab_path = os.path.join(ev_path,'src/projects/eos/tab/pwpfits')

        ## Fill necessary values:

        # ID_HEADER
        self.id_header['EIDdataReader_exe'] = initial_data.id_exe
        self.id_header['EIDdataReader_checkpoint'] = os.path.join(os.path.join(initial_data.ou.outpath,initial_data.ou.hr_path),'checkpoint.dat')
        eos = initial_data.parfile.pardic['NS_EoS_description']#[-3:]
        self.id_header['eos_tab_file'] = os.path.join(eos_tab_path,'eos_'+eos.lower()+'.pwp')
        if os.path.exists(self.id_header['eos_tab_file'])==False:
            print('ERROR: EOS tab file not found!')

        # GRID_SETUP
        self.grid_setup['amr_lmax'] = grid_params['amr_lmax']
        if lmax>lmax2:
            self.grid_setup['amr_lmax2'] = grid_params['amr_lmax2']
        self.grid_setup['amr_move_nxyz'] = grid_params['amr_move_nxyz']
        self.grid_setup['nxyz']     = grid_params['nxyz']     
        self.grid_setup['amr_move_lcube'] = grid_params['amr_move_lcube']
        self.grid_setup['dxyz']          = grid_params['dxyz']
        #self.grid_setup['amr_bo_dxmax']  = grid_params['amr_bo_dxmax']
        
        # HYDRO
        if flux=='LLF':
            self.hydro['hrsc_rec_metric'] = 'LAG6'
            self.hydro['hrsc_flux']  = 'LLF'
        elif flux=='EFL':
            self.hydro['hrsc_rec'] = 'CENO3'
            self.hydro['hrsc_rec_metric'] = 'LAG4'
            self.hydro['hrsc_flux']  = 'HOEV2'
            self.hydro = {**self.hydro, **self.entropy_viscosity}
        
        # EVOLUTION
        self.evolution['z4_shiftdriver'] = grid_params['z4_shiftdriver']

        # OUTPUT
        self.output['0douttime'] = grid_params['douttime0']
        self.output['1douttime'] = grid_params['douttime1']
        self.output['2douttime'] = grid_params['douttime2']
        self.output['3douttime'] = grid_params['douttime2']*4
        if flux=='EFL':
            self.output['2doutput'] += ' grhd_entro grhd_viscD' 

        # OTHERS
        self.ahmod['AHmod_time'] = self.output['3douttime']
        self.invariants['Invariants_output_time'] = self.output['2douttime']
        # Radii
        radii_string = self.invariants['invariants_modes_r']
        radii = [int(num) for num in radii_string.split()]
        try:
            add_radii = int(grid_params['invariants_modes_r'])
            radii.append(add_radii)
        except:
            add_radii = [int(num) for num in grid_params['invariants_modes_r'].split()]
            radii.extend(add_radii)
        ext_radii = ' '.join(str(num) for num in sorted(radii))
        self.invariants['invariants_modes_r'] = ext_radii #grid_params['invariants_modes_r']
        self.invariants['invariants_energy_r'] = ext_radii #grid_params['invariants_energy_r']
        self.adm_mass['ADM_mass_r'] = ext_radii #grid_params['ADM_mass_r']


        # Make parfile
        self.make_parfile()
        EV_PARDIC = {**self.id_header, **self.grid_setup, **self.hydro, **self.evolution, **self.output, **self.boundary_and_gauge, **self.ahmod, **self.invariants, **self.hydroanalysis, **self.adm_mass}
        
        return EV_PARDIC


    def get_grid_params(self, initial_data, resolution, lmax, lmax2, flux, save=True):
        print('==> Setting up grid parameters for parfile ..')

        iddic = initial_data.ou.id_dic
        bhmass = float(iddic['BH_irreducible_mass_current'])
        
        if bhmass < 4:
            factor = 2
            if lmax==lmax2:
                factor = 1.5
            else:
                factor = 1.7
        elif bhmass < 5.5:
            factor = 2.5
        else:
            factor = 3

        grid_params = {}
        grid_params['amr_move_nxyz'] = resolution
        grid_params['nxyz'] = factor*grid_params['amr_move_nxyz']
        grid_params['amr_lmax'] = lmax
        grid_params['amr_lmax2'] = lmax2
        grid_params['amr_move_lcube'] = 2
        grid_params['dtfac'] = 0.25

        assert(grid_params['amr_lmax'] >= grid_params['amr_lmax2'])

        # read these params:
        params=["BH_Christodoulou_mass_current","BH_max_radius","NS_ADM_mass","NS_max_radius","BHNS_separation"]

        # read all file and get the parameters
        id_file = initial_data.ou.txt_path 
        param_dict = {}
        with open(id_file,"r") as fp:
            line = fp.readline()
            while line:
                line = fp.readline()
                for p in params:
                    if re.search(r"^{}[\W]".format(p),line):
                        # trim new line and white spaces
                        line = line.replace('\n', "")
                        line = line.replace(" ", "")
                        param_dict[p] = (float(line.split('=')[1]))
        
        # diameter:
        BH_d = 2*param_dict["BH_max_radius"]
        NS_d = 2*param_dict["NS_max_radius"]
        NS_m = param_dict["NS_ADM_mass"]
        BH_m = param_dict["BH_Christodoulou_mass_current"]
        init_s= param_dict["BHNS_separation"]

        ## compute:
        ## note: to output correctly you should compute things from finest level
        dxyz_fine = round((NS_d*(1.15))/(grid_params['amr_move_nxyz']),6) # NOTE: 15% bigger of NS diameter
        grid_params['dxyz'] = dxyz_fine * (2**grid_params['amr_lmax2']) # coarsest grid space
        grid_params['amr_bo_dxmax'] = (2.4/grid_params['amr_move_nxyz'])*64
        grid_params['douttime0']    = grid_params['dxyz']/(2**(grid_params['amr_lmax'])) * 2**int(grid_params['amr_lmax']/3)/grid_params['dtfac']
        grid_params['douttime1']    = grid_params['dxyz']/(2**(grid_params['amr_lmax'])) * 2**int(grid_params['amr_lmax']/3)/grid_params['dtfac']
        grid_params['douttime2']    = grid_params['dxyz']/(2**(grid_params['amr_lmax'])) * 2**int(grid_params['amr_lmax']/3)/grid_params['dtfac']*4
        grid_params['z4_shiftdriver'] = 2.0/(BH_m+NS_m)

        NS_box_len = grid_params['dxyz']/(2**(grid_params['amr_lmax2'])) * grid_params['amr_move_nxyz']
        BH_box_len = grid_params['dxyz']/(2**(grid_params['amr_lmax'])) * grid_params['amr_move_nxyz']

        BH_box_len2= "not defined, since there is only one moving box"
        if (grid_params['amr_lmax'] - grid_params['amr_move_lcube'] > 1):
            BH_box_len2  = grid_params['dxyz']/(2**(grid_params['amr_lmax']-1)) * grid_params['amr_move_nxyz']
        else:
            BH_box_len2  = grid_params['dxyz']/(2**(grid_params['amr_lmax']-1)) * grid_params['nxyz']

        # radii for GW, ADM etc.
        radii = ""
        outer_r_int = int(grid_params['dxyz']*grid_params['nxyz']/2)
        if outer_r_int > 1500:
            for r in range(300,1500,200):
                radii += "{} ".format(r)
        if outer_r_int > 6000:
            for r in range(1500,6000,500):
                radii += "{} ".format(r)
        for r in range(6000,outer_r_int,1000):
            radii += "{} ".format(r)
        radii += "{} ".format(outer_r_int)

        grid_params['invariants_modes_r'] = radii
        grid_params['invariants_energy_r'] = radii
        grid_params['ADM_mass_r'] = radii
        if save:
            file = open(os.path.join(self.path,'grid_setup.log'), 'a')
            file.write("==> Grid setup: \n")
            file.write("dxyz           = {} # NS_diameter*(1.15))/(amr_move_nxyz)*(2**amr_lmax2) \n".format((grid_params['dxyz'])))
            file.write("amr_bo_dxmax   = {} # (2.4/amr_move_nxyz)*64  \n".format((grid_params['amr_bo_dxmax'])))
            file.write("z4_shiftdriver = {} # 2.0/(BH_m+NS_m) \n".format((grid_params['z4_shiftdriver'])))
            file.write("0douttime    = {} # dxyz/(2**(amr_lmax)) * 2**int(amr_lmax/3)/dtfac \n".format((grid_params['douttime0'])))
            file.write("1douttime    = {} # dxyz/(2**(amr_lmax)) * 2**int(amr_lmax/3)/dtfac \n".format((grid_params['douttime1'])))
            file.write("2douttime    = {} # dxyz/(2**(amr_lmax)) * 2**int(amr_lmax/3)/dtfac*4 \n".format((grid_params['douttime2'])))
            file.write("AHmod_time   = {} # 2douttime \n".format((grid_params['douttime2'])))
            file.write("Invariants_output_time = {} # 2douttime \n".format((grid_params['douttime2'])))
            file.write("extraction radii = {} ".format(radii))

            file.write("\n\n## Configs: \n")

            ## calculate moving boxes consistency:
            ## the size of the smallest non moving level (box)
            small_non_move_box = dxyz_fine* 2**(grid_params['amr_lmax2'] - grid_params['amr_move_lcube'])* grid_params['nxyz']
            ## the biggest moving box
            big_move_box = dxyz_fine * 2**(grid_params['amr_lmax2'] - grid_params['amr_move_lcube']-1)*grid_params['amr_move_nxyz']

            ## 2*big_move_box, coeffs 2 in order to account for buffers
            if (small_non_move_box - (2*big_move_box + init_s) > 0):
                file.write("# MOVING BOXES ARE OK. \n")
            else:
                file.write("# WARNING: MOVING BOXES DO NOT FIT INTO FIXED LEVEL! \n")

            file.write("# BHNS separation = {} \n".format(init_s))
            file.write("# dxyz_finest_NS  = {} # dxyz in finest level around the NS after roundoff \n".format(dxyz_fine))
            file.write("# dxyz_finest_BH  = {} # dxyz in finest level around the BH after roundoff \n".format(dxyz_fine/2**(grid_params['amr_lmax']-grid_params['amr_lmax2'])))
            file.write("# BH_m = {} \n".format(BH_m))
            file.write("# NS_m = {} \n".format(NS_m))
            file.write("# BH_diameter = {} \n".format(BH_d))
            file.write("# NS_diameter = {} \n".format(NS_d))
            file.write("# NS_finest_box_len = {} \n".format(NS_box_len))
            file.write("# BH_finest_box_len = {} \n".format(BH_box_len))
            if (grid_params['amr_lmax'] - grid_params['amr_move_lcube'] > 1):
                file.write("# BH_next_finest_box_len = {} > {} ???"
                        " should BH_next_finest_box_len > (1.25)*BH_diameter \n".
                        format(BH_box_len2,(1.25)*BH_d))

                file.write("# BH_m/h_min|finest     = {} > 20 ??? "
                    " see https://arxiv.org/pdf/1007.4789.pdf  \n".
                    format(BH_m/(grid_params['dxyz']/2**grid_params['amr_lmax'])))

                file.write("# BH_m/h_min|2nd_finest = {} > 20 ??? "
                    " see https://arxiv.org/pdf/1007.4789.pdf  \n".
                    format(BH_m/(grid_params['dxyz']/2**(grid_params['amr_lmax']-1))))
                    
                file.write("# outer boundary radius ~ {} > {} ???"
                    " should dxyz*nxyz/2 > 5*init_separation \n".
                    format(grid_params['dxyz'] * grid_params['nxyz']/2,5*init_s))
                
            file.close()
    
        return grid_params
                
    def make_parfile(self):
        # Write par file
        with open(os.path.join(self.path,'bam_evo.par'), 'w') as f:
            f.write("\n############################################################################# \n# ID_HEADER \n")
            for key, value in self.id_header.items():
                f.write('%s =   %s\n' % (key, value))
            f.write("\n############################################################################# \n# GRID_SETUP \n")
            for key, value in self.grid_setup.items():
                f.write('%s =   %s\n' % (key, value))    
            f.write("\n############################################################################# \n# HYDRO \n")
            for key, value in self.hydro.items():
                f.write('%s =   %s\n' % (key, value))
            f.write("\n############################################################################# \n# EVOLUTION \n")
            for key, value in self.evolution.items():
                f.write('%s =   %s\n' % (key, value))
            f.write("\n############################################################################# \n# OUTPUT \n")
            for key, value in self.output.items():
                f.write('%s =   %s\n' % (key, value))
            f.write("\n############################################################################# \n# BOUNDARY_AND_GAUGE \n")
            for key, value in self.boundary_and_gauge.items():
                f.write('%s =   %s\n' % (key, value))
            f.write("\n############################################################################# \n# AHMOD \n")
            for key, value in self.ahmod.items():
                f.write('%s =   %s\n' % (key, value))
            f.write("\n############################################################################# \n# INVARIANTS \n")
            for key, value in self.invariants.items():
                f.write('%s =   %s\n' % (key, value))
            f.write("\n############################################################################# \n# HYDROANALYSIS \n")
            for key, value in self.hydroanalysis.items():
                f.write('%s =   %s\n' % (key, value))
            f.write("\n############################################################################# \n# ADM_MASS \n")
            for key, value in self.adm_mass.items():
                f.write('%s =   %s\n' % (key, value))
