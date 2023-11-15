import os

PARDIC = {
    'Project' : 'BH_NS_binary_initial_data',
    'BHNS_separation'                  : '@@', # fill data
    'BHNS_angular_velocity'           : 'auto',
    'BHNS_infall_velocity'             : '0',
    'BHNS_start_off'                   : 'parameter_file',
    'BHNS_observe_ADM_P'               : 'S_obj1+S_obj2,default',
    'BHNS_observe_ADM_J'               : 'S_obj1+S_obj2,default',
    'BHNS_P_ADM_control_method'        : 'adjust(x_CM,y_CM)',
    'BHNS_P_ADM_control_update_weight' : '0.(x10)->0.2',
    'BHNS_P_ADM_control_tolerance'     : '1E-8',
    'BHNS_P_ADM_control_threshold'     : '10',
    'BH_irreducible_mass'        : '@@', # fill data
    'BH_chi_x'                   : '@@', # fill data
    'BH_chi_y'                   : '@@', # fill data
    'BH_chi_z'                   : '@@', # fill data
    'BH_boost_Vx'                : 'off',
    'BH_Eq_inner_BC_fields'      : 'XCTS',
    'BH_Eq_inner_BC_beta'        : 'alpha+Omega*r',
    'BH_Eq_inner_BC_alpha'       : 'none',
    'BH_start_off'               : 'IsoSchild',
    'BH_surface_type'            : 'perfect_s2',
    'BH_surface_Ylm_max_l'       : '1',
    'BH_tune_BH_radius_criteria' : 'fix_irreducible_mass',
    'BH_mass_tolerance'          : '1E-8',
    'BH_radius_update_weight'    : '0.(x10)->0.2(x390)->0.',
    'BH_spin_update_weight'      : '0.(x10)->0.2(x390)->0.',
    'BH_spin_tolerance'          : '1E-3',
    'NS_baryonic_mass'               : '@@', # fill data
    'NS_EoS_description'            : '@@', # 'K', 'SLy'
    'NS_EoS_type'                    : '@@', # 'polytropic', 'piecewise_polytropic'
    'NS_EoS_unit'                    : 'geo', # 'geo', 'geo'
    'NS_EoS_K0'                      : '@@', # '[92.12]', '[8.9493e-02]'
    'NS_EoS_Gamma'                   : '@@', # '[2]', '[1.3569e+00, 3.0050e+00, 2.9880e+00, 2.8510e+00]'
    'NS_EoS_rho0_th'                 : '@@', # [0, 2.3674e-04, 8.1147e-04, 1.6191e-03] only for PWP EOS
    'NS_Omega_x'             	       : '0.', 
    'NS_Omega_y'                     : '0.', 
    'NS_Omega_z'                     : '0.', 
    'NS_surface_type'                : 'perfect_s2->topology_s2',
    'NS_surface_finder'              : 'bisection',
    'NS_surface_change_threshold'    : '0.0',
    'NS_surface_Ylm_max_l'           : '10',
    'NS_enthalpy_allowed_residual'   : '1E-8',
    'NS_enthalpy_update_weight'      : '0.1',
    'NS_Euler_const_criteria'        : 'fix_baryonic_mass',
    'NS_Euler_const_update_weight'   : '1.',
    'NS_force_balance_equation'      : 'none(x4)->adjust(d/dy:Omega)',
    'NS_force_balance_update_weight' : '0.2',
    'NS_adjust_center_method'        : 'taylor_expansion',
    'NS_adjust_center_update_weight' : '1.',
    'NS_extrapolate_matter_fields'   : 'inverse_r_expmr',
    'NS_Eq_phi_polish'               : '0.1',
    'NS_start_off'                   : 'TOV',
    'SYS_initialize'        : 'TOV+IsoSchild',
    'SYS_initialize_fields' : 'XCTS',
    'Free_data_conformal_metric'             : 'flat',
    'Free_data_conformal_Christoffel_symbol' : 'flat',
    'Free_data_conformal_Ricci'              : 'flat',
    'Free_data_trK'                          : 'zero',
    'ADM_constraints_method'     : 'from_scratch',
    'ADM_B1I_form'               : 'inspiral',
    'ADM_compute_adm_Kuu_method' : 'use_AIJ',
    'Tij_NS_decomposition' : 'XCTS',
    'Tij_NS_gConf'         : 'general',
    'checkpoint_every' : '0h',
    'Derivative_Method'               : 'Spectral',
    'Interpolation_Method'            : 'Spectral',
    'Fourier_Transformation_Method'   : 'RFT',
    'dF/du_for_Newton_Method'         : 'Spectral',
    'grid_kind'                    : 'SplitCubedSpherical(BH+NS)',
    'grid_set_NS'                  : 'left',
    'grid_set_BH'                  : 'right,excised',
    'grid_NS_central_box_length'   : 'auto',
    'grid_BH_central_box_length'   : 'auto',
    'grid_outermost_radius'        : '1E5',
    'grid_verbose'                 : 'no',
    'n_a'                          : '8(x200) ->10(x100)->12(x100)->14(x100)->16(x100)->18(x100)->20(x100)->22(x100)',
    'n_b'                          : '8(x200) ->10(x100)->12(x100)->14(x100)->16(x100)->18(x100)->20(x100)->22(x100)',
    'n_c'                         : '12(x200)->14(x100)->16(x100)->18(x100)->20(x100)->22(x100)->24(x100)->26(x100)',
    'grid_SplitCS_max_n_a'         : '40',
    'grid_SplitCS_max_n_b'         : '40',
    'grid_SplitCS_max_n_c'         : '40',
    'Eq_type'                         : 'Elliptic',
    'Eq_elliptic_test'                : 'no',
    'Eq_phi'                          : 'XCTS_curve_Type3_DDM, NS',
    'Eq_psi'                          : 'XCTS_curve_excision_Type1_DDM, .*',
    'Eq_alphaPsi'                     : 'XCTS_curve_excision_Type2_DDM, .*',
    'Eq_B0_U0'                        : 'XCTS_flat_excision_Type1_DDM , .*',
    'Eq_B0_U1'                        : 'XCTS_flat_excision_Type1_DDM , .*',
    'Eq_B0_U2'                        : 'XCTS_flat_excision_Type1_DDM , .*',
    'Eq_update_method'                : 'relaxed_scheme',
    'Eq_update_weight_phi'            : '0.2',
    'Eq_update_weight_psi'            : '0.2',
    'Eq_update_weight_alphaPsi'       : '0.2',
    'Eq_update_weight_B0_U0'          : '0.2',
    'Eq_update_weight_B0_U1'          : '0.2',
    'Eq_update_weight_B0_U2'          : '0.2',
    'solve_Order'                     : 'psi,alphaPsi,B0_U0,B0_U1,B0_U2,phi',
    'solve_Newton_Update_Weight'      : '1.',
    'solve_residual'                  : '1E-10',
    'solve_residual_factor'           : '1E-5',
    'solve_Max_Iteration'             : '1',
    'solve_Max_Newton_Step'           : '1',
    'solve_Method'                    : 'DDM_Schur_Complement',
    'solve_UMFPACK_refinement_step'   : '0',
    'solve_UMFPACK_size'              : '1',
    'txt_output_0d'  : 'ham,mom,eq_residual',
    'txt_output_1d'  : '^phi,^psi,^alphaPsi,^B0,^beta,eq_residual,ham,mom',
    'txt_output_1d_line' : '(X,0.5,0.5),(0.5,Y,0.5),(0.5,0.5,Z)'
}

########################################
# Parameter File class
########################################

class Parameter_File():
    """
    ------------------
    Initialization:
    ------------------
    path    : where the initial data should be produced
    params  : dictionary with the basic parameters of the binary
    """
    def __init__(self, path, params):
        self.path = path
        if params==None:
            PARDICR = {}
            parfile = os.path.join(self.path,[i for i in os.listdir(self.path) if i.endswith('.par')][0])
            with open(parfile) as f:
                for line in f:
                    if line.startswith('#') or len(line)==0:
                        continue
                    else:
                        try:
                            key, value = line.strip().split('=')
                            if '#' in value:
                                value = value.split('#')[0].strip()
                            try:
                                float(value)
                            except ValueError:
                                value = value.split('->')[-1]
                            key = key.strip()
                            PARDICR[key] = value
                        except:
                            print('skip..')
            self.pardic = PARDICR
        else:
            self.user_params = params
            self.pardic = self.parameter_dictionary()

    def parameter_dictionary(self):
        params = self.user_params
        PARDIC['BHNS_separation'] = params['binary_separation']
        # black hole
        PARDIC['BH_chi_x'] = 0.
        PARDIC['BH_chi_y'] = 0.
        PARDIC['BH_chi_z'] = params['bh_chi_z']
        PARDIC['BH_irreducible_mass'] = params['bh_mass']
        # neutron star
        PARDIC['NS_baryonic_mass'] = params['ns_mass']
        if params['eos'] == 'SLy':
            PARDIC['NS_EoS_description'] = 'SLy'
            PARDIC['NS_EoS_type'] = 'piecewise_polytropic'
            PARDIC['NS_EoS_K0'] = '[8.95133496e-02]'
            PARDIC['NS_EoS_Gamma'] = '[1.35692395e+00,3.00500000e+00,2.98800000e+00,2.85100000e+00]'
            PARDIC['NS_EoS_rho0_th'] = '[0.00000000e+00,2.36719491e-04,8.11322219e-04,1.61880065e-03]'
        elif params['eos'] == 'MS1b':
            PARDIC['NS_EoS_description'] = 'MS1b'
            PARDIC['NS_EoS_type'] = 'piecewise_polytropic'
            PARDIC['NS_EoS_K0'] = '[8.94989354e-02]'
            PARDIC['NS_EoS_Gamma'] = '[1.3569e+00, 3.4560e+00, 3.0110e+00, 1.4250e+00]'
            PARDIC['NS_EoS_rho0_th'] = '[0.00000000e+00, 1.8401e-04, 8.1147e-04, 1.6191e-03]'
        else:
            print('===> EoS not valid !')
        
        return PARDIC
