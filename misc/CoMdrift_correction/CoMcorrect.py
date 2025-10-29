from nrtools.initialdata.initialdata import *
from nrtools.evolution.evolution import *
from watpy.wave.gwutils import *
from watpy.wave.wave import wave
import os
import argparse
import h5py

def com_motion(time,mbh,mb,px_ns, py_ns, pz_ns, px_bh, py_bh, pz_bh):
        t = time#[time,time]
        x_A = [px_bh,py_bh,pz_bh]
        x_B = [px_ns,py_ns,pz_ns]
        m_A = mbh # Christodolou mass
        m_B = mb # rest mass NS
        m = m_A + m_B
        CoM = ((m_A * x_A) + (m_B * x_B)) / m

        # The waveform will be in units of the total mass, so we need to convert for compatibility
        t /= m
        CoM /= m

        return t, CoM


def estimate_avg_com_motion(time,mbh,mb,px_ns, py_ns, pz_ns, px_bh, py_bh, pz_bh,
    skip_beginning_fraction=0.01,
    skip_ending_fraction=0.10,
    plot=False,
    fit_acceleration=False,
):
    """Calculate optimal translation and boost from Horizons.h5

    This returns the optimal initial position and velocity such that the CoM is
    best approximated as having these initial values.  If the coordinate system
    is transformed by these values, the new CoM motion will be as close to the
    origin as possible (in the sense of squared distance from the origin
    integrated over time).

    The translation to be applied to the data should be calculated given the
    values returned by this function as

        np.array([x_i+v_i*t_j+0.5*a_i*t_j**2 for t_j in t])


    Parameters
    ----------
    path_to_horizons_h5: string, optional [default is 'Horizons.h5']
        Absolute or relative path to 'Horizons.h5'
    skip_beginning_fraction : float, optional
        Exclude this portion from the beginning of the data.  Note that this is
        a fraction, rather than a percentage.  The default value is 0.01,
        meaning the first 1% of the data will be ignored.
    skip_ending_fraction : float, optional
        Exclude this portion from the end of the data.  Note that this is a
        fraction, rather than a percentage.  The default value is 0.10, meaning
        the last 10% of the data will be ignored.
    plot : bool, optional
        If True, save plot showing CoM tracks before and after offset and
        translation, to file `CoM_before_and_after_translation.pdf` in the same
        directory as Horizons.h5.  Default: False.
    fit_acceleration: bool, optional
        If True, allow for an acceleration in the fit, and return as third
        parameter.  Default: False.
    path_to_matter_h5: None or string [default is None]
        Absolute or relative path to 'Matter.h5' for neutron star trajectories.
        If None, we assume that only black holes are involved, so that both A
        and B trajectories will be in Horizons.h5.
    m_A: None or float [default is None]
        If None, the mass will be read from Horizons.h5 (for BHs) or Matter.h5
        (for NSs).
    m_B: None or float [default is None]
        If None, the mass will be read from Horizons.h5 (for BHs) or Matter.h5
        (for NSs).

    Returns
    -------
    x_i : length-3 array of floats
        Best-fit initial position of the center of mass
    v_i : length-3 array of floats
        Best-fit initial velocity of the center of mass
    a_i : length-3 array of floats
        Best-fit initial acceleration of the center of mass [only if
        `fit_acceleration=True` is in input arguments]
    t_i : float
        Initial time used.  This is determined by the `skip_beginning_fraction`
        input parameter.
    t_f : float
        Final time used.  This is determined by the `skip_ending_fraction`
        input parameter.

    """
    import os.path
    import numpy as np
    from scipy.integrate import simps

    t, com = com_motion(time,mbh,mb,px_ns, py_ns, pz_ns, px_bh, py_bh, pz_bh)
    print(f"com shape: {com.shape}")
    print(f"t shape: {t.shape}")

    # We will be skipping the beginning and end of the data;
    # this gives us the initial and final indices
    t_i, t_f = t[0] + (t[-1] - t[0]) * skip_beginning_fraction, t[-1] - (t[-1] - t[0]) * skip_ending_fraction
    i_i, i_f = np.argmin(np.abs(t - t_i)), np.argmin(np.abs(t - t_f))

    # Find the optimum analytically
    tt = np.array([t,t])
    print(f"t shape: {tt.shape}")

    CoM_0 = simps(com[:,i_i : i_f + 1], t[i_i : i_f + 1], axis=1)
    CoM_1 = simps((t[None, :] * com)[:,i_i : i_f + 1], t[i_i : i_f + 1], axis=1)
    if fit_acceleration:
        CoM_2 = simps((t[:, np.newaxis] ** 2 * com)[i_i : i_f + 1], t[i_i : i_f + 1], axis=0)
        x_i = (
            3
            * (
                CoM_0
                * (3 * t_f ** 4 + 12 * t_f ** 3 * t_i + 30 * t_f ** 2 * t_i ** 2 + 12 * t_f * t_i ** 3 + 3 * t_i ** 4)
                - 12 * CoM_1 * (t_f + t_i) * (t_f ** 2 + 3 * t_f * t_i + t_i ** 2)
                + CoM_2 * (10 * t_f ** 2 + 40 * t_f * t_i + 10 * t_i ** 2)
            )
            / (t_f - t_i) ** 5
        )
        v_i = (
            12
            * (
                -3 * CoM_0 * (t_f + t_i) * (t_f ** 2 + 3 * t_f * t_i + t_i ** 2)
                + CoM_1 * (16 * t_f ** 2 + 28 * t_f * t_i + 16 * t_i ** 2)
                + CoM_2 * (-15 * t_f - 15 * t_i)
            )
            / (t_f - t_i) ** 5
        )
        a_i = (
            60
            * (CoM_0 * (t_f ** 2 + 4 * t_f * t_i + t_i ** 2) + CoM_1 * (-6 * t_f - 6 * t_i) + 6 * CoM_2)
            / (t_f - t_i) ** 5
        )
    else:
        x_i = 2 * (CoM_0 * (2 * t_f ** 3 - 2 * t_i ** 3) + CoM_1 * (-3 * t_f ** 2 + 3 * t_i ** 2)) / (t_f - t_i) ** 4
        v_i = 6 * (CoM_0 * (-t_f - t_i) + 2 * CoM_1) / (t_f - t_i) ** 3
        a_i = 0.0

    # If desired, save the plots
    if plot:
        import matplotlib as mpl

        #mpl.use("Agg")  # Must come after importing mpl, but before importing plt
        import matplotlib.pyplot as plt

        delta_x = np.array([x_i + v_i * t_j + 0.5 * a_i * t_j ** 2 for t_j in t])
        comprm = com.T - delta_x
        max_displacement = np.linalg.norm(delta_x, axis=1).max()
        max_d_color = min(1.0, 10 * max_displacement)
        fig = plt.figure(figsize=(10, 7))
        plt.plot([t,t], com, alpha=0.25, lw=1.5)
        plt.gca().set_prop_cycle(None)
        lineObjects = plt.plot(t, comprm, lw=2)
        plt.xlabel(r"Coordinate time")
        plt.ylabel(r"CoM coordinate values")
        plt.legend(iter(lineObjects), ("x", "y", "z"), loc="upper left")
        #plt.show()
        plt.savefig("CoM_before_and_after_transformation.pdf")

        # plot tracks
        plt.plot(com)
        plt.savefig("puncture_tracks_CoMcorrect.pdf")

    #print("Optimal x_i: [{}, {}, {}]".format(*x_i.T))
    print("Optimal x_i: {}".format(x_i))
    #print("Optimal v_i: [{}, {}, {}]".format(*v_i.T))
    print("Optimal v_i: {}".format(v_i))
    if fit_acceleration:
        print("Optimal a_i: [{}, {}, {}]".format(*a_i))
    print(f"t_i, t_f: {t_i}, {t_f}")

    if fit_acceleration:
        return x_i, v_i, a_i, t_i, t_f
    else:
        return x_i, v_i, t_i, t_f
    
def remove_avg_com_motion(
    name="simulation",
    path_to_waveform_h5="rhOverM_Asymptotic_GeometricUnits.h5/Extrapolated_N2.dir",
    skip_beginning_fraction=0.01, skip_ending_fraction=0.10,plot=False,file_write_mode="w",
    ev_output=None, m_A=None,m_B=None,file_format="NRAR"):
    """Rewrite waveform data in center-of-mass frame

    This simply uses `estimate_avg_com_motion`, and then transforms to that
    frame as appropriate.  Most of the options are simply passed to that
    function.  Note, however, that the path to the Horizons.h5 file defaults to
    the directory of the waveform H5 file.

    Additional parameters
    ---------------------
    name : str, name of simulation
    path_to_waveform_h5 : str, optional
        Absolute or relative path to SpEC waveform file, including the
        directory within the H5 file, if appropriate.  Default value is
        'rhOverM_Asymptotic_GeometricUnits.h5/Extrapolated_N2.dir'.
    path_to_horizons_h5: None or string [default is None]
        Absolute or relative path to 'Horizons.h5'.  If None, this will be
        inferred from the waveform path.
    path_to_matter_h5: None or string [default is None]
        Absolute or relative path to 'Matter.h5' for neutron star trajectories.
        If None, this will be inferred from the waveformpath.  If the file does
        not exist, we assume that only black holes are involved, so that both A
        and B trajectories will be in Horizons.h5.
    m_A: None or float [default is None]
        If None and there is no Matter.h5 file, the mass will be read from
        Horizons.h5; otherwise, the mass will be read from the metadata
        `reference_mass1`.
    m_B: None or float [default is None]
        If None and there is no Matter.h5 file, the mass will be read from
        Horizons.h5; otherwise, the mass will be read from the metadata
        `reference_mass2`.
    file_format: 'NRAR' or 'RPXM' [default is 'NRAR']
        The file format of the waveform data H5 file. The 'NRAR' format is the
        default file format found in the SXS Catalog. The 'RPXM' format is data
        compressed with the rotating_paired_xor_multishuffle_bzip2 scheme.

    Returns
    -------
    w_m : WaveformModes object
        This is the transformed object in the new frame

    """
    import os.path
    import pathlib
    import re
    import numpy as np
    from scri.SpEC.file_io import read_from_h5, write_to_h5, rotating_paired_xor_multishuffle_bzip2
    from scri import h, hdot, psi4, psi3, psi2, psi1, psi0

    directory = os.path.dirname(os.path.abspath(path_to_waveform_h5.split(".h5", 1)[0] + ".h5"))
    subdir = os.path.basename(path_to_waveform_h5.split(".h5", 1)[1])

    # Read the waveform data in
    w_m = read_from_h5(path_to_waveform_h5)

    # Get the CoM motion from Horizons.h5
    px_ns, py_ns, pz_ns, px_bh, py_bh, pz_bh, time = ev_output.extract_objects_tracks(xyz=True)
    x_0, v_0, t_0, t_f = estimate_avg_com_motion(time,mbh,mb,px_ns, py_ns, pz_ns, px_bh, py_bh, pz_bh, skip_beginning_fraction=skip_beginning_fraction, skip_ending_fraction=skip_ending_fraction, 
                                                 plot=False, fit_acceleration=False)

    q = m_A[0]/m_B[0]
    nu = q_to_nu(q)
    # Set up the plot and plot the original data
    if plot:
        import matplotlib as mpl
        import scri.plotting

        try:
            if isinstance(w_m.metadata.alternative_names, list):
                SXS_BBH = "\n" + ", ".join(w_m.metadata.alternative_names)
            else:
                SXS_BBH = "\n" + w_m.metadata.alternative_names
        except:
            SXS_BBH = ""
        t_merger = w_m.max_norm_time() - 300.0
        t_ringdown = w_m.max_norm_time() + 100.0
        t_final = w_m.t[-1]
        delta_x = np.array([x_0 + v_0 * t_j for t_j in w_m.t])
        max_displacement = np.linalg.norm(delta_x, axis=1).max()
        max_d_color = min(1.0, 9 * max_displacement)
        LM_indices1 = [[2, 2], [2, 1], [3, 2], [3, 3], [4, 4]]
        LM_indices1 = [[ell, m] for ell, m in LM_indices1 if [ell, m] in w_m.LM.tolist()]
        indices1 = [(ell * (ell + 1) - 2 ** 2 + m) for ell, m in LM_indices1]
        LM_indices2 = [[2, 2], [2, 1], [3, 2], [3, 3], [4, 4]]
        LM_indices2 = [[ell, m] for ell, m in LM_indices2 if [ell, m] in w_m.LM.tolist()]
        indices2 = [(ell * (ell + 1) - 2 ** 2 + m) for ell, m in LM_indices2]
        fig1 = plt.figure(1, figsize=(10, 7))
        plt.gca().set_xscale("merger_zoom", t_merger=t_merger, t_ringdown=t_ringdown, t_final=t_final)
        lines1 = plt.semilogy(w_m.t, abs(w_m.data[:, indices1])*nu, alpha=0.35, lw=1.5)
        fig2 = plt.figure(2, figsize=(10, 7))
        plt.gca().set_xscale("merger_zoom", t_merger=t_merger, t_ringdown=t_ringdown, t_final=t_final)
        lines2 = plt.semilogy(w_m.t, abs(w_m.data[:, indices2])*nu, alpha=0.35, lw=1.5)

    # Import auxilliary waveforms if we are transforming Psi3, Psi2, Psi1, or Psi0. To apply a CoM
    # correction to these data types, information is required from all higher ranked Weyl scalars,
    # e.g. Psi2 requires information from Psi3 and Psi4.
    aux_waveforms = {}

    # Transform the mode data
    w_m = w_m.transform(space_translation=x_0, boost_velocity=v_0, **aux_waveforms)

    # Write the data to the new file
    if file_format.lower() == "nrar":
        path_to_new_waveform_h5 = re.sub(
            w_m.descriptor_string + "_",
            "",  # Remove 'rhOverM_', 'rMPsi4_', or whatever
            path_to_waveform_h5.replace(".h5", "_CoM.h5", 1),  # Add '_CoM' once
            flags=re.I,
        )  # Ignore case of 'psi4'/'Psi4', etc.
        write_to_h5(
            w_m,
            path_to_new_waveform_h5,
            file_write_mode=file_write_mode,
            attributes={"space_translation": x_0, "boost_velocity": v_0},
        )
    elif file_format.lower() == "rpxm":
        path_to_new_waveform_h5 = path_to_waveform_h5.replace(".h5", "_CoM.h5", 1)
        w_m.boost_velocity = v_0
        w_m.space_translation = x_0
        rotating_paired_xor_multishuffle_bzip2.save(w_m, path_to_new_waveform_h5)

    # Finish by plotting the new data and save to PDF
    if plot:
        plt.figure(1)
        for line, index, (ell, m) in zip(lines1, indices1, LM_indices1):
            plt.semilogy(w_m.t, abs(w_m.data[:, index])*nu, color=plt.getp(line, "color"), lw=1.5, label=f"({ell}, {m})")
        plt.figure(2)
        for line, index, (ell, m) in zip(lines2, indices2, LM_indices2):
            plt.semilogy(w_m.t, abs(w_m.data[:, index])*nu, color=plt.getp(line, "color"), lw=1.5, label=f"({ell}, {m})")
        for fig, num in [(fig1, 1), (fig2, 2)]:
            plt.figure(num)
            plt.axvline(t_merger, color="black", lw=2, alpha=0.125)
            plt.axvline(t_ringdown, color="black", lw=2, alpha=0.125)
            plt.xlabel(r"Coordinate time")
            plt.ylabel(r"Mode amplitudes")
            plt.xlim(0.0, w_m.t[-1])
            plt.ylim(1e-4, 1)
            plt.grid(axis="y")
            fig.text(
                0.5,
                0.91,
                name
                + SXS_BBH
                + "\n$x_0$ = [{}]".format(", ".join([str(tmp) for tmp in x_0]))
                + "\n$v_0$ = [{}]".format(", ".join([str(tmp) for tmp in v_0])),
                fontsize=8,
                ha="center",
            )
            #fig.text(
            #    0.004,
            #    0.004,
            #    str(max_displacement),
            #    fontsize=24,
            #    ha="left",
            #    va="bottom",
            #    bbox=dict(facecolor=mpl.cm.jet(max_d_color), alpha=max_d_color, boxstyle="square,pad=0.2"),
            #)
            plt.legend(loc="upper left", framealpha=0.9)
            fig.savefig(os.path.join(directory, "Modes_{}_{}.pdf".format(name, num)))
            plt.close()

    return w_m


def save_arrays_to_h5(wm, filename='rhOverM_Asymptotic_GeometricUnits.h5'):
    """
    Saves three arrays as columns in a dataset located in an HDF5 group.

    HDF5 structure:
    - Group: 'Extrapolated_N2.dir'
        - Dataset: 'Y_l2_m2.dat'

    Parameters:
    - multipolar waveform object from watpy
    - filename: name of the output .h5 file
    """
    # Write to HDF5 using h5py
    with h5py.File(filename, 'w') as h5file:
        group = h5file.create_group('Extrapolated_N2.dir')
        for mode in wm.modes:
            l, m = mode
            radii = wm.radii
            #print(radii,len(radii))
            M = wm.mass

            radii = radii[4:] ### THING TO TUNE CASE BY CASE
            #print(radii)

            ys = []
            yt = []
            for rad in radii:
                w22 = wm.get(l=l,m=m,r=rad)
                p22 = w22.p4
                ys.append(np.abs(p22))
                yt.append(np.unwrap(np.angle(p22)))
            
            time = w22.time_ret()/M

            method = 'Lousto' # 'polynomial'

            if method=='polynomial':
                #### with polynomial 1/R^K
                K = 3   ### THING TO TUNE CASE BY CASE
                yinf1 = radius_extrap_polynomial(ys,radii,K) # amplitude
                yinf2 = radius_extrap_polynomial(yt,radii,K) # phase
                p22_inf = yinf1 * np.exp(1j*yinf2) #yinf1 + 1j*yinf2

            elif method=='Lousto':
                #### Lousto method
                madm = float(dic['BHNS_ADM_mass'])
                ind = 2
                yyy = ys[ind] * np.exp(1j*yt[ind])
                p22_inf = radius_extrap(time,yyy, radii[ind], l=2, m=2, m0=madm)

            else:
                print('Method not implemented')

            f0 = w22.prop['init.frequency']
            fcut = 2 * f0 / max(1,abs(w22.prop['mmode']))
            dt = w22.time[1] - w22.time[0]
            h22_inf = 1. * fixed_freq_int_2( 1. * p22_inf, fcut, dt = dt)

            array1, array2, array3 = map(np.asarray, (time, np.real(h22_inf), np.imag(h22_inf)))
            data = np.column_stack((array1, array2, array3))
            group.create_dataset('Y_l'+str(l)+'_m'+str(m)+'.dat', data=data)


def main():
    ## Read args
    parser = argparse.ArgumentParser(description='Evolution Diagnostics, e.g. `python check_ev.py -i id_name -c cluster_name`, either in ARA, DRACO or HAWK')
    # add arguments
    parser.add_argument('-i', '--input', type=str, required=True, help='simulation name')
    parser.add_argument('-e', '--evo', type=str, required=False, help='evolution name (optional)')
    parser.add_argument('-c', '--cluster', type=str, required=True, help='cluster')

    # parse the arguments
    args = parser.parse_args()

    # access the arguments
    simname = args.input
    evol = args.evo
    cluster = args.cluster

    return simname, evol, cluster

if __name__ == '__main__':
    simname, evol, cluster = main()


    if cluster=='ARA':
        basedir = '/beegfs/mo63poy/BHNS_Elliptica'
        id_exe = '/home/mo63poy/BHNS_Elliptica/Elliptica/Exe/elliptica'
        ev_path = '/home/mo63poy/BNS_BAM/BAM'
    elif cluster=='DRACO':
        basedir = '/work/mo63poy/BHNS_INITIAL_DATA'
        id_exe = '/home/mo63poy/Elliptica/Exe/elliptica'
        ev_path = '/home/mo63poy/BAM'
    elif cluster=='HAWK': # modify HERE!
        batchsys = 'pbs'
        basedir = '/zhome/academic/HLRS/xuj/xujapigo/BHNS_INITIAL_DATA'
        id_exe = '/zhome/academic/HLRS/xuj/xujapigo/Elliptica/Exe/elliptica'
        bam_path = '/zhome/academic/HLRS/xuj/xujapigo/BAM'
        workspace = '/lustre/hpe/ws10/ws10.3/ws/xujapigo-bhns'

        evo_path = os.path.join(workspace,simname) # where the simulation is running (inside our hawk workspace)
    else:
        print('Error: cluster name unknown')

    simpath = os.path.join(basedir,simname)

    ## Do evolution(s)
    # Initialize ID object
    bhns_id = Initial_Data(path=simpath,params=None, id_exe=id_exe)
    id_output = bhns_id.ou
    dic = id_output.id_dic
    lmax2 = 6
    resolution, lmax, flux = ev_folder_parser(evol)
    bhns_ev = Evolution(path=os.path.join(simpath,evol), ev_path=ev_path, initial_data=bhns_id, resolution=int(resolution), lmax=int(lmax), lmax2=lmax2,flux=flux)
    ev_output = bhns_ev.ou
    #ev_output.plot_moving_puncture()

    px_ns, py_ns, px_bh, py_bh, time = ev_output.extract_objects_tracks()
    mb = float(dic['NS_baryonic_mass_current'])*np.ones(len(px_ns))
    mbh = float(dic['BH_Christodoulou_mass_current'])*np.ones(len(px_bh))
    
    core_out = bhns_ev.core_out
    wm = bhns_ev.get_core_wm_object(ignore_negative_m=False)
    #print(wm.modes)


    #plt.plot(time,np.real(p22_inf)*400,label='psi4 (rescaled)')
    #plt.plot(time,np.real(h22_inf),label='h')
    #plt.legend()
    #plt.show()

    '''
    # check phase
    yphi = np.unwrap(np.angle(p22_inf))

    pphi_last = np.unwrap(np.angle(ys[-1] * np.exp(1j*yt[-1])))-5*np.pi
    pphi_first = np.unwrap(np.angle(ys[0] * np.exp(1j*yt[0])))-1*np.pi

    plt.plot(time,pphi_last,label='last rad')
    plt.plot(time,pphi_first,label='first rad')
    plt.plot(time,yphi,label='extr')
    plt.legend()
    plt.show()

    plt.plot(time,np.abs(yphi-pphi_last),label='inf-last')
    plt.plot(time,np.abs(yphi-pphi_first),label='inf-first')
    plt.yscale('log')
    plt.xlim([0,time[-1]])
    plt.legend()
    plt.show()
    '''

    # the magic happens below:
    save_arrays_to_h5(wm)
    w_m = remove_avg_com_motion(simname,skip_beginning_fraction=0.01, skip_ending_fraction=0.10,plot=True,file_write_mode="w",ev_output=ev_output, m_A=mbh,m_B=mb,file_format="NRAR") 

    # save back to CoRe format
    new_name = simname+"_CoM"
    if os.path.exists(new_name):
        print(simname+' dircetory exists')
    else:
        os.mkdir(new_name)


    Mbh, Mg, M = id_output.get_msun_masses()
    q = Mbh/Mg
    nu = q_to_nu(q)
    waveforms = {}
    with h5py.File('rhOverM_Asymptotic_GeometricUnits_CoM.h5', 'r') as h5file:
        group = h5file['Extrapolated_N2.dir']

        for name, dataset in group.items():
            # Skip 0D datasets (like History.txt)
            if len(dataset.shape) == 0:
                continue

            data = dataset[:]
            time, h_real, h_imag = data.T
            h = h_real + 1j * h_imag
            
            # Extract (l, m) from the name, e.g., 'Y_l2_m2.dat'
            parts = name.replace('.dat', '').split('_')
            l = int(parts[1][1:])  
            m = int(parts[2][1:])  
            
            waveforms[(l, m)] = {'time': time, 'h': h}

    for l,m in list(waveforms.keys()):
        wv = wave(path=new_name,code='core',mass=M,f0=wm.f0)
        wv.prop['var'] = 'h'
        wv.prop['lmode']           = l
        wv.prop['mmode']           = m
        wv.prop['detector.radius'] = -1
        wv.h = waveforms[(l,m)]['h']*nu
        wv.time = waveforms[(l,m)]['time']
        wv.write_to_txt(var='h',path=new_name)

    os.system('rm *.h5')
