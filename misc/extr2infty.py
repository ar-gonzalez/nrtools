import argparse
import os
from nrtools.initialdata.initialdata import *
from nrtools.evolution.evolution import *
from nrtools.utils.utils import get_rad
from watpy.wave.wave import wfile_get_mass, wfile_get_detrad
from watpy.wave.gwutils import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['text.usetex']= True
matplotlib.rcParams['font.serif']= 'Palatino' 
matplotlib.rcParams['font.size']= 15 #28

lvl = {
    6: '#9ebcda',
    8: '#8c96c6',
    10: '#8856a7',
    11: '#810f7c'
}

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='compare different refinement levels of a BAM simulation')

    # Add arguments
    parser.add_argument('-i','--input', help='ID name')
    parser.add_argument('-e','--evolution', help='Evolution name')
    parser.add_argument('-c','--cluster', help='Cluster where it is stored')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Access the values of the arguments
    simname = args.input
    evo = args.evolution
    cluster = args.cluster

    return simname, evo, cluster

if __name__ == '__main__':
    simname, evo, cluster = main()

    if cluster=='ARA':
        basedir = '/beegfs/mo63poy/BHNS_Elliptica'
        id_exe = '/home/mo63poy/BHNS_Elliptica/Elliptica/Exe/elliptica'
        ev_path = '/home/mo63poy/BNS_BAM/BAM'
    elif cluster=='DRACO':
        basedir = '/home/mo63poy/BHNS_INITIAL_DATA'
        id_exe = '/home/mo63poy/Elliptica/Exe/elliptica'
        ev_path = '/home/mo63poy/BAM'
    elif cluster=='home': # my PC
        basedir = '/home/agonzalez/Documents/PhD/NR/BHNS_evs/evolutions'
    else:
        print('Error: cluster name unknown')

    simpath = os.path.join(basedir,simname)
    
    psi_ds = []
    psii_ds = []
    t_ds = []

    # Get Data
    radii = []

    resolution, lmax, flux = ev_folder_parser(evo)
    core_out = os.path.join(os.path.join(simpath,evo),'CoReDB')
    file_list = sorted([i for i in os.listdir(core_out) if i.startswith('Rpsi4_l2_m1_r')])
    radius = int(wfile_get_detrad(os.path.join(core_out,file_list[-1])))

    print("Extraction radius: ",radius)

    file = [i for i in file_list if i.endswith(str(radius)+'.txt')][0]
    uM, rpsi4, ipsi4, momg, aM, phi1, t = np.loadtxt(fname=os.path.join(core_out,file), comments='#', usecols=(0,1,2,3,4,5,6), unpack=True)
    #psi_ds.append(rpsi4)
    #psii_ds.append(ipsi4)
    #t_ds.append(uM)
    phi1 = -1*np.unwrap(np.angle(rpsi4 + 1j*ipsi4))
    amp = np.abs(rpsi4 + 1j*ipsi4)
    tmrg = uM[np.argmax(amp)]


   
    ############################
    # Extrapolation w/polynomial of order K
    ############################
    file_list = [i for i in os.listdir(core_out) if i.startswith('Rpsi4_l2_m1_r')]
    file_list = sorted(file_list)
    ys = []
    yt = []
    ya = []
    yp = []
    radii = []
    for i,file1 in enumerate(file_list):
        rad1 = get_rad(file1) 
        uy, rp4, ip4, _, ay, py, ty = np.loadtxt(fname=os.path.join(core_out,file1), comments='#', usecols=(0,1,2,3,4,5,6), unpack=True)
        if rad1>200:
            ys.append(rp4) 
            yt.append(ip4)
            ya.append(ay)
            yp.append(py)
            radii.append(rad1)
        if i==0:
            aay = ay
            ppy=py
            tty = uy
            rad2 = rad1

    K = 3
    yinf1 = radius_extrap_polynomial(ya,radii,K)
    yinf2 = radius_extrap_polynomial(yp,radii,K)
    yinf = ya*np.exp(yinf2)#yinf1 + 1j*yinf2
    yphi = yinf2#np.unwrap(np.angle(yinf))
    yamp = yinf1#np.abs(yinf)

    yphi= -1*yphi - np.pi*2
    ppy = ppy #- np.pi*8
    '''
    plt.plot(uM,phi1,label='phi1')
    plt.plot(uy,yphi,label='yphi')
    plt.plot(tty,ppy,label='ppy')
    plt.legend()
    plt.show()
    '''
    
    ppy = np.interp(uy, tty, ppy)
    dp_r = yphi-ppy # (wrt to lowest radius)
    yphi = np.interp(uM, uy, yphi)
    dp_rf = yphi-phi1 # uncertainty due to finite radius (wrt highest radius)
    print('Highest radius: ',radii[-1], 'lowest: ', rad2)
    print(radii)
    '''
    plt.plot(uy,yamp,label='extr')
    plt.xlim([0,uy[-1]])
    plt.plot(tty,aay,label='first rad')
    plt.plot(uy,ay,label='last rad')
    plt.legend()
    plt.show()
    '''
    ######################
    # Lousto method
    ######################
    bhns_id = Initial_Data(path=simpath,params=None, id_exe=id_exe)
    id_output = bhns_id.ou
    dic = id_output.id_dic
    madm = float(dic['BHNS_ADM_mass'])
    tty=uy
    yinf = radius_extrap(tty, ys[-1]+1j*yt[-1], radii[-1], l=2, m=1, m0=madm)

    plt.plot(tty,np.abs(yinf),label='extr')
    plt.plot(tty,aay,'--',label='first rad')
    plt.plot(uy,ay,'--',label='last rad')
    plt.legend()
    plt.show()

    lousto_phi = -1*np.unwrap(np.angle(yinf))

    plt.plot(tty,lousto_phi,label='extr')
    plt.plot(tty,ppy,'--',label='first rad')
    plt.plot(uy,py,'--',label='last rad')
    plt.legend()
    plt.show()

    # get strain
    _, _, mtot = id_output.get_msun_masses()
    _, Momg22 = id_output.get_gw_freqs()
    f0 = Momg22 / (2*np.pi) / mtot

    mmode = 1
    fcut = 2 * f0 / max(1,abs(mmode))
    dt = tty[1] - tty[0]
    strain = fixed_freq_int_2(yinf, fcut, dt = dt)

    bhns_ev = Evolution(path=os.path.join(simpath,evo), ev_path=ev_path, initial_data=bhns_id, resolution=int(resolution), lmax=int(lmax), lmax2=6,flux=flux)
    #wm = bhns_ev.get_core_wm_object()
    #w21_300 = wm.get(l=2,m=1,r=300)
    #w21_800 = wm.get(l=2,m=1,r=800)
    uM, rpsi4, ipsi4, momg, aM, phi1, t = np.loadtxt(fname=os.path.join(core_out,'Rh_l2_m1_r00800.txt'), comments='#', usecols=(0,1,2,3,4,5,6), unpack=True)
    plt.plot(tty,np.abs(strain),label='extr')
    plt.plot(tty,np.real(strain))
    #plt.plot(tty,np.abs(w21_300.h),'--',label='r=300')
    #plt.plot(uM,aM,'--',label='r=800')
    #plt.legend()
    plt.show()


    quit()
    ############################
    # total error budget
    ############################

    rinf_color = '#d95f02'
    rich_color = '#7570b3'
    bkg_color = '#1b9e77'

    fig = plt.subplots(figsize=(7,4))
    plt.plot(tty,np.log10(np.abs(dp_r)), label=r'$R(\infty)$ - $R('+str(int(rad2))+'), K=3$', color=rinf_color)
    plt.plot(uM,np.log10(np.abs(dp_rf)), label=r'$R(\infty)$ - $R('+str(int(rad1))+'), K=3$', color=rinf_color,linestyle='--')
    #plt.fill_betweenx(y=np.arange(-20.,20.,0.1),x1=rich_mrg,x2=tmrg2,color='gray',alpha=0.4)
    plt.xlim([0,uM[-1]])
    plt.ylim([-4,1.5])
    plt.xlabel(r'$u/M$')
    plt.ylabel(r'$\log_{10}|\Delta\phi_{22}|$')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()
