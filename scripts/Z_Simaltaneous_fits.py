#!/usr/bin/env python
import uproot
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
from json import load,dump
import argparse
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from scipy.interpolate import CubicSpline
from matplotlib.colors import ListedColormap
from iminuit import Minuit
import sqlite3
from matplotlib.patches import Ellipse
from datetime import datetime
#where is the main function? :(

#This allows for me to pass arguments in to the terminal to change important paramters without changing the code
parser = argparse.ArgumentParser(description='some string')
#Run with --Callibration="TRUE" to use callibration
parser.add_argument('--Calibration',default="FALSE",type=str)
parser.add_argument('--Smear_error',default="NONE",type=str,choices=["PLUS", "MINUS", "NONE"])
parser.add_argument('--C_ratio_error',default="NONE",type=str,choices=["PLUS", "MINUS", "NONE"])
args = parser.parse_args()

DATAIR="/storage/epp2/phshgg/Public/MPhysProject_2025_2026/tuples/1/"
with uproot.open(f"{DATAIR}/DecayTree__Z__Z__d13600GeV_24c4.root:DecayTree") as t:

    #momentasim = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    simdatam=t.arrays(["mum_eta","mum_phi","mum_pt","mup_eta" ,"mup_phi" ,"mup_pt","true_boson_mass"],library="np")  #i ahve ztrue and i ahve z reconstructed simulation
    
    tmass=simdatam["true_boson_mass"]


with uproot.open(f"{DATAIR}/DecayTree__Z__DATA__d13600GeV_24c4.root:DecayTree") as tt:

    #momentasim = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    datam=tt.arrays(["mum_eta","mum_phi","mum_pt","mup_eta" ,"mup_phi" ,"mup_pt"],library="np")

if (args.Calibration).lower() == "true":
    try:
        with open("Calibration_output.json",) as InputFile:
            json_data=load(InputFile)
        C_rat,C_err =  json_data["C_ratio"]
        Smear_factor,Smear_factor_error = json_data["Smear_factor"]
        print(f'Woo we imported C_ratio to be {C_rat} ± {C_err}')
        print(f'Woo we imported smear factor to be {Smear_factor} ± {Smear_factor_error}')
        fname = f"mass-width_values_and_error.txt"
        if (args.Smear_error).lower() == "plus":
            Smear_factor=Smear_factor+Smear_factor_error
            fname = f"smear_plus.txt"
        elif (args.Smear_error).lower() == "minus":
            Smear_factor=Smear_factor-Smear_factor_error
            fname = f"smear_minus.txt"
        if (args.C_ratio_error).lower() == "plus":
            C_rat=C_rat+C_err
            fname = f"C_rat_plus.txt"
        elif (args.C_ratio_error).lower() == "minus":
            C_rat=C_rat-C_err
            fname = f"C_rat_minus.txt"
        
        rng = np.random.default_rng(seed=10)
        Norm_rand = rng.normal(0,1,size=len(tmass))
        alpha = 1 - C_rat
        calibration_factor = 1+alpha+Norm_rand*Smear_factor

        print(f'Woo the correction is now  {calibration_factor}')

    except FileNotFoundError:
        print("Please run the script Upsilon.py first to get a calibration value")
        #This is where I would return to stop the script if we were inside a main func to avoid error
        calibration_factor=1
else:
    calibration_factor = 1
    fname = f"useless.txt"



def conaeq(reconmass,cal):
    mpe=reconmass["mup_eta"]
    mppt=reconmass["mup_pt"] 
    mpphi=reconmass["mup_phi"] 
    mme=reconmass["mum_eta"] 
    mmpt=reconmass["mum_pt"] 
    mmphi=reconmass["mum_phi"]
    # calibration being applied here
    mpx=mppt*np.cos(mpphi)*cal
    mmx=mmpt*np.cos(mmphi)*cal
    mpy=mppt*np.sin(mpphi)*cal
    mmy=mmpt*np.sin(mmphi)*cal
    mpz=mppt*np.sinh(mpe)*cal
    mmz=mmpt*np.sinh(mme)*cal

    m=[]
    Ep=[]
    Em=[]
    M_mass=105.658*10**(-3)
    EP=np.sqrt(M_mass**2 +mpx**2 +mpy**2+mpz**2)   
    EN=np.sqrt(M_mass**2 +mmx**2 +mmy**2+mmz**2)
    m=np.sqrt((EP+EN)**2-((mpx+mmx)**2+(mpy+mmy)**2+(mpz+mmz)**2))
    
    return m

def fiteq(x,a,b,m,w,scale):
    gamma=np.sqrt(m**2*(m**2+w**2))
    k=(2*np.sqrt(2)*m*gamma*w)/(np.pi*np.sqrt(m**2+gamma))
    f= ((k/((x**2-m**2)**2+(m**2)*w**2))*scale  + a*np.exp(b*x))  # i can just do this manually
    return f

def sim_fits(tmass,simdatam,datam,calibration_factor,use_diagram):
    c=calibration_factor #error in calibration
    dataHist,databinn,_d=plt.hist(conaeq(datam,1), bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-data-reconstructed",linewidth=1) #the data recosntuction data #for recosntructed simulation 
    trmassHist,binn,_t=plt.hist(tmass, bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-true",density=True,linewidth=1) #for true mass
    centers=0.5*(binn[1:]+binn[:-1])
    w=3
    m=91.1876
    fitParam,_tt = curve_fit(fiteq,centers,trmassHist,p0=[5,-0.7,m,w,0.4],bounds=([0.0,-1.0,60,0,0],[100.0,0,120,10,1]),maxfev=10000)# this uses the true mass

    # i want 9 differnet parameters

    ratio905_1=fiteq(tmass,fitParam[0],fitParam[1],90.5,1,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights905_1,binnweights905_1 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio905_1)
    ratio905_2=fiteq(tmass,fitParam[0],fitParam[1],90.5,2,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights905_2,binnweights905_2 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio905_2)
    ratio905_3=fiteq(tmass,fitParam[0],fitParam[1],90.5,3,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights905_3,binnweights905_3 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio905_3)

    ratio91_2=fiteq(tmass,fitParam[0],fitParam[1],91,2,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights91_2,binnweights91_2 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio91_2)
    ratio91_1=fiteq(tmass,fitParam[0],fitParam[1],91,1,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights91_1,binnweights91_1 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio91_1)
    ratio91_3=fiteq(tmass,fitParam[0],fitParam[1],91,3,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights91_3,binnweights91_3 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio91_3)

    ratio915_1=fiteq(tmass,fitParam[0],fitParam[1],91.5,1,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights915_1,binnweights915_1 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio915_1)
    ratio915_2=fiteq(tmass,fitParam[0],fitParam[1],91.5,2,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights915_2,binnweights915_2 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio915_2)
    ratio915_3=fiteq(tmass,fitParam[0],fitParam[1],91.5,3,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights915_3,binnweights915_3 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio915_3)

    sim_scale_factor905_1= np.sum(dataHist) / np.sum(simmassHistweights905_1)
    simHist_scaled905_1 = simmassHistweights905_1 * sim_scale_factor905_1
    sim_scale_factor905_2= np.sum(dataHist) / np.sum(simmassHistweights905_2)
    simHist_scaled905_2 = simmassHistweights905_2 * sim_scale_factor905_2
    sim_scale_factor905_3= np.sum(dataHist) / np.sum(simmassHistweights905_3)
    simHist_scaled905_3 = simmassHistweights905_3 * sim_scale_factor905_3
    
    sim_scale_factor91_1= np.sum(dataHist) / np.sum(simmassHistweights91_1)
    simHist_scaled91_1 = simmassHistweights91_1 * sim_scale_factor91_1
    sim_scale_factor91_2= np.sum(dataHist) / np.sum(simmassHistweights91_2)
    simHist_scaled91_2 = simmassHistweights91_2 * sim_scale_factor91_2
    sim_scale_factor91_3= np.sum(dataHist) / np.sum(simmassHistweights91_3)
    simHist_scaled91_3 = simmassHistweights91_3 * sim_scale_factor91_3

    sim_scale_factor915_1= np.sum(dataHist) / np.sum(simmassHistweights915_1)
    simHist_scaled915_1 = simmassHistweights915_1 * sim_scale_factor915_1
    sim_scale_factor915_2= np.sum(dataHist) / np.sum(simmassHistweights915_2)
    simHist_scaled915_2 = simmassHistweights915_2 * sim_scale_factor915_2
    sim_scale_factor915_3= np.sum(dataHist) / np.sum(simmassHistweights915_3)
    simHist_scaled915_3 = simmassHistweights915_3 * sim_scale_factor915_3



    plt.plot(centers, simHist_scaled905_1, '-', linewidth=2, label='scaled-weighted 90.5_1 simulation')
    plt.plot(centers, simHist_scaled91_1, '-', linewidth=2, label='scaled-weighted 91_1 simulation')
    plt.plot(centers, simHist_scaled915_1, '-', linewidth=2, label='scaled-weighted 91.5_1 simulation')

    plt.plot(centers, simHist_scaled905_2, '-', linewidth=2, label='scaled-weighted 90.5_2 simulation')
    plt.plot(centers, simHist_scaled91_2, '-', linewidth=2, label='scaled-weighted 91_2 simulation')
    plt.plot(centers, simHist_scaled915_2, '-', linewidth=2, label='scaled-weighted 91.5_2 simulation')

    plt.plot(centers, simHist_scaled905_3, '-', linewidth=2, label='scaled-weighted 90.5_3 simulation')
    plt.plot(centers, simHist_scaled91_3, '-', linewidth=2, label='scaled-weighted 91_3 simulation')
    plt.plot(centers, simHist_scaled915_3, '-', linewidth=2, label='scaled-weighted 91.5_3 simulation')
    plt.title(f"Zreconstructed weight sim ({'real' if use_diagram==True else 'magnet' if use_diagram== 'pos_dipole' else 'neg-magnet' if use_diagram== 'neg_dipole' else 'dont-use'})")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig(f"transient/Zreconstructed weight sim-combined_mess ({'real' if use_diagram==True else 'pos-magnet' if use_diagram== 'pos_dipole' else 'neg-magnet' if use_diagram== 'neg_dipole' else 'dont-use'}).pdf")
    plt.clf()
    plt.figure()
    chi905_1=np.sum(((dataHist-simHist_scaled905_1)**2)/dataHist)
    chi905_2=np.sum(((dataHist-simHist_scaled905_2)**2)/dataHist)
    chi905_3=np.sum(((dataHist-simHist_scaled905_3)**2)/dataHist)

    chi91_1=np.sum(((dataHist-simHist_scaled91_1)**2)/dataHist)
    chi91_2=np.sum(((dataHist-simHist_scaled91_2)**2)/dataHist)
    chi91_3=np.sum(((dataHist-simHist_scaled91_3)**2)/dataHist)

    chi915_1=np.sum(((dataHist-simHist_scaled915_1)**2)/dataHist)
    chi915_2=np.sum(((dataHist-simHist_scaled915_2)**2)/dataHist)
    chi915_3=np.sum(((dataHist-simHist_scaled915_3)**2)/dataHist)

    x=np.array([1,2,3,1,2,3,1,2,3])
    y=np.array([90.5,90.5,90.5,91,91,91,91.5,91.5,91.5])
    chi_values=np.array([chi905_1,chi905_2,chi905_3,chi91_1,chi91_2,chi91_3,chi915_1,chi915_2,chi915_3])

    t_width=np.array([1,2,3,1,2,3,1,2,3])
    t_mass=np.array([90.5,90.5,90.5,91,91,91,91.5,91.5,91.5])
    t_data=np.array([simHist_scaled905_1,simHist_scaled905_2,simHist_scaled905_3,simHist_scaled91_1,simHist_scaled91_2,simHist_scaled91_3,simHist_scaled915_1,simHist_scaled915_2,simHist_scaled915_3])
    coords = np.column_stack([t_mass,t_width])
    templates=np.array(t_data)

    def interpolate_chi_func(mass_value,width_value):
        interpolate_template=griddata(coords,templates,(mass_value,width_value),method='linear')
        return np.sum(((dataHist-interpolate_template)**2)/interpolate_template) 
    
#///////////////////////////////////////////////////////////////////////////////// calulcations for min mass and width
        # the Minuit doesnt seem to work without passign a fucntion through it 
    m = Minuit(interpolate_chi_func, 90.6, 2.4)
    m.migrad()
    m.hesse()
    mass_result=m.values["mass_value"]
    width_result=m.values["width_value"]
    mass_error=m.errors["mass_value"]
    width_error=m.errors["width_value"]
    covariance_matrix=m.covariance
    corelation_coefficient=(m.covariance["mass_value", "width_value"]) /(width_error*mass_error)
    print(f" the mass is {mass_result} ± {mass_error}" )
    print("Chi^2 mini:", m.fval)
    print(f"the width is {width_result} ± {width_error}")
    print(covariance_matrix)
    print(f"corelation is {corelation_coefficient}")
    ndf=len(binnweights905_1)-2
    chimin=np.min(chi_values)
    print(chimin)
    if fname != "dont_use":
        with open(fname, "w") as f:
            f.write(f"{mass_result}\n")
            f.write(f"{width_result}")
            if fname==f"mass-width_values_and_error.txt":
                f.write(f"\n{mass_error}\n")
                f.write(f"{width_error}\n")
                f.write(f"{chimin}\n")
                f.write(f"{ndf}\n")
                f.write(f"{corelation_coefficient}")

        f.close()

    if use_diagram != True:
        with open(use_diagram, "w") as f:
            f.write(f"{mass_result}\n")
            f.write(f"{width_result}")
            f.write(f"\n{mass_error}\n")
            f.write(f"{width_error}\n")
            f.write(f"{chimin}\n")
            f.write(f"{ndf}\n")
            f.write(f"{corelation_coefficient}")



    # plotting the stack graph with similtanous fits
    #the top half.

    dataerrors=np.sqrt(dataHist)

    plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "STIX"],    # this is making the graphs look like PRL
    "mathtext.fontset": "stix",
    "font.size": 12,
})
    plt.rcParams.update({
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
    "axes.linewidth": 0.8,
})
    plt.figure()
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 6),
                                   gridspec_kw={'hspace': 0.05})
    
    #top
   # ax1.scatter(centers,dataHist ,label='data ',color="black",s=10)
    ax1.errorbar(centers, dataHist, yerr=dataerrors,label='Data ',color="black",fmt=".",markersize=2.5 )
    ax1.step(centers, simHist_scaled91_2, '-', linewidth=2,where='mid', label='91 Gev mass 2 Gev width template')
    ax1.set_ylabel("Frequency")
    ax1.legend(loc='upper left',frameon=True, fontsize=8)
    ax1.set_ylim(bottom=0)
    # bottom
    ratio_91_3_91_2=simHist_scaled91_3/simHist_scaled91_2
    ratio_91_2_91_2=simHist_scaled91_2/simHist_scaled91_2
    ratio_915_91_2=simHist_scaled915_2/simHist_scaled91_2

    ratio_data_91_2=dataHist/simHist_scaled91_2
    ratioerror=ratio_data_91_2*dataerrors/dataHist
    #ax2.scatter(centers,ratio_data_91,label='data ',color="black",s=10)
    ax2.errorbar(centers, ratio_data_91_2, yerr=ratioerror,label='Data/91_2 Gev template ',color="black",fmt=".",markersize=2.5 )
    ax2.step(centers, ratio_91_3_91_2, '-', linewidth=2,where='mid', label='91_3/91_2 Gev template ratio')
    ax2.step(centers, ratio_91_2_91_2, '-', linewidth=2,where='mid', label='91_2/91_2_2 Gev template ratio')
    ax2.step(centers, ratio_915_91_2, '-', linewidth=2,where='mid', label='91.5_2/91_2 Gev template ratio')
    ax2.set_ylabel("Ratio")
    ax2.set_xlabel("Dimuon mass GeV")
    ax2.set_ylim(bottom=0.8)
    ax2.legend(loc='upper left',frameon=True, fontsize=8)
    plt.savefig(f"transient/Z-stack_similtaneous ({'real' if use_diagram==True else 'pos-magnet' if use_diagram== 'pos_dipole' else 'neg-magnet' if use_diagram== 'neg_dipole' else ''}).pdf")
                
    print(chi_values)

    # graph elipse of mass width
    plt.figure()
    eigenvalues , eigenvectors  = np.linalg.eigh(covariance_matrix)
    order = eigenvalues.argsort()[::-1] 
    eigenvalues = eigenvalues[order]
    eigenvectors=eigenvectors[:, order]  #making the bigegst eignecalue and coresponding vector be first in order
    major_axis_vector = eigenvectors[:,0]
    vx, vy = major_axis_vector[0], major_axis_vector[1]
    theta = np.degrees(np.arctan2(vy, vx))
    major_axis_length = 2 * np.sqrt(eigenvalues[0]) 
    minor_axis_length = 2 * np.sqrt(eigenvalues[1])
    fig, ax = plt.subplots()
    ellipse = Ellipse( (mass_result, width_result), width=major_axis_length, height=minor_axis_length, angle=theta, edgecolor='red', facecolor='none', linewidth=2)
    ax.add_patch(ellipse)
    ax.scatter(mass_result, width_result, color='blue', label='Measurement')
    ax.set_xlabel("Mass") 
    ax.set_ylabel("Width") 
    ax.set_title("Mass Width error Ellipse") 
    ax.legend(loc="upper left") 
    ax.set_xlim(91.13, 91.155) 
    ax.set_ylim(1.965,2.035) 
    plt.savefig(f"transient/Z-Error-Ellipse ({'real' if use_diagram==True else 'pos-magnet' if use_diagram== 'pos_dipole' else 'neg-magnet' if use_diagram== 'neg_dipole' else ''}).pdf")

sim_fits(tmass,simdatam,datam,calibration_factor,True)


# this is to test for magnetic monople 
#real means standard ie can use, watch out that adign the erros will chnage this so real is only trustworthy using general calibration and the run all script
# make sure that calibration ="true" is used and nothign else or if using other arguenmtns make sure jsut calibration=ture is used last
if fname == f"mass-width_values_and_error.txt":
    with uproot.open(f"{DATAIR}/DecayTree__Z__DATA__d13600GeV_24c4.root:DecayTree") as tt:
        magnet_pol=tt.arrays(["yearpol"],library="np")
    magnet_pol=np.array(magnet_pol["yearpol"])
    print(magnet_pol)
    pos_mask= magnet_pol > 0
    neg_mask= magnet_pol < 0
    pos_datam= {k: v[pos_mask] for k  , v in datam.items()}
    neg_datam= {k: v[neg_mask] for k, v in datam.items()}

    # with uproot.open(f"{DATAIR}/DecayTree__Z__Z__d13600GeV_24c4.root:DecayTree") as t:
    #     magnet_pol=t.arrays(["yearpol"],library="np")
    # magnet_pol=np.array(magnet_pol["yearpol"])
    # print(magnet_pol)
    # pos_mask= magnet_pol > 0
    # neg_mask= magnet_pol < 0
    # simpos_datam= {k: v[pos_mask] for k  , v in simdatam.items()}
    # simneg_datam= {l: w[neg_mask] for l, w in simdatam.items()}
    # truemass_pos=tmass[pos_mask]
    # truemass_neg=tmass[neg_mask]

    # calibration_slice_plus=calibration_factor[pos_mask]
    # calibration_slice_minus=calibration_factor[neg_mask]
    magnet="pos_dipole" 
    fname="dont_use"

    sim_fits(tmass,simdatam,pos_datam,calibration_factor,magnet)
    magnet="neg_dipole" 
    sim_fits(tmass,simdatam,neg_datam,calibration_factor,magnet)
# can add  more here qutie easily 