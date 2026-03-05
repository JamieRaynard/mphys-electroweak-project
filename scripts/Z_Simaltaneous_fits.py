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
    datam=tt.arrays(["mum_eta","mum_phi","mum_pt","mup_eta" ,"mup_phi" ,"mup_pt","yearpol","V_PT"],library="np")

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
        
        calibration_factor = (C_rat,Smear_factor)

        print(f'Woo the correction is now  {calibration_factor}')

    except FileNotFoundError:
        print("Please run the script Upsilon.py first to get a calibration value")
        #This is where I would return to stop the script if we were inside a main func to avoid error
        calibration_factor=None

else:
    calibration_factor = None
    fname = f"useless.txt"



def conaeq(reconmass,cal):
    mpe=reconmass["mup_eta"]
    mppt=reconmass["mup_pt"] 
    mpphi=reconmass["mup_phi"] 
    mme=reconmass["mum_eta"] 
    mmpt=reconmass["mum_pt"] 
    mmphi=reconmass["mum_phi"]
    # calibration being applied here
    mpx=mppt*np.cos(mpphi)
    mmx=mmpt*np.cos(mmphi)
    mpy=mppt*np.sin(mpphi)
    mmy=mmpt*np.sin(mmphi)
    mpz=mppt*np.sinh(mpe)
    mmz=mmpt*np.sinh(mme)

    mup_P = np.array([mpx,mpy,mpz])
    mum_P = np.array([mmx,mmy,mmz])

    if cal:
        #x = cal[1]/(np.mean([mup_P,mum_P]))
        #cal[1] = x
        #mup_P = 1/(cal[0]*(1/mup_P)+cal[1])
        #mum_P = -1/(cal[0]*(-1/mum_P)+cal[1])
        rng_p = np.random.default_rng(seed=10)
        rng_m = np.random.default_rng(seed=11)
        Norm_rand_p = rng_p.normal(0,1,size=len(mpx))
        Norm_rand_m = rng_m.normal(0,1,size=len(mpx))
        mup_P = mup_P*cal[0]*(1+mup_P*Norm_rand_p*cal[1])
        mum_P = mum_P*cal[0]*(1+mum_P*Norm_rand_m*cal[1])
        mpx,mpy,mpz = mup_P
        mmx,mmy,mmz = mum_P

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

def sim_fits(tmass,simdatam,datam,calibration_factor,use_diagram,bin_number):
    c=calibration_factor #error in calibration
    bin_space=np.linspace(80.0,100.0,bin_number)
    dataHist,databinn,_d=plt.hist(conaeq(datam,None), bins=bin_space, histtype="step",label="Z-data-reconstructed",linewidth=1) #the data recosntuction data #for recosntructed simulation 
    trmassHist,binn,_t=plt.hist(tmass, bins=bin_space, histtype="step",label="Z-true",density=True,linewidth=1) #for true mass
    centers=0.5*(binn[1:]+binn[:-1])
    w=3
    m=91.1876
    fitParam,_tt = curve_fit(fiteq,centers,trmassHist,p0=[5,-0.7,m,w,0.4],bounds=([0.0,-1.0,60,0,0],[100.0,0,120,10,1]),maxfev=10000)# this uses the true mass

    plt.figure()#-------------------------true mass function graph
    trmassHist,binn,_t=plt.hist(tmass, bins=bin_space, histtype="step",label="Z-true",density=True,linewidth=1)
    fitParam,_tt = curve_fit(fiteq,centers,trmassHist,p0=[5,-0.7,m,w,0.4],bounds=([0.0,-1.0,60,0,0],[100.0,0,120,10,1]),maxfev=10000)
    x_fit = np.linspace(80, 100, 500)
    y_fit=fiteq(x_fit,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    plt.plot(x_fit,y_fit)
    plt.xlabel("Dimuon Mass GeV")
    plt.ylabel("Events")
    plt.savefig(f"transient/True_mass_fit ({'real' if use_diagram==True else use_diagram}).pdf")

    # i want 9 differnet parameters

    ratio905_1=fiteq(tmass,fitParam[0],fitParam[1],90.5,1,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])#ratio--
    simmassHistweights905_1,binnweights905_1 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio905_1)#----------------------------widghted sim ie tempalte
    ratio905_2=fiteq(tmass,fitParam[0],fitParam[1],90.5,2,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights905_2,binnweights905_2 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio905_2)
    ratio905_3=fiteq(tmass,fitParam[0],fitParam[1],90.5,3,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights905_3,binnweights905_3 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio905_3)

    ratio91_2=fiteq(tmass,fitParam[0],fitParam[1],91,2,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights91_2,binnweights91_2 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio91_2)
    ratio91_1=fiteq(tmass,fitParam[0],fitParam[1],91,1,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights91_1,binnweights91_1 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio91_1)
    ratio91_3=fiteq(tmass,fitParam[0],fitParam[1],91,3,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights91_3,binnweights91_3 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio91_3)

    ratio915_1=fiteq(tmass,fitParam[0],fitParam[1],91.5,1,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights915_1,binnweights915_1 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio915_1)
    ratio915_2=fiteq(tmass,fitParam[0],fitParam[1],91.5,2,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights915_2,binnweights915_2 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio915_2)
    ratio915_3=fiteq(tmass,fitParam[0],fitParam[1],91.5,3,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights915_3,binnweights915_3 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio915_3)

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
    plt.title(f"Zreconstructed weight sim ({'real' if use_diagram==True else use_diagram })")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig(f"transient/Zreconstructed weight sim-combined_mess ({'real' if use_diagram==True else use_diagram}).pdf")
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
        return np.sum(((dataHist-interpolate_template)**2)/interpolate_template) #------------------------ch^2 function
    
#///////////////////////////////////////////////////////////////////////////////// calulcations for min mass and width
        # the Minuit doesnt seem to work without passign a fucntion through it 
    m = Minuit(interpolate_chi_func, 90.6, 2.4)#--------------------------------------minimise
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
    #c1,c2 variables to be changed in order to make graph convey right effect-ajust the tempaltes
    c1=mass_result+0.3 
    c2=width_result
    ratio_a=fiteq(tmass,fitParam[0],fitParam[1],c1,c2,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights_a,binnweights915_3 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio_a)

    sim_scale_factor_a= np.sum(dataHist) / np.sum(simmassHistweights_a)
    simHist_scaled_a = simmassHistweights_a * sim_scale_factor_a

    c1=mass_result-0.3 
    c2=width_result
    ratio_b=fiteq(tmass,fitParam[0],fitParam[1],c1,c2,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights_b,binnweights915_3 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio_b)

    sim_scale_factor_b= np.sum(dataHist) / np.sum(simmassHistweights_b)
    simHist_scaled_b = simmassHistweights_b * sim_scale_factor_b

    c1=mass_result
    c2=width_result+0.5
    ratio_c=fiteq(tmass,fitParam[0],fitParam[1],c1,c2,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights_c,binnweights915_3 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio_c)

    sim_scale_factor_c= np.sum(dataHist) / np.sum(simmassHistweights_c)
    simHist_scaled_c = simmassHistweights_c * sim_scale_factor_c

    c1=mass_result
    c2=width_result-0.5
    ratio_d=fiteq(tmass,fitParam[0],fitParam[1],c1,c2,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights_d,binnweights915_3 = np.histogram(conaeq(simdatam,c), bins=bin_space, weights=ratio_d)

    sim_scale_factor_d= np.sum(dataHist) / np.sum(simmassHistweights_d)
    simHist_scaled_d = simmassHistweights_d * sim_scale_factor_d




    best_interpolate_template=griddata(coords,templates,(mass_result,width_result),method='linear')#-----------------------best fit

    ax1.errorbar(centers, dataHist, yerr=dataerrors,label='Data ',color="black",fmt=".",markersize=2.5 )
    ax1.step(centers, best_interpolate_template, '-', linewidth=2,where='mid', label='best fit')
    ax1.set_ylabel("Counts")
    ax1.legend(loc='upper left',frameon=True, fontsize=8)
    ax1.set_ylim(bottom=0)
    # bottom
    ratio_a=simHist_scaled_a/best_interpolate_template
    
    ratio_b=simHist_scaled_b/best_interpolate_template

    ratio_fit=best_interpolate_template/best_interpolate_template

    ratio_c=simHist_scaled_c /best_interpolate_template

    ratio_d=simHist_scaled_d /best_interpolate_template

    ratio_data_91_2=dataHist/best_interpolate_template
    ratioerror=ratio_data_91_2*dataerrors/dataHist
    #ax2.scatter(centers,ratio_data_91,label='data ',color="black",s=10)
    ax2.errorbar(centers, ratio_data_91_2, yerr=ratioerror,label='Data ',color="black",fmt=".",markersize=2.5 )
    ax2.step(centers, ratio_a, '-', linewidth=1,where='mid', label='mass',color="red")
    ax2.step(centers, ratio_fit, '-', linewidth=1,where='mid',color="black")
    ax2.step(centers, ratio_b, '-', linewidth=1,where='mid', color="red")
    ax2.step(centers, ratio_c, '-', linewidth=1,where='mid', label='width ',color="blue")
    ax2.step(centers, ratio_d, '-', linewidth=1,where='mid', color="blue")
    ax2.set_ylabel("Ratio/best fit")
    ax2.set_xlabel("Mass / GeV")
    ax2.set_ylim(bottom=0.8)
    ax1.set_xlim(86, 96)
    ax2.legend(loc='upper left',frameon=True, fontsize=8)
    plt.savefig(f"transient/Z-stack_similtaneous ({'real' if use_diagram==True else use_diagram}).pdf")
                
    print(chi_values)

    # graph elipse of mass width
    plt.figure()
    eigenvalues , eigenvectors  = np.linalg.eigh(covariance_matrix) #data results---------------------
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
    ax.scatter(mass_result, width_result, color='red', label='Experimental Measurement')
    ax.set_xlabel("Mass") 
    ax.set_ylabel("Width") 
    ax.set_title("Mass Width error Ellipse") 
    
    ax.set_xlim(91.11, 91.22) 
    ax.set_ylim(2,2.56)

    theory_corelation=0.29342  #Theory------------------------------
    theory_mass=91.2047
    theory_massE=0.0088
    theory_width=2.49423
    theory_widthE=0.0023
    covariance_matrix = np.array([
    [theory_massE**2                             , theory_corelation*theory_massE*theory_widthE],  # 0.01 is the estiamte of the corelation coefficnet for theroretical cahnged to what a paper gave J blas electoweak fit
    [theory_corelation*theory_massE*theory_widthE, theory_widthE**2 ]                   ]) 
    eigenvalues,eigenvectors  = np.linalg.eigh(covariance_matrix)
    order = eigenvalues.argsort()[::-1] 
    eigenvalues = eigenvalues[order]
    eigenvectors=eigenvectors[:, order]  #making the bigegst eignecalue and coresponding vector be first in order
    major_axis_vector = eigenvectors[:,0]
    vx, vy = major_axis_vector[0], major_axis_vector[1]
    theta = np.degrees(np.arctan2(vy, vx))
    major_axis_length = 2 * np.sqrt(eigenvalues[0]) 
    minor_axis_length = 2 * np.sqrt(eigenvalues[1])
    ellipse2 = Ellipse( (theory_mass, theory_width), width=major_axis_length, height=minor_axis_length, angle=theta, edgecolor='blue', facecolor='none', linewidth=2)
    ax.add_patch(ellipse2)
    ax.scatter(theory_mass, theory_width, color='blue', label='standard model theory')
    #LEP elipse--------------------------------------
    theory_corelation= 0.038  
    theory_mass=91.1875
    theory_massE=0.0021
    theory_width=2.4952
    theory_widthE=0.0023
    covariance_matrix = np.array([
    [theory_massE**2                             , theory_corelation*theory_massE*theory_widthE],  # 0.01 is the estiamte of the corelation coefficnet for theroretical cahnged to what a paper gave J blas electoweak fit
    [theory_corelation*theory_massE*theory_widthE, theory_widthE**2 ]                   ]) 
    eigenvalues,eigenvectors  = np.linalg.eigh(covariance_matrix)
    order = eigenvalues.argsort()[::-1] 
    eigenvalues = eigenvalues[order]
    eigenvectors=eigenvectors[:, order]  #making the bigegst eignecalue and coresponding vector be first in order
    major_axis_vector = eigenvectors[:,0]
    vx, vy = major_axis_vector[0], major_axis_vector[1]
    theta = np.degrees(np.arctan2(vy, vx))
    major_axis_length = 2 * np.sqrt(eigenvalues[0]) 
    minor_axis_length = 2 * np.sqrt(eigenvalues[1])
    ellipse3 = Ellipse( (theory_mass, theory_width), width=major_axis_length, height=minor_axis_length, angle=theta, edgecolor='green', facecolor='none', linewidth=2)
    ax.add_patch(ellipse3)
    ax.scatter(theory_mass, theory_width, color='green', label='LEP')

    #---------------------------------------------other mesurments 
    #lhcb
    Z_mass_lhcb=91.1857
    lhcb_E=0.0083
    x_vals = np.linspace(*ax.get_xlim(), 500)
    ax.axvline(x=Z_mass_lhcb,label="LHCb prior measurment",color='black')
    ax.axvline(x=Z_mass_lhcb + lhcb_E, linestyle='--',color='black')
    ax.axvline(x=Z_mass_lhcb - lhcb_E, linestyle='--',color='black')
    #ax.fill_betweenx(y=ax.get_ylim(),x1=lhcb_E - lhcb_E,x2=lhcb_E + lhcb_E,alpha=0.3,color="grey")

 

    ax.legend(loc="upper left")
    plt.savefig(f"transient/Z-Error-Ellipse ({'real' if use_diagram==True else use_diagram}).pdf")
    
    #---------------------mass graph only
    fig, ax = plt.subplots()
    y=[5]
    x=[mass_result]
    x_er=[mass_error]
    ax.errorbar(x, y, xerr=x_er , fmt='o', color='red', capsize=3,label="Our result")
        
    y=[7] 
    x=[Z_mass_lhcb]
    x_er=[lhcb_E]
    ax.errorbar(x, y, xerr=x_er , fmt='o', color='black', capsize=3,label="Prior LHCb")
    sm_theory_mass=91.2047
    sm_theory_massE=0.0088
    y=[1]   
    x=[sm_theory_mass]
    x_er=[sm_theory_massE]
    ax.errorbar(x, y, xerr=x_er , fmt='o', color='blue', capsize=3,label="Standard model prediction")
    
    y=[3]   
    x=[theory_mass]
    x_er=[theory_massE]
    ax.errorbar(x, y, xerr=x_er , fmt='o', color='green', capsize=3,label="LEP")

    Z_mass_CDF=91.1943
    CDF_E=0.0138
    y=[9]   
    x=[Z_mass_CDF]
    x_er=[CDF_E]
    ax.errorbar(x, y, xerr=x_er , fmt='o', color='yellow', capsize=3,label="CDF")
    ax.legend(loc="upper left")
    ax.set_xlabel("Mass / GeV")
    ax.set_yticks([])
    plt.savefig(f"transient/Z-mass ({'real' if use_diagram==True else use_diagram}).pdf")
    
#selection cuts------------------------

muon_pt_p=datam["mup_pt"]
muon_pt_n=datam["mum_pt"]

muon_eta_p=datam["mup_eta"]
muon_eta_n=datam["mum_eta"]

pos_mask= (muon_eta_p < 4.4) & (muon_eta_n < 4.4) & (muon_eta_p > 2.2) & (muon_eta_n > 2.2) & (muon_pt_n > 20) & (muon_pt_p > 20)
pos_datam= {k: v[pos_mask] for k  , v in datam.items()}
seperate_Pt_psuedo_selection=True

sim_fits(tmass,simdatam,pos_datam,calibration_factor,True,50) #-All DATA goes through the cuts



if fname == f"mass-width_values_and_error.txt":
    fname="dont_use"
    # this is to test changing bin numbers------------------------------------------
    sim_fits(tmass,simdatam,pos_datam,calibration_factor,"30-bins",30)
    sim_fits(tmass,simdatam,pos_datam,calibration_factor,"80-bins",80)

# this is to test for magnetic monople 
#real means standard ie can use, watch out that adign the erros will chnage this so real is only trustworthy using general calibration and the run all script
# make sure that calibration ="true" is used and nothign else or if using other arguenmtns make sure jsut calibration=ture is used last
    magnet_pol=pos_datam["yearpol"]
    #print(f"magnet valeus are { magnet_pol}")
    pos_mask2= magnet_pol > 0
    neg_mask2= magnet_pol < 0
    pos_datam2= {k: v[pos_mask2] for k  , v in pos_datam.items()}
    neg_datam2= {k: v[neg_mask2] for k, v in pos_datam.items()}
    magnet="pos_dipole" 
    sim_fits(tmass,simdatam,pos_datam2,calibration_factor,magnet,50)
    magnet="neg_dipole" 
    sim_fits(tmass,simdatam,neg_datam2,calibration_factor,magnet,50)
# can add  more here qutie easily 
# this is splitting up agnles azimuthal of positive muon in detector

    muP_thi=pos_datam["mup_phi"]
    pos_mask2= muP_thi > 0
    neg_mask2= muP_thi < 0
    pos_datam2= {k: v[pos_mask2] for k  , v in pos_datam.items()}
    neg_datam2= {k: v[neg_mask2] for k, v in pos_datam.items()}

    angle="0__π" 
    sim_fits(tmass,simdatam,pos_datam2,calibration_factor,angle,50)
    angle="-0__π" 
    sim_fits(tmass,simdatam,neg_datam2,calibration_factor,angle,50)
# splitting by transeverse muon momentum 
    Dimuon_P=pos_datam["V_PT"]
    #print(np.mean(Dimuon_P))
    pos_mask2= Dimuon_P > 18.3
    neg_mask2=Dimuon_P <= 18.3
    pos_datam2= {k: v[pos_mask2] for k  , v in pos_datam.items()}
    neg_datam2= {k: v[neg_mask2] for k, v in pos_datam.items()}
    dimuon="Higher_Pt"
    sim_fits(tmass,simdatam,pos_datam2,calibration_factor,dimuon,50)
    dimuon="Lower_Pt"
    sim_fits(tmass,simdatam,neg_datam2,calibration_factor,dimuon,50)



# this splitting happens earlier anyway but is used as a check and to see what indilvidual selection cuts do
# splitting by individual transverse momenta-----more of a selection cut then a comaprison
    pos_mask= (muon_pt_n > 20) & (muon_pt_p > 20)
    pos_datam= {k: v[pos_mask] for k  , v in datam.items()}
    seperate_Pt_selection="seperate_Pt_selection"
    sim_fits(tmass,simdatam,pos_datam,calibration_factor,seperate_Pt_selection,50)
#splitting by individual pseudorapidity
    pos_mask= (muon_eta_p < 4.4) & (muon_eta_n < 4.4) & (muon_eta_p > 2.2) & (muon_eta_n > 2.2)
    pos_datam= {k: v[pos_mask] for k  , v in datam.items()}
    psuedo="psuedo"
    sim_fits(tmass,simdatam,pos_datam,calibration_factor,psuedo,50)
#splitting by both pseudorapidity and transverse momenta
    pos_mask= (muon_eta_p < 4.4) & (muon_eta_n < 4.4) & (muon_eta_p > 2.2) & (muon_eta_n > 2.2) & (muon_pt_n > 20) & (muon_pt_p > 20)
    pos_datam= {k: v[pos_mask] for k  , v in datam.items()}
    seperate_Pt_psuedo_selection="seperate_Pt_psuedo_selection"
    sim_fits(tmass,simdatam,pos_datam,calibration_factor,seperate_Pt_psuedo_selection,50)