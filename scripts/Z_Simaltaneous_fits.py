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
#where is the main function? :(

#This allows for me to pass arguments in to the terminal to change important paramters without changing the code
parser = argparse.ArgumentParser(description='some string')
#Run with --Callibration="TRUE" to use callibration
parser.add_argument('--Calibration',default="FALSE",type=str)
args = parser.parse_args()
if (args.Calibration).lower() == "true":
    try:
        with open("C_ratio.json",) as InputFile:
            C_rat,C_err =  load(InputFile)["C_ratio"]
        print(f'Woo we imported {C_rat} ± {C_err}')
    except FileNotFoundError:
        print("Please run the script Upsilon.py first to get a calibration value")
        #This is where I would return to stop the script if we were inside a main func to avoid error
        C_rat,C_err = 1,0
else:
    C_rat = 1
    C_err = 0


DATAIR="/storage/epp2/phshgg/Public/MPhysProject_2025_2026/tuples/1/"
with uproot.open(f"{DATAIR}/DecayTree__Z__Z__d13600GeV_24c4.root:DecayTree") as t:

    #momentasim = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    simdatam=t.arrays(["mum_eta","mum_phi","mum_pt","mup_eta" ,"mup_phi" ,"mup_pt","true_boson_mass"],library="np")  #i ahve ztrue and i ahve z reconstructed simulation
    
    tmass=simdatam["true_boson_mass"]


with uproot.open(f"{DATAIR}/DecayTree__Z__DATA__d13600GeV_24c4.root:DecayTree") as tt:

    #momentasim = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    datam=tt.arrays(["mum_eta","mum_phi","mum_pt","mup_eta" ,"mup_phi" ,"mup_pt"],library="np")


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

def sim_fits(tmass,simdatam,datam,calibrate,err):
    c=calibrate #error in calibration
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
    plt.title(f"Zreconstructed weight sim ({'real' if err==True else 'dont-use'})")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig(f"transient/Zreconstructed weight sim-combined_mess ({'real' if err==False else 'dont-use'}).pdf")
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
    z=np.array([chi905_1,chi905_2,chi905_3,chi91_1,chi91_2,chi91_3,chi915_1,chi915_2,chi915_3])
    #xi = np.linspace(x.min()-0.1, x.max()+0.1, 100)
   # yi = np.linspace(y.min()-0.1, y.max()+0.1, 100)
    #X, Y = np.meshgrid(xi, yi)
    #Z = griddata((x, y), z, (X, Y), method='cubic')
    t = np.arange(len(x))  
    cs_x = CubicSpline(t, x)
    cs_y = CubicSpline(t, y)
    cs_z = CubicSpline(t, z)
    t_fine = np.linspace(t.min(), t.max(), 200)
    x_fine = cs_x(t_fine)
    y_fine = cs_y(t_fine)
    z_fine = cs_z(t_fine)
    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')
    #ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)
    ax.plot(x_fine, y_fine, z_fine, color='blue', label='Cubic Spline line')
    ax.scatter(x, y, z,color="red",label="data")
    ax.set_xlabel("target width")
    ax.set_ylabel("target mass")
    ax.set_zlabel("chi²")
    ax.legend()
    plt.savefig(f"transient/Z width+mass CHI({'real' if err==False else 'dont-use'}).pdf")

    xi = np.linspace(x.min(), x.max(), 200)
    yi = np.linspace(y.min(), y.max(), 200)
    Xi, Yi = np.meshgrid(xi, yi)


    Zi = griddata((x, y),z,(Xi, Yi),method='linear')

    cmap =ListedColormap(["blue","green"])
    plt.figure()
    plt.pcolormesh(Xi, Yi, Zi,edgecolors='face', shading='auto',cmap='Blues_r')
    plt.xlabel("tempalte width")
    plt.ylabel("template mass")
    plt.colorbar(label="chi^2")
    plt.show()
    plt.savefig(f"transient/Z heatmap({'real' if err==False else 'dont-use'}).pdf")

sim_fits(tmass,simdatam,datam,C_rat+C_err,False)