#!/usr/bin/env python
import uproot
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
from json import load,dump
import argparse

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