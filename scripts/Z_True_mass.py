#!/usr/bin/env python
import uproot
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
def fiteq(x,a,b,roe,scale):
    k=(2*np.sqrt(2)*roe**2*np.sqrt(roe**2+1))
    f= (((k/(((x**2-roe**2)**2+roe**2)*(np.pi*np.sqrt(roe**2+roe*np.sqrt(roe**2+1))))))  + a*np.exp(b*x))*scale  # i can just do this manually
    return f

DATADIRsim="/storage/epp2/phshgg/Public/MPhysProject_2025_2026/tuples/0/"

with uproot.open(f"{DATADIRsim}/MCDecayTree__Fiducial__Z__d13600GeV_24c4.root:MCDecayTree") as t:
    #print(t.keys())
    #momentasim = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    #print(momentasim)
    #DecayTree__Z__DATA__d13600GeV_24c4.root
    #DecayTree__Z__Z__d13600GeV_24c4.root
    for branch in t.keys():
        print(branch)
    massdatun=t.arrays(["true_bare_mass"],library="np")
    #for i in range(0,1000):
       # print(massdatun[i])

    masking=np.isfinite(massdatun["true_bare_mass"])
    #massdata=massdatun[masking]
    massdat=massdatun["true_bare_mass"][masking]

def plot(massdat):
    massHist,binn,_=plt.hist(massdat, bins=300, histtype="step",range=(65,115),label="Z mass",density=True)
    centers=0.5*(binn[1:]+binn[:-1])
    width=1
    M0=91.1876
    roe=M0/width
    fitParam,_ = curve_fit(fiteq,centers,massHist,p0=[5,-0.7,roe,0.4],bounds=([0.0,-1.0,10,0],[100.0,0,100,1]),maxfev=10000)
    print(fitParam)
    model = fiteq(centers,fitParam[0],fitParam[1],fitParam[2],fitParam[3])
    plt.plot(centers,model)
    #plt.plot(centers,fiteq(centers,1,1,30)*1e4,'g')
    plt.title("Z invariant true mass")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.savefig("Ztruegraph+model.pdf")
    plt.clf()
    print("gr")
plot(massdat)