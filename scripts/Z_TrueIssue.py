#!/usr/bin/env python
import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import rel_breitwigner
from scipy.optimize import curve_fit
def fiteq(x,a,b,roe,w):

    f=a*np.exp(b*x) +  rel_breitwigner.pdf(x,roe,scale=w)     # i can just do this manually
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
    massHist,binn,_=plt.hist(massdat, bins=300, histtype="step",range=(58,122),label="Z mass",density=True)
    centers=0.5*(binn[1:]+binn[:-1])
    w=5
    M0=90
    roe=M0/w
    fitParam,_ = curve_fit(fiteq,centers,massHist,p0=[0.5,2,roe,w],bounds=([0,-np.inf,2,1],[100,np.inf,80,20]))
    print(fitParam)
    model = fiteq(centers,fitParam[0],fitParam[1],fitParam[2],fitParam[3])
    plt.plot(centers,model)
    plt.title("Z invariant true mass")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.savefig("Ztruegraph+model.pdf")
    plt.clf()
    print("gr")
plot(massdat)