#!/usr/bin/env python
import uproot
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
def fiteq(x,a,b,m,w,scale):
    gamma=np.sqrt(m**2*(m**2+w**2))
    k=(2*np.sqrt(2)*m*gamma*w)/(np.pi*np.sqrt(m**2+gamma))
    f= ((k/((x**2-m**2)**2+(m**2)*w**2))*scale  + a*np.exp(b*x))  # i can just do this manually
    return f

DATADIRt="/storage/epp2/phshgg/Public/MPhysProject_2025_2026/tuples/0/"

with uproot.open(f"{DATADIRt}/MCDecayTree__Fiducial__Z__d13600GeV_24c4.root:MCDecayTree") as t:
    #print(t.keys())
    #momentasim = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    #print(momentasim)
    #DecayTree__Z__DATA__d13600GeV_24c4.root
    #DecayTree__Z__Z__d13600GeV_24c4.root
    #for branch in t.keys():
        #print(branch)
    massdatun=t.arrays(["true_bare_mass"],library="np")
    #for i in range(0,1000):
       # print(massdatun[i])

    masking=np.isfinite(massdatun["true_bare_mass"])
    #massdata=massdatun[masking]
    massdat=massdatun["true_bare_mass"][masking]

DATADIRrecon="/storage/epp2/phshgg/Public/MPhysProject_2025_2026/tuples/0/"

with uproot.open(f"{DATADIRrecon}/DecayTree__Z__Z__d13600GeV_24c4.root:DecayTree") as tt:
    #print(t.keys())
    #momentasim = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    #print(momentasim)
    #DecayTree__Z__DATA__d13600GeV_24c4.root
    #DecayTree__Z__Z__d13600GeV_24c4.root
    #for branch in tt.keys():
       # print(branch)
    reconmass=tt.arrays(["mum_eta","mum_phi","mum_pt","mup_eta" ,"mup_phi" ,"mup_pt"],library="np")
    #for i in range(0,1000):
       # print(massdatun[i])

    #masking=np.isfinite(massdatun["true_bare_mass"])
    #massdat=massdatun["true_bare_mass"][masking]
def conaeq(reconmass):
    mpe=reconmass["mup_eta"]
    mppt=reconmass["mup_pt"] 
    mpphi=reconmass["mup_phi"] 
    mme=reconmass["mum_eta"] 
    mmpt=reconmass["mum_pt"] 
    mmphi=reconmass["mum_phi"]
    mpx=mppt*np.cos(mpphi)
    mmx=mmpt*np.cos(mmphi)
    mpy=mppt*np.sin(mpphi)
    mmy=mmpt*np.sin(mmphi)
    mpz=mppt*np.sinh(mpe)
    mmz=mmpt*np.sinh(mme)

    m=[]
    Ep=[]
    Em=[]
    M_mass=105.658*10**(-3)
    EP=np.sqrt(M_mass**2 +mpx**2 +mpy**2+mpz**2)   
    EN=np.sqrt(M_mass**2 +mmx**2 +mmy**2+mmz**2)
    v=np.sqrt((EP+EN)**2-((mpx+mmx)**2+(mpy+mmy)**2+(mpz+mmz)**2))
    m.append(v)
    return m

recmass=conaeq(reconmass)
def plot(massdat,recmass):
    massHist,binn,_=plt.hist(massdat, bins=300, histtype="step",range=(65,115),label="Z mass true",density=True,linewidth=1)
    massHistret,binnret,_ret=plt.hist(recmass, bins=300, histtype="step",range=(65,115),label="Z mass reconstructed",density=True,linewidth=1)
    centers=0.5*(binn[1:]+binn[:-1])
    w=3
    m=91.1876
    fitParam,_ = curve_fit(fiteq,centers,massHist,p0=[5,-0.7,m,w,0.4],bounds=([0.0,-1.0,60,0,0],[100.0,0,120,10,1]),maxfev=10000)
    print(fitParam)
    model = fiteq(centers,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    plt.plot(centers,model)
    #plt.plot(centers,fiteq(centers,1,1,30)*1e4,'g')
    plt.title("Z invariant true mass")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig("Ztrue-ret-fit.pdf")
    plt.clf()
    print("gr")
plot(massdat,recmass)