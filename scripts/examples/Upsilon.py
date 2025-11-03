#!/usr/bin/env python

import uproot
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import argparse
from scipy.stats import crystalball

def ResFit(x,total,mean,sd,A,B):
    term = -0.5*((x-mean)**2 / sd**2)
    return A*np.exp(-B*x) + total / (np.sqrt(2*np.pi)*sd) * np.exp(term) #+D

    #from scipy.stats documentation crystallball(x,beta,m): x = (x-mean)/sd
def CrystalBallFit(x,beta,m,loc,scale,A,B):
    return crystalball.pdf(x,beta,m,loc=loc,scale=scale) + A*np.exp(-B*x)

def main():

    parser = argparse.ArgumentParser(description='some string')
    parser.add_argument('--Source',default="d",type=str)
    parser.add_argument('--Smearing',default="off",type=str)
    args = parser.parse_args()
    if (args.Source).lower() == "d" or (args.Source).lower() == "data":
        loc = "DATA"
    elif (args.Source).lower() == "s" or (args.Source).lower() =="sim":
        loc = "U1S"
    else:
        print("Invalid source given")
        return 1
    
    DATADIR="/storage/epp2/phshgg/Public/MPhysProject_2025_2026/tuples/0/"
    MUON_MASS = 0.1057
    with uproot.open(f"{DATADIR}/DecayTree__U1S__{loc}__d13600GeV_24c4.root:DecayTree") as t:
        #data = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX","mum_PY","mum_PZ"],library="np")
        data = t.arrays(["mup_P","mup_pt","mup_eta","mup_phi","mum_P","mum_pt","mum_eta","mum_phi"],library="np")

    data["mup_PX"] = data["mup_pt"]*np.cos(data["mup_phi"])
    data["mup_PY"] = data["mup_pt"]*np.sin(data["mup_phi"])
    data["mup_PZ"] = data["mup_pt"]*np.sinh(data["mup_eta"])
    data["mup_P"] = data["mup_pt"]*np.cosh(data["mup_eta"])
    data["mum_PX"] = data["mum_pt"]*np.cos(data["mum_phi"])
    data["mum_PY"] = data["mum_pt"]*np.sin(data["mum_phi"])
    data["mum_PZ"] = data["mum_pt"]*np.sinh(data["mum_eta"])
    data["mum_P"] = data["mum_pt"]*np.cosh(data["mum_eta"])

    mup_P,mum_P = np.array([data["mup_PX"],data["mup_PY"],data["mup_PZ"]]),np.array([data["mum_PX"],data["mum_PY"],data["mum_PZ"]])
    if (args.Smearing).lower() == "on" and loc == "U1S":
        mup_P,mum_P = mup_P*(1+np.random.normal(0,0.0055**2)),mum_P*(1+np.random.normal(0,0.55**2))
        smear = "_SmearingOn"
    elif loc =="U1S":
        smear = "_SmearingOff"
    else:
        smear = ""
    mup_E,mum_E = np.sqrt(data["mup_P"]**2+MUON_MASS**2),np.sqrt(data["mum_P"]**2+MUON_MASS**2)
    tot_E = mup_E + mum_E
    tot_PX = mup_P[0] + mum_P[0]
    tot_PY = mup_P[1] + mum_P[1]
    tot_PZ = mup_P[2] + mum_P[2]
    tot_P = np.sqrt(tot_PX**2+tot_PY**2+tot_PZ**2)
    mass = np.sqrt(tot_E**2 - tot_P**2)
    
    massHist,bins,_ = plt.hist(mass,bins=100,range=(9.200,9.750),histtype='step',label="Upsilon mass",density=True)
    
    binwidth = bins[1] - bins[0]
    binlist = [bins[0]+0.5*binwidth]
    for i in range(1,(len(bins)-1)):
        binlist.append(binlist[-1]+binwidth)
    bincenters = np.array(binlist)
    
    #Gaussian fit:
    #fitParam,_ = curve_fit(ResFit,bincenters,massHist,p0=[1.0, 9.46, 0.04, 0.5, 0.001],bounds=([0, 9.4, 0.018, 0, 0],[30.0, 9.5, 0.1, 10.0, 1e3]),maxfev=10000)
    #model = ResFit(bincenters,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    #plt.plot(bincenters,model,label=r'$ Ae^{-bx}+\mathrm{Gauss}(T,\bar{x},\sigma) $')
    fitParam = (curve_fit(CrystalBallFit,bincenters,massHist,p0=[1.75,2.5,9.45,1e-3,0.0,0.0],bounds=([0.0,0.0,9.4,1e-5,0.0,0.0],[10.0,5.0,9.5,2.0,30.0,10.0]),maxfev=10000))[0]
    print(fitParam)
    model = CrystalBallFit(bincenters,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4],fitParam[5])
    plt.plot(bincenters,model,label="Crystall Ball function\nWith Background")
    plt.legend()
    plt.xlabel("Mass / GeV")
    plt.ylabel("Frequency Density")
    plt.title(f"Reconstructed Upsilon {loc}, Smearing:{smear}")
    plt.savefig(f"Upsilon_mass_{loc}{smear}.pdf")
    plt.clf()

if __name__ == '__main__':
    main()