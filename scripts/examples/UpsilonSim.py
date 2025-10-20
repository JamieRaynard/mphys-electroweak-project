#!/usr/bin/env python

import uproot
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

def ResFit(x,total,mean,sd,A,B,D):
    term = -0.5*((x-mean)**2 / sd**2)
    return A*np.exp(-B*x) + total / (np.sqrt(2*np.pi)*sd) * np.exp(term) +D

def main():
    DATADIR="/storage/epp2/phshgg/Public/ew_analyses/run3_tuples/11/"
    MUON_MASS = 105.7
    with uproot.open(f"{DATADIR}/d13600GeV_block7_Down_U1S_Sim10f_18112001_Turbo.root:U1S/DecayTree") as t:
        data = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX","mum_PY","mum_PZ"],library="np")
    
    mup_P,mum_P = np.array([data["mup_PX"],data["mup_PY"],data["mup_PZ"]]),np.array([data["mum_PX"],data["mum_PY"],data["mum_PZ"]])
    mup_E,mum_E = np.sqrt(mup_P[0]**2+mup_P[1]**2+mup_P[2]**2+MUON_MASS**2),np.sqrt(mum_P[0]**2+mum_P[1]**2+mum_P[2]**2+MUON_MASS**2)
    tot_E = mup_E + mum_E
    tot_PX = mup_P[0] + mum_P[0]
    tot_PY = mup_P[1] + mum_P[1]
    tot_PZ = mup_P[2] + mum_P[2]
    tot_P = np.sqrt(tot_PX**2+tot_PY**2+tot_PZ**2)
    mass = np.sqrt(tot_E**2 - tot_P**2)

    massHist,bins,_ = plt.hist(mass,bins=100,range=(9200,9750),histtype='step',label="Upsilon mass",density=True)
    binwidth = bins[1] - bins[0]
    binlist = [bins[0]+0.5*binwidth]
    for i in range(1,(len(bins)-1)):
        binlist.append(binlist[-1]+binwidth)
    bincenters = np.array(binlist)
    d_y = np.sqrt(massHist)
    fitParam,_ = curve_fit(ResFit,bincenters,massHist,p0=[0.5,9450,20,0.8,0.001,0],bounds=([0,9400,18,0,0,-0.03],[30,9500,100,10,1e3,0.03]),sigma=d_y,absolute_sigma=True)
    print(fitParam)
    model = ResFit(bincenters,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4],fitParam[5])
    plt.plot(bincenters,model)

    plt.legend()
    plt.xlabel("Mass / MeV")
    plt.ylabel("Frequency Density")
    plt.title("Upsilon invariant mass (Simulation)")
    plt.savefig("Upsilon_mass_sim.pdf")
    plt.clf()

if __name__ == '__main__':
    main()