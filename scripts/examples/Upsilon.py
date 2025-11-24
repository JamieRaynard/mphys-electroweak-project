#!/usr/bin/env python

import uproot
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit,minimize_scalar
import argparse
from scipy.stats import crystalball
from json import dump

#This is the Gaussian fit that I tried before the Crytall Ball, it is better at capturing the right tail
def ResFit(x,total,mean,sd,A,B):
    term = -0.5*((x-mean)**2 / sd**2)
    return A*np.exp(-B*x) + total / (np.sqrt(2*np.pi)*sd) * np.exp(term) #+D

#This Crystal Ball Fuction better encapsulates the QED radiative tail on the left side of the curve
#from scipy.stats documentation crystallball(x,beta,m): x = (x-mean)/sd
#Needed the loc and scale to act analogous to the mean and sd in a Gaussian
def CrystalBallFit(x,beta,m,loc,scale,A,B):
    return crystalball.pdf(x,beta,m,loc=loc,scale=scale) + A*np.exp(-B*x)

def GetBranches(loc):
    DATADIR="/storage/epp2/phshgg/Public/MPhysProject_2025_2026/tuples/0/"
    with uproot.open(f"{DATADIR}/DecayTree__U1S__{loc}__d13600GeV_24c4.root:DecayTree") as t:
        Rawdata = t.arrays(["mup_pt","mup_eta","mup_phi","mum_pt","mum_eta","mum_phi"],library="np")
    data = ConvertCoords(Rawdata)
    return data

def ConvertCoords(data):
    #File no longer contains px,py,pz but pt,eta,phi so needs to convert to reconstruct the mass
    #REMEMBER THESE ARE IN GEV!!
    data["mup_PX"] = data["mup_pt"]*np.cos(data["mup_phi"])
    data["mup_PY"] = data["mup_pt"]*np.sin(data["mup_phi"])
    data["mup_PZ"] = data["mup_pt"]*np.sinh(data["mup_eta"])
    data["mup_P"] = data["mup_pt"]*np.cosh(data["mup_eta"])
    data["mum_PX"] = data["mum_pt"]*np.cos(data["mum_phi"])
    data["mum_PY"] = data["mum_pt"]*np.sin(data["mum_phi"])
    data["mum_PZ"] = data["mum_pt"]*np.sinh(data["mum_eta"])
    data["mum_P"] = data["mum_pt"]*np.cosh(data["mum_eta"])
    return data

def Reconstruct(mup_P,mum_P,mup_E,mum_E):
    tot_E = mup_E + mum_E
    tot_PX = mup_P[0] + mum_P[0]
    tot_PY = mup_P[1] + mum_P[1]
    tot_PZ = mup_P[2] + mum_P[2]
    tot_P = np.sqrt(tot_PX**2+tot_PY**2+tot_PZ**2)
    mass_sq = np.maximum(tot_E**2 - tot_P**2,0)
    mass = np.sqrt(mass_sq)
    return mass

#loc and smear are just variables to determine the file name the graph will be saved under
def PlotHistogram(mass,loc,smear="",Output=None):
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

    #Crystal Ball fit:
    fitParam,cov = (curve_fit(CrystalBallFit,bincenters,massHist,p0=[1.75,2.5,9.45,1e-3,0.0,0.0],bounds=([0.0,0.0,9.4,1e-5,0.0,0.0],[10.0,5.0,9.5,2.0,30.0,10.0]),maxfev=10000))
    err = np.sqrt(np.diag(cov))
    print(fitParam,"\n",err)
    model = CrystalBallFit(bincenters,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4],fitParam[5])
    plt.plot(bincenters,model,label="Crystall Ball function\nWith Background")
    
    plt.legend()
    plt.xlabel("Mass / GeV")
    plt.ylabel("Frequency Density")
    plt.title(f"Reconstructed Upsilon {loc} {smear}")
    plt.savefig(f"transient/Upsilon_mass_{loc}{smear}.pdf")
    plt.clf()

    if Output == "alpha":
        return CalcAlpha(fitParam[2],err[2])
    elif Output == "width":
        return (fitParam[3],err[3])
    else:
        return 0

def Alpha(m,m_pdg):
    return (m/m_pdg) - 1

def CalcAlpha(m,m_err):
    M_PDG = 9.46040
    M_PDG_ERR = 0.00013
    alph = Alpha(m,M_PDG)
    err_m = Alpha(m+m_err,M_PDG) - alph
    err_mpdg = Alpha(m,M_PDG+M_PDG_ERR) - alph
    err_tot = np.sqrt(err_m**2+err_mpdg**2)
    return (alph,err_tot)

def c_ratio(alph_s,alph_d):
    return (1+alph_s)/(1+alph_d)

def CalcC(alpha_s,alpha_d):
    c = c_ratio(alpha_s[0],alpha_d[0])
    err_s = c_ratio(alpha_s[0]+alpha_s[1],alpha_d[0]) - c
    err_d = c_ratio(alpha_s[0],alpha_d[0]+alpha_d[1]) - c
    err_tot = np.sqrt(err_s**2+err_d**2)
    return (c,err_tot)

def width_chi2(sigma,width_data,width_data_err,mup_P_orig,mum_P_orig,mup_E,mum_E):
    sd_dif = np.sqrt(0.0455171**2 - 0.0400880351**2)
    factor = 1+(np.random.normal(0,sd_dif,size=len(mum_E)))*sigma
    mup_P = mup_P_orig*factor
    mum_P = mum_P_orig*factor
    mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    width_sim = PlotHistogram(mass,"U1S",Output="width")
    return (width_sim[0] - width_data)**2 / width_data_err**2

def CalcSmearFactor():
    MUON_MASS = 0.1057
    data = GetBranches("DATA")
    mup_P,mum_P = np.array([data["mup_PX"],data["mup_PY"],data["mup_PZ"]]),np.array([data["mum_PX"],data["mum_PY"],data["mum_PZ"]])
    mup_E,mum_E = np.sqrt(data["mup_P"]**2+MUON_MASS**2),np.sqrt(data["mum_P"]**2+MUON_MASS**2)
    mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    width = PlotHistogram(mass,"DATA",Output="width")
    best_sigma = minimize_scalar(width_chi2,args=(width[0],width[1],mup_P,mum_P,mup_E,mum_E),bounds=(0.0,0.1),method="bounded")
    return best_sigma.x

def main():
    MUON_MASS = 0.1057
    #This allows for me to pass arguments in to the terminal to change important paramters without changing the code
    parser = argparse.ArgumentParser(description='some string')
    #Run with --Source="s" to switch to simulation
    parser.add_argument('--Source',default="s",type=str)
    #Run with --Smearing="on" to smear the momentum; will only work for simulation
    parser.add_argument('--Smearing',default="off",type=str)
    parser.add_argument('--Calibration',default="off",type=str)
    parser.add_argument('--FullOutput',default="FALSE",type=str)
    args = parser.parse_args()

    if (args.Source).lower() == "d" or (args.Source).lower() == "data":
        loc = "DATA"
    elif (args.Source).lower() == "s" or (args.Source).lower() =="sim":
        loc = "U1S"
    else:
        print("Invalid source given")
        return 1
    if (args.FullOutput).lower() == "true" or (args.Calibration).lower() == "on":
        loc = "U1S"

    data = GetBranches(loc)
    output = {}
    
    mup_P,mum_P = np.array([data["mup_PX"],data["mup_PY"],data["mup_PZ"]]),np.array([data["mum_PX"],data["mum_PY"],data["mum_PZ"]])
    mup_E,mum_E = np.sqrt(data["mup_P"]**2+MUON_MASS**2),np.sqrt(data["mum_P"]**2+MUON_MASS**2)
    
    if ((args.Smearing).lower() == "on" and loc == "U1S") or (args.FullOutput).lower() == "true":
        #This applies a Gaussian smearing to the simulated momenta to try to make them more like the real data
        sigma = CalcSmearFactor()
        print(f'WOOO got a scale variable: {sigma}')
        sd_dif = np.sqrt(0.0455171**2 - 0.0400880351**2)
        factor = 1+(np.random.normal(0,sd_dif,size=len(mum_E)))*sigma
        mup_P *= factor
        mum_P *= factor
        #This is just convinient for the file name
        smear = "_SmearingOn"
        output["Smear_factor"] = sigma
    elif loc =="U1S":
        smear = "_SmearingOff"
    else:
        smear = ""
    
    mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)

    if (args.FullOutput).lower() == "true" or (args.Calibration).lower() == "on":
        #This is a tuple of format (value,uncertainty)
        alpha_s = PlotHistogram(mass,loc,Output="alpha")
        data = GetBranches("DATA")
        mass = Reconstruct(data)
        alpha_d = PlotHistogram(mass,"DATA",Output="alpha",smear=smear)
        c = CalcC(alpha_s,alpha_d)
        output["C_ratio"] =  c
    else:
        PlotHistogram(mass,loc,smear=smear)

    with open("Calibration_output.json","w") as OutputFile:
            dump(output,OutputFile,indent=2)

    return 0

        

if __name__ == '__main__':
    main()