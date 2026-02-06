#!/usr/bin/env python

import uproot
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit,minimize_scalar
import argparse
from scipy.stats import crystalball
from json import dump,load

#This is the Gaussian fit that I tried before the Crytall Ball, it is better at capturing the right tail
def ResFit(x,total,mean,sd,A,B):
    term = -0.5*((x-mean)**2 / sd**2)
    return A*np.exp(-B*x) + total / (np.sqrt(2*np.pi)*sd) * np.exp(term) #+D

#This Crystal Ball Fuction better encapsulates the QED radiative tail on the left side of the curve
#from scipy.stats documentation crystallball(x,beta,m): x = (x-mean)/sd
#Needed the loc and scale to act analogous to the mean and sd in a Gaussian
def CrystalBallFit(x,beta,m,loc,scale,N,A,B,binwidth):
    return (N*crystalball.pdf(x,beta,m,loc=loc,scale=scale) + A*np.exp(-B*x))*binwidth

def GetBranches(loc,calibration_factor=1):
    MUON_MASS = 0.1057
    DATADIR="/storage/epp2/phshgg/Public/MPhysProject_2025_2026/tuples/0/"
    with uproot.open(f"{DATADIR}/DecayTree__U1S__{loc}__d13600GeV_24c4.root:DecayTree") as t:
        Rawdata = t.arrays(["mup_pt","mup_eta","mup_phi","mum_pt","mum_eta","mum_phi"],library="np")
    data = ConvertCoords(Rawdata) #root file is in eta,phi but we want px,py,pz
    mup_P,mum_P = np.array([data["mup_PX"],data["mup_PY"],data["mup_PZ"]]),np.array([data["mum_PX"],data["mum_PY"],data["mum_PZ"]])
    mup_P = mup_P*calibration_factor
    mum_P = mum_P*calibration_factor
    mup_P_mag = np.sqrt(mup_P[0]**2+mup_P[1]**2+mup_P[2]**2)
    mum_P_mag = np.sqrt(mum_P[0]**2+mum_P[1]**2+mum_P[2]**2)
    mup_E,mum_E = np.sqrt(mup_P_mag**2+MUON_MASS**2),np.sqrt(mum_P_mag**2+MUON_MASS**2)
    #For now I am only using these 4 values
    return mup_P,mum_P,mup_E,mum_E

def ConvertCoords(data):
    #File no longer contains px,py,pz but pt,eta,phi so needs to convert to reconstruct the mass
    #REMEMBER THESE ARE IN GEV!!
    data["mup_PX"] = data["mup_pt"]*np.cos(data["mup_phi"])
    data["mup_PY"] = data["mup_pt"]*np.sin(data["mup_phi"])
    data["mup_PZ"] = data["mup_pt"]*np.sinh(data["mup_eta"])
    data["mum_PX"] = data["mum_pt"]*np.cos(data["mum_phi"])
    data["mum_PY"] = data["mum_pt"]*np.sin(data["mum_phi"])
    data["mum_PZ"] = data["mum_pt"]*np.sinh(data["mum_eta"])
    return data

def Reconstruct(mup_P,mum_P,mup_E,mum_E):
    tot_E = mup_E + mum_E
    tot_PX = mup_P[0] + mum_P[0]
    tot_PY = mup_P[1] + mum_P[1]
    tot_PZ = mup_P[2] + mum_P[2]
    tot_P = np.sqrt(tot_PX**2+tot_PY**2+tot_PZ**2)
    #This just stops any issues of having E^2 < P^2, which shouldn't happen now I fixed calibration anyway
    mass_sq = np.maximum(tot_E**2 - tot_P**2,0)
    mass = np.sqrt(mass_sq)
    return mass

#loc and smear are just variables to determine the file name the graph will be saved under
def PlotHistogram(mass,filename,Output=None):
    #xmassHist, xbins, x = plt.hist(mass,bins=100,range=(9.2,9.75),histtype='step')
    massHist,bins = np.histogram(mass,bins=100,range=(9.200,9.750))
    binwidth = bins[1] - bins[0]
    binlist = [bins[0]+0.5*binwidth]
    for i in range(1,(len(bins)-1)):
        binlist.append(binlist[-1]+binwidth)
    bincenters = np.array(binlist)
    N_tot = np.sum(massHist)
    plt.scatter(bincenters,massHist,s=5 ,c='black')
    plt.errorbar(bincenters,massHist,yerr=np.sqrt(massHist),fmt='none',ecolor='black',elinewidth=1,capsize=2)
    
    #Crystal Ball fit:
    p0 = [1.19698532e+00,1.33208227e+00,9.45914418e+00,4.94449266e-02,N_tot,0.5*N_tot,9.9e-5]
    bounds = ([0.5, 1.0, 9.40, 0.005, 0.0, 0.0, 0.0], [5.0, 10.0, 9.50, 0.10,  N_tot, 10*N_tot, 1.0])
    fitParam, cov = curve_fit(lambda x, beta, m, loc, scale, Ns, A, B: CrystalBallFit(x, beta, m, loc, scale, Ns, A, B, binwidth), bincenters, massHist, p0=p0, bounds=bounds, maxfev=100000)
    err = np.sqrt(np.diag(cov))
    print(fitParam,"\n",err)
    model = CrystalBallFit(bincenters,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4],fitParam[5],fitParam[6],binwidth)
    plt.plot(bincenters,model,label="Crystall Ball function\nWith Background")
    #A more accurate fit could be a double tailed crystal ball
    
    #print(f'Saving plot to transient/Upsilon_mass_{filename}.pdf')

    plt.legend()
    plt.xlabel("Mass / GeV")
    plt.ylabel("Counts")
    plt.ylim(bottom=0)
    #plt.title(f"Reconstructed Upsilon {filename}")
    plt.savefig(f"transient/Upsilon_mass_{filename}.png")
    plt.clf()

    if Output:
        #At the moment these are the only values I use but can add others (e.g. for A and B for background) easily here
        outputvalues = {
            "mass": (fitParam[2],err[2]),
            "width": (fitParam[3],err[3]),
            "A": (fitParam[5],err[5]),
            "B": (fitParam[6],err[6])
        }
        return (outputvalues)
    else:
        return 0

def CompareHistograms(data_mass,unscaled_sim_mass,scaled_sim_mass):
    data_massHist, bins = np.histogram(data_mass, bins=100, range=(9.25,9.75))
    binwidth = bins[1] - bins[0]
    binlist = [bins[0]+0.5*binwidth]
    for i in range(1,(len(bins)-1)):
        binlist.append(binlist[-1]+binwidth)
    bincenters = np.array(binlist)

    #(Probably too) simplistic model for the background of just a flat uniform dist.
    #background = np.min(data_massHist)

    fitParam = PlotHistogram(data_mass,'DATA_fit',Output=True)
    background = fitParam["A"][0]*np.exp(-1*float(fitParam["B"][0])*bincenters)*binwidth

    unscaled_sim_massHist,bins = (np.histogram(unscaled_sim_mass, bins=100, range=(9.25,9.75)))
    scaled_sim_massHist,bins = (np.histogram(scaled_sim_mass, bins = 100, range = (9.25,9.75)))

    #Original
    # unscaled_sim_massHist = unscaled_sim_massHist + background*(np.sum(unscaled_sim_massHist)/np.sum(data_massHist))
    # scaled_sim_massHist = scaled_sim_massHist + background*(np.sum(scaled_sim_massHist)/np.sum(data_massHist))
    # unscaled_sim_massHist = unscaled_sim_massHist * (np.sum(data_massHist)/np.sum(unscaled_sim_massHist))
    # scaled_sim_massHist = scaled_sim_massHist * (np.sum(data_massHist)/np.sum(scaled_sim_massHist))

    #Idea:
    data_massHist_noBG = data_massHist - background
    data_massHist_noBG[data_massHist_noBG < 0.0] = 0.0
    unscaled_sim_massHist = unscaled_sim_massHist * (np.sum(data_massHist_noBG)/np.sum(unscaled_sim_massHist))
    scaled_sim_massHist = scaled_sim_massHist * (np.sum(data_massHist_noBG)/np.sum(scaled_sim_massHist))
    

    plt.bar(bincenters, background, width=binwidth, label="Background", color="lightgray", align="center")
    plt.step(bincenters, unscaled_sim_massHist+background,where="mid",label="Sim without smearing",color="blue",zorder=2)
    plt.step(bincenters, scaled_sim_massHist+background, where="mid", label="Sim with smearing",color="orange",zorder=2)
    plt.scatter(bincenters, data_massHist, label = "Data", s=2 ,c='black',zorder=3)
    plt.errorbar(bincenters, data_massHist, yerr=np.sqrt(data_massHist),fmt='none')

    # unscaled_sim_massHist = unscaled_sim_massHist + background
    # scaled_sim_massHist = scaled_sim_massHist + background
    # plt.step(bincenters,unscaled_sim_massHist,label="Sim without smearing")
    # plt.step(bincenters,scaled_sim_massHist, label = "sim with smearing")

    plt.legend()
    plt.xlabel("Mass / GeV")
    plt.ylabel("Counts")
    plt.ylim(bottom=0)
    #plt.title(r"Comparing the effect of momentum smearing")
    plt.savefig(f"transient/Upsilon_mass_comparisson.png")
    plt.clf()

    return 0

def Comparing():
    try:
        with open("Calibration_output.json",) as InputFile:
            Calibration = load(InputFile)
        alpha = 1 - Calibration["C_ratio"][0]
        Smear_factor = Calibration["Smear_factor"][0]
    except FileNotFoundError:
        print('Please run the script with --FullOutput="TRUE" first to get calibration information')
        return 1
    except KeyError:
        print('Please run the script with --FullOutput="TRUE" first to get all required calibration information')
        return 1
    
    mup_P,mum_P,mup_E,mum_E = GetBranches("DATA")
    data_mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    mup_P,mum_P,mup_E,mum_E = GetBranches("U1S",calibration_factor=1+alpha)
    unscaled_sim_mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    rng = np.random.default_rng(seed=10)
    Norm_rand = rng.normal(0,1,size=len(mum_E))
    calibration_factor = 1+alpha+Norm_rand*Smear_factor
    mup_P,mum_P,mup_E,mum_E = GetBranches("U1S",calibration_factor=calibration_factor)
    scaled_sim_mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    CompareHistograms(data_mass,unscaled_sim_mass,scaled_sim_mass)
    return 0

def CalcScaling():
    mup_P,mum_P,mup_E,mum_E = GetBranches("U1S")
    mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    sim_results = PlotHistogram(mass,"U1S",Output=True)
    alpha_s = CalcAlpha(sim_results["mass"][0],sim_results["mass"][1])
    mup_P,mum_P,mup_E,mum_E = GetBranches("DATA")
    mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    data_results = PlotHistogram(mass,"DATA",Output=True)
    alpha_d = CalcAlpha(data_results["mass"][0],data_results["mass"][1])
    return CalcC(alpha_s,alpha_d)


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
    return (1+alph_d)/(1+alph_s)
    #I am not 100% sure which way around these should be

def CalcC(alpha_s,alpha_d):
    c = c_ratio(alpha_s[0],alpha_d[0])
    err_s = c_ratio(alpha_s[0]+alpha_s[1],alpha_d[0]) - c
    err_d = c_ratio(alpha_s[0],alpha_d[0]+alpha_d[1]) - c
    err_tot = np.sqrt(err_s**2+err_d**2)
    return (c,err_tot)

def width_chi2(sigma,width_data,width_data_err,mup_P_orig,mum_P_orig,mup_E,mum_E):
    #Outdated: need to recalculate the mup_E and mum_E with smearing
    rng = np.random.default_rng(seed=10)
    Norm_rand = rng.normal(0,1,size=len(mum_E))
    factor = 1+Norm_rand*sigma
    mup_P = mup_P_orig*factor
    mum_P = mum_P_orig*factor
    mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    width_sim = PlotHistogram(mass,"U1S",Output=True)["width"]
    return (width_sim[0] - width_data)**2 / width_data_err**2

def CalcSmearFactorByMinimise():
    #Currently unused
    mup_P,mum_P,mup_E,mum_E = GetBranches("DATA")
    mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    data_width = PlotHistogram(mass,"DATA",Output=True)["width"]
    mup_P,mum_P,mup_E,mum_E = GetBranches("U1S")
    best_sigma = minimize_scalar(width_chi2,args=(data_width[0],data_width[1],mup_P,mum_P,mup_E,mum_E),bounds=(0.0,0.1),method="bounded")
    return best_sigma.x

def SmearFactor(sim_width,sim_mass,data_width,data_mass):
    return np.sqrt((data_width/data_mass)**2-(sim_width/sim_mass)**2)

def CalcSmearFactor():
    mup_P,mum_P,mup_E,mum_E = GetBranches("U1S")
    sim_mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    sim_results = PlotHistogram(sim_mass,"U1S",Output=True)

    mup_P,mum_P,mup_E,mum_E = GetBranches("DATA")
    data_mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    data_results = PlotHistogram(data_mass,"DATA",Output=True)
    
    print(f'Data width: {data_results["width"][0]} ± {data_results["width"][1]} \nUnsmeared sim width: {sim_results["width"][0]} ± {sim_results["width"][1]}')
    sigma = SmearFactor(sim_results["width"][0],sim_results["mass"][0],data_results["width"][0],data_results["mass"][0])
    err_due_sim_width = SmearFactor(sim_results["width"][0]+sim_results["width"][1],sim_results["mass"][0],data_results["width"][0],data_results["mass"][0]) - sigma
    err_due_sim_mass = SmearFactor(sim_results["width"][0],sim_results["mass"][0]+sim_results["mass"][1],data_results["width"][0],data_results["mass"][0]) - sigma
    err_due_data_width = SmearFactor(sim_results["width"][0],sim_results["mass"][0],data_results["width"][0]+data_results["width"][1],data_results["mass"][0]) - sigma
    err_due_data_mass = SmearFactor(sim_results["width"][0],sim_results["mass"][0],data_results["width"][0],data_results["mass"][0]+data_results["mass"][1]) - sigma
    err_sigma = np.sqrt(err_due_sim_width**2+err_due_sim_mass**2+err_due_data_width**2+err_due_data_mass**2)
    return (sigma,err_sigma)

def main():
    #This allows for me to pass arguments in to the terminal to change important paramters without changing the code
    parser = argparse.ArgumentParser(description='some string')
    #Run with --Source="s" to switch to simulation
    parser.add_argument('--Source',default="s",type=str)
    #Run with --Smearing="on" to smear the momentum; will only work for simulation
    parser.add_argument('--Smearing',default="off",type=str)
    parser.add_argument('--Calibration',default="off",type=str)
    parser.add_argument('--FullOutput',default="FALSE",type=str)
    parser.add_argument('--Compare',default="FALSE",type=str)
    args = parser.parse_args()

    #If you run compare, that overwrites all other arguments
    if (args.Compare).lower() == "true":
        success = Comparing()
        return success

    if (args.Source).lower() == "d" or (args.Source).lower() == "data":
        loc = "DATA"
    elif (args.Source).lower() == "s" or (args.Source).lower() =="sim":
        loc = "U1S"
    else:
        print("Invalid source given")
        return 1
    if (args.FullOutput).lower() == "true":
        loc = "U1S"

    output = {}
    filename=loc
    factor=1

    #Smearing
    if (((args.Smearing).lower() == "on"  or (args.Smearing).lower() == "true") and loc == "U1S") or (args.FullOutput).lower() == "true":
        #This applies a Gaussian smearing to the simulated momenta to try to make them more like the real data
        sigma,sigma_err = CalcSmearFactor()
        print(f'WOOO got a smearing variable: {sigma} ± {sigma_err}')
        rng = np.random.default_rng(seed=10)
        mup_P,mum_P,mup_E,mum_E = GetBranches(loc,calibration_factor=factor) #temp solution: needs improving
        Norm_rand = rng.normal(0,1,size=len(mum_E))
        factor += (Norm_rand*sigma)
        output["Smear_factor"] = (sigma,sigma_err)
        filename = filename+"_Smeared"

    #scaling
    if  ((args.Calibration).lower()) == "on" or ((args.Calibration).lower() == "true") or ((args.FullOutput).lower() == "true"):
        c,c_err = CalcScaling()
        print(f'WOOO got a scaling variable: {c} ± {c_err}')
        factor += c
        output["C_ratio"] = (c,c_err)
        filename= filename+"_Scaled"
    
    #I don't think this is needed
    # mup_P,mum_P,mup_E,mum_E = GetBranches(loc,calibration_factor=factor)
    # mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    # PlotHistogram(mass,filename)
    # if (args.FullOutput).lower() == "true":
    #     mup_P,mum_P,mup_E,mum_E = GetBranches("DATA",calibration_factor=factor)
    #     mass = Reconstruct(mup_P,mum_P,mup_E,mum_E)
    #     PlotHistogram(mass,"DATA_Smeared_Scaled")

    with open("Calibration_output.json","w") as OutputFile:
            dump(output,OutputFile,indent=2)

    return 0

        

if __name__ == '__main__':
    main()