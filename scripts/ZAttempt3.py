#!/usr/bin/env python
import uproot
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
DATAIR="/storage/epp2/phshgg/Public/MPhysProject_2025_2026/tuples/1/"
with uproot.open(f"{DATAIR}/DecayTree__Z__Z__d13600GeV_24c4.root:DecayTree") as t:

    #momentasim = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    simdatam=t.arrays(["mum_eta","mum_phi","mum_pt","mup_eta" ,"mup_phi" ,"mup_pt","true_boson_mass"],library="np")  #i ahve ztrue and i ahve z reconstructed simulation
    
    tmass=simdatam["true_boson_mass"]


with uproot.open(f"{DATAIR}/DecayTree__Z__DATA__d13600GeV_24c4.root:DecayTree") as tt:

    #momentasim = t.arrays(["mup_PX","mup_PY","mup_PZ","mum_PX" ,"mum_PY" ,"mum_PZ"],library="np")
    datam=tt.arrays(["mum_eta","mum_phi","mum_pt","mup_eta" ,"mup_phi" ,"mup_pt"],library="np")


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
    m=np.sqrt((EP+EN)**2-((mpx+mmx)**2+(mpy+mmy)**2+(mpz+mmz)**2))
    
    return m

def fiteq(x,a,b,m,w,scale):
    gamma=np.sqrt(m**2*(m**2+w**2))
    k=(2*np.sqrt(2)*m*gamma*w)/(np.pi*np.sqrt(m**2+gamma))
    f= ((k/((x**2-m**2)**2+(m**2)*w**2))*scale  + a*np.exp(b*x))  # i can just do this manually
    return f

def quad(da,a,b):
    f=a*da**2 +b
    return f

def TrueZfit_allplots(tmass,simdatam,datam):
    dataHist,databinn,_d=plt.hist(conaeq(datam), bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-data-reconstructed",density=True,linewidth=1) #the data recosntuction data
    remassHist,rebinn,_s=plt.hist(conaeq(simdatam), bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-sim reconstructed",density=True,linewidth=1) #for recosntructed simulation 
    trmassHist,binn,_t=plt.hist(tmass, bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-true",density=True,linewidth=1) #for true mass
    centers=0.5*(binn[1:]+binn[:-1])
    w=3
    m=91.1876
    fitParam,_tt = curve_fit(fiteq,centers,trmassHist,p0=[5,-0.7,m,w,0.4],bounds=([0.0,-1.0,60,0,0],[100.0,0,120,10,1]),maxfev=10000)# this uses the true mass
    #print(fitParam)
    model = fiteq(centers,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    plt.plot(centers,model)
    plt.title("Z-true fit+reconstructed")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig("Z-true fit +reconstructed.pdf")
    plt.clf()

def Zdata_with_fit(tmass,simdatam,datam):
    dataHist,databinn = np.histogram(conaeq(datam), bins=np.linspace(80.0,100.0,50),density=True)
    
    centersdat=0.5*(databinn[1:]+databinn[:-1])

    plt.scatter(centersdat,dataHist,label="Z-data")
    plt.title("Z-data with fit")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig("Z-data with fit.pdf")
    plt.clf()

def mikasuggestions(tmass,simdatam,datam):
    dataHist,databinn=np.histogram(conaeq(datam), bins=np.linspace(80.0,100.0,50))
    simHist,databinn=np.histogram(conaeq(simdatam), bins=np.linspace(80.0,100.0,50))
    chi90=np.sum(((dataHist- simHist)**2)/dataHist)
    print(chi90)

    centers=0.5*(databinn[1:]+databinn[:-1])
    sim_scale_factor = np.sum(dataHist) / np.sum(simHist) 
    simHist_scaled = simHist * sim_scale_factor
    dataHist,databinn,_d=plt.hist(conaeq(datam), bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-data-reconstructed",linewidth=1)
    plt.step(centers, simHist_scaled, where='mid', label="Z-sim scaled")
    plt.title("Z-data with scale sim")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency")
    plt.legend(loc='upper right')
    plt.savefig("Z-data-sim_scale.pdf")
    plt.clf()


def plot(tmass,simdatam,datam):
    dataHist,databinn,_d=plt.hist(conaeq(datam), bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-data-reconstructed",linewidth=1) #the data recosntuction data #for recosntructed simulation 
    trmassHist,binn,_t=plt.hist(tmass, bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-true",density=True,linewidth=1) #for true mass
    centers=0.5*(binn[1:]+binn[:-1])
    w=3
    m=91.1876
    fitParam,_tt = curve_fit(fiteq,centers,trmassHist,p0=[5,-0.7,m,w,0.4],bounds=([0.0,-1.0,60,0,0],[100.0,0,120,10,1]),maxfev=10000)# this uses the true mass
    print(fitParam)
    
    # ratio only has the t mass in it, the hsitogram suses this wieghtign witht the recosntructed in it
    #scaling the hisogram


    ratio90=fiteq(tmass,fitParam[0],fitParam[1],90,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights90,binnweights90 = np.histogram(conaeq(simdatam), bins=np.linspace(80.0,100.0,50), weights=ratio90)

    ratio91=fiteq(tmass,fitParam[0],fitParam[1],91,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights91,binnweights91 = np.histogram(conaeq(simdatam), bins=np.linspace(80.0,100.0,50), weights=ratio91)

    ratio92=fiteq(tmass,fitParam[0],fitParam[1],92,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights92,binnweights92 = np.histogram(conaeq(simdatam), bins=np.linspace(80.0,100.0,50), weights=ratio92)
    sim_scale_factor90= np.sum(dataHist) / np.sum(simmassHistweights90)
    simHist_scaled90 = simmassHistweights90 * sim_scale_factor90
    
    sim_scale_factor91= np.sum(dataHist) / np.sum(simmassHistweights91) 
    simHist_scaled91 = simmassHistweights91 * sim_scale_factor91

    sim_scale_factor92 = np.sum(dataHist) / np.sum(simmassHistweights92) 
    simHist_scaled92 = simmassHistweights92 * sim_scale_factor92

    plt.plot(centers, simHist_scaled90, '-', linewidth=2, label='scaled-weighted 90 simulation')
    plt.plot(centers, simHist_scaled91, '-', linewidth=2, label='scaled-weighted 91 simulation')
    plt.plot(centers, simHist_scaled92, '-', linewidth=2, label='scaled-weighted 92 simulation')

    plt.title("Z-reconstucted weighted sim")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig("weights graph.pdf")
    plt.clf()
    plt.figure() # now for the chi ^2 plot of from the last 3  // now this will ahve to be recosntucted daat from real life not sim

    chi90=np.sum(((dataHist-simHist_scaled90)**2)/dataHist)
    chi91=np.sum(((dataHist-simHist_scaled91)**2)/dataHist)
    chi92=np.sum(((dataHist-simHist_scaled92)**2)/dataHist)

    plt.title("Z-chi")
    y=np.array([chi90,chi91,chi92])
    print(y)
    x=np.array([90,91,92])
    plt.xlabel("target Zmass")
    plt.ylabel("Chi^2")

    plt.scatter(x,y)
    qfit = np.polyfit(x, y, deg=2)
    functionpfit = np.poly1d(qfit) 
    print(qfit)
    xx = np.linspace((89), (93), 300)
    plt.plot(xx,functionpfit(xx))
    plt.savefig("z-chi.pdf")

    min=-qfit[1]/(2*qfit[0])
    f=qfit[0]*min**2+qfit[1]*min+qfit[2]
    print("min ch^2 is",f)
    print("mz mass is", min)
    xerror_p=(-qfit[1]+np.sqrt(qfit[1]**2 -4*qfit[0]*(qfit[2]-1-f)))/(2*qfit[0])
    xerror_n=(-qfit[1]-np.sqrt(qfit[1]**2 -4*qfit[0]*(qfit[2]-1-f)))/(2*qfit[0])
    #print(xerror_p)
    #print(xerror_n)
    #print(min-xerror_p,min-xerror_n)
    print("the mass error is",min-xerror_n)
    plt.figure()


plot(tmass,simdatam,datam)
#TrueZfit_allplots(tmass,simdatam,datam)
#Zdata_with_fit(tmass,simdatam,datam)
#mikasuggestions(tmass,simdatam,datam)