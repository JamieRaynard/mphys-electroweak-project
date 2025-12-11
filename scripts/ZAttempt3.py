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


def plotmass(tmass,simdatam,datam,calibrate,err):
    c=calibrate #error in calibration
    dataHist,databinn,_d=plt.hist(conaeq(datam,1), bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-data-reconstructed",linewidth=1) #the data recosntuction data #for recosntructed simulation 
    trmassHist,binn,_t=plt.hist(tmass, bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-true",density=True,linewidth=1) #for true mass
    centers=0.5*(binn[1:]+binn[:-1])
    w=3
    m=91.1876
    fitParam,_tt = curve_fit(fiteq,centers,trmassHist,p0=[5,-0.7,m,w,0.4],bounds=([0.0,-1.0,60,0,0],[100.0,0,120,10,1]),maxfev=10000)# this uses the true mass
    #print(fitParam)
    
    # ratio only has the t mass in it, the hsitogram suses this wieghtign witht the recosntructed in it
    #scaling the hisogram

    ratio905=fiteq(tmass,fitParam[0],fitParam[1],90.5,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights905,binnweights905 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio905)

    ratio91=fiteq(tmass,fitParam[0],fitParam[1],91,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights91,binnweights91 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio91)

    ratio915=fiteq(tmass,fitParam[0],fitParam[1],91.5,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights915,binnweights915 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio915)
    sim_scale_factor905= np.sum(dataHist) / np.sum(simmassHistweights905)
    simHist_scaled905 = simmassHistweights905 * sim_scale_factor905
    
    sim_scale_factor91= np.sum(dataHist) / np.sum(simmassHistweights91) 
    simHist_scaled91 = simmassHistweights91 * sim_scale_factor91

    sim_scale_factor915 = np.sum(dataHist) / np.sum(simmassHistweights915) 
    simHist_scaled915 = simmassHistweights915 * sim_scale_factor915

    plt.plot(centers, simHist_scaled905, '-', linewidth=2, label='scaled-weighted 90.5 simulation')
    plt.plot(centers, simHist_scaled91, '-', linewidth=2, label='scaled-weighted 91 simulation')
    plt.plot(centers, simHist_scaled915, '-', linewidth=2, label='scaled-weighted 91.5 simulation')
    plt.title(f"Zreconstructed weight sim ({'real' if err==True else 'dont-use'})")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig(f"transient/Zreconstructed weight sim ({'real' if err==False else 'dont-use'}).pdf")
    plt.clf()
    plt.figure() # now for the chi ^2 plot of from the last 3  // now this will ahve to be recosntucted daat from real life not sim

    chi905=np.sum(((dataHist-simHist_scaled905)**2)/dataHist)
    chi91=np.sum(((dataHist-simHist_scaled91)**2)/dataHist)
    chi915=np.sum(((dataHist-simHist_scaled915)**2)/dataHist)

    plt.title("Z-chi")
    y=np.array([chi905,chi91,chi915])
    #print(y)
    x=np.array([90.5,91,91.5])
    plt.xlabel("target Zmass")
    plt.ylabel("Chi^2")

    plt.scatter(x,y)
    qfit = np.polyfit(x, y, deg=2)
    functionpfit = np.poly1d(qfit) 
    #print(qfit)
    xx = np.linspace((90.25), (91.75), 300)
    plt.plot(xx,functionpfit(xx))
    plt.savefig(f"transient/Z-chi({'real' if err==False else 'dont-use'}).pdf")

    min=-qfit[1]/(2*qfit[0])
    f=qfit[0]*min**2+qfit[1]*min+qfit[2]
    if err== False:   # if it is calibration with 0 error to get the values of mass
        print(fitParam)
        print(y)
        print(qfit)
        print("min ch^2 is",f)
        print("mz mass is", min)
        xerror_p=(-qfit[1]+np.sqrt(qfit[1]**2 -4*qfit[0]*(qfit[2]-1-f)))/(2*qfit[0])
        xerror_n=(-qfit[1]-np.sqrt(qfit[1]**2 -4*qfit[0]*(qfit[2]-1-f)))/(2*qfit[0])
    #print(xerror_p)
    #print(xerror_n)
    #print(min-xerror_p,min-xerror_n)
        print("the mass error is",min-xerror_n)
    plt.figure()
    return min

def Zstack_plot(tmass,simdatam,datam):
    c=C_rat #error in calibration
    dataHist,databinn,_d=plt.hist(conaeq(datam,1), bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-data-reconstructed",linewidth=1) #the data recosntuction data #for recosntructed simulation 
    trmassHist,binn,_t=plt.hist(tmass, bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-true",density=True,linewidth=1) #for true mass
    centers=0.5*(binn[1:]+binn[:-1])
    w=3
    m=91.1876
    fitParam,_tt = curve_fit(fiteq,centers,trmassHist,p0=[5,-0.7,m,w,0.4],bounds=([0.0,-1.0,60,0,0],[100.0,0,120,10,1]),maxfev=10000)# this uses the true mass
    print(fitParam)
    
    # ratio only has the t mass in it, the hsitogram suses this wieghtign witht the recosntructed in it
    #scaling the hisogram
    

    ratio905=fiteq(tmass,fitParam[0],fitParam[1],90.5,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights905,binnweights905 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio905)

    ratio91=fiteq(tmass,fitParam[0],fitParam[1],91,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights91,binnweights91 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio91)

    ratio915=fiteq(tmass,fitParam[0],fitParam[1],91.5,fitParam[3],fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights915,binnweights915 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio915)
    sim_scale_factor905= np.sum(dataHist) / np.sum(simmassHistweights905)
    simHist_scaled905 = simmassHistweights905 * sim_scale_factor905
    
    sim_scale_factor91= np.sum(dataHist) / np.sum(simmassHistweights91) 
    simHist_scaled91 = simmassHistweights91 * sim_scale_factor91

    sim_scale_factor915 = np.sum(dataHist) / np.sum(simmassHistweights915) 
    simHist_scaled915 = simmassHistweights915 * sim_scale_factor915

    plt.plot(centers, simHist_scaled905, '-', linewidth=2, label='scaled-weighted 90.5 simulation')
    plt.plot(centers, simHist_scaled91, '-', linewidth=2, label='scaled-weighted 91 simulation')
    plt.plot(centers, simHist_scaled915, '-', linewidth=2, label='scaled-weighted 91.5 simulation')

    plt.title("Z-reconstucted weighted sim")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig("transient/weights graph.pdf")
    plt.clf()

    plt.figure()
    plt.scatter(centers,dataHist, label='data ',color="black",s=10)
    plt.step(centers, simHist_scaled91, '-', linewidth=2,where='mid', label='scaled-weighted 91 simulation')
    plt.title("Z-data with 91 sim")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency ")
    plt.legend(loc='upper left',frameon=True, fontsize=8)
    
    plt.ylim(bottom=0)
    plt.savefig("transient/z-data-scatter.pdf")
    plt.show()
    #plt.clf()
    plt.figure()
    ratio_905_91=simHist_scaled905/simHist_scaled91
    ratio_91_91=simHist_scaled91/simHist_scaled91
    ratio_915_91=simHist_scaled915/simHist_scaled91
    ratio_data_91=dataHist/simHist_scaled91
    
    plt.scatter(centers,ratio_data_91,label='data ',color="black",s=10)
    plt.step(centers, ratio_905_91, '-', linewidth=2,where='mid', label='90.5-91 ratio')
    plt.step(centers, ratio_91_91, '-', linewidth=2,where='mid', label='91-91 ratio')
    plt.step(centers, ratio_915_91, '-', linewidth=2,where='mid', label='91.5-91 ratio')
    plt.title("Z-tempalte and data/template 91 GeV ratios")
    plt.xlabel("Mass_Gev")
    plt.ylabel("ratio")
    plt.legend(loc='upper left',frameon=True, fontsize=8)
    plt.savefig("transient/z-ratios.pdf")
    plt.show()

    #stack stuff
    dataerrors=np.sqrt(dataHist)
    plt.figure()
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 6),
                                   gridspec_kw={'hspace': 0.05})
    
    #top
   # ax1.scatter(centers,dataHist ,label='data ',color="black",s=10)
    ax1.errorbar(centers, dataHist, yerr=dataerrors,label='data ',color="black",fmt=".",markersize=2.5 )
    ax1.step(centers, simHist_scaled91, '-', linewidth=2,where='mid', label='scaled-weighted 91 simulation')
    ax1.set_ylabel("frequency")
    ax1.legend(loc='upper left',frameon=True, fontsize=8)

    # bottom
    ratioerror=ratio_data_91*dataerrors/dataHist
    #ax2.scatter(centers,ratio_data_91,label='data ',color="black",s=10)
    ax2.errorbar(centers, ratio_data_91, yerr=ratioerror,label='data/91Gev template ',color="black",fmt=".",markersize=2.5 )
    ax2.step(centers, ratio_905_91, '-', linewidth=2,where='mid', label='90.5-91 ratio')
    ax2.step(centers, ratio_91_91, '-', linewidth=2,where='mid', label='91-91 ratio')
    ax2.step(centers, ratio_915_91, '-', linewidth=2,where='mid', label='91.5-91 ratio')
    ax2.set_ylabel("ratio")
    ax2.set_xlabel("Dimuon mass GeV")
    ax2.legend(loc='upper left',frameon=True, fontsize=8)
    plt.savefig("transient/Z-stack.pdf")

    #plt.tight_layout()
    plt.show()

def Widthchi(tmass,simdatam,datam,calibrate,err):
    c=calibrate #error in calibration
    dataHist,databinn,_d=plt.hist(conaeq(datam,1), bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-data-reconstructed",linewidth=1) #the data recosntuction data #for recosntructed simulation 
    trmassHist,binn,_t=plt.hist(tmass, bins=np.linspace(80.0,100.0,50), histtype="step",label="Z-true",density=True,linewidth=1) #for true mass
    centers=0.5*(binn[1:]+binn[:-1])
    w=3
    m=91.1876
    fitParam,_tt = curve_fit(fiteq,centers,trmassHist,p0=[5,-0.7,m,w,0.4],bounds=([0.0,-1.0,60,0,0],[100.0,0,120,10,1]),maxfev=10000)# this uses the true mass
    #print(fitParam)
    
    # ratio only has the t mass in it, the hsitogram suses this wieghtign witht the recosntructed in it
    #scaling the hisogram
    
    
    ratio249=fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],1,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights249,binnweights249 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio249)

    ratio250=fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],2,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights250,binnweights250 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio250)

    ratio251=fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],3,fitParam[4])/fiteq(tmass,fitParam[0],fitParam[1],fitParam[2],fitParam[3],fitParam[4])
    simmassHistweights251,binnweights251 = np.histogram(conaeq(simdatam,c), bins=np.linspace(80.0,100.0,50), weights=ratio251)

    sim_scale_factor249= np.sum(dataHist) / np.sum(simmassHistweights249)
    simHist_scaled249 = simmassHistweights249 * sim_scale_factor249
    
    sim_scale_factor250= np.sum(dataHist) / np.sum(simmassHistweights250) 
    simHist_scaled250 = simmassHistweights250 * sim_scale_factor250

    sim_scale_factor251= np.sum(dataHist) / np.sum(simmassHistweights251) 
    simHist_scaled251 = simmassHistweights251 * sim_scale_factor251

    plt.plot(centers, simHist_scaled249, '-', linewidth=2, label='scaled-weighted 2.4 simulation')
    plt.plot(centers, simHist_scaled250, '-', linewidth=2, label='scaled-weighted 2.50 simulation')
    plt.plot(centers, simHist_scaled251, '-', linewidth=2, label='scaled-weighted 2.6 simulation')
    plt.title(f"Zreconstructed width sim ({'real' if err==True else 'dont-use'})")
    plt.xlabel("Mass_Gev")
    plt.ylabel("Frequency Density")
    plt.legend(loc='upper right')
    plt.savefig(f"transient/Zreconstructed width sim ({'real' if err==False else 'dont-use'}).pdf")
    plt.clf()
    plt.figure() # now for the chi ^2 plot of from the last 3  // now this will ahve to be recosntucted daat from real life not sim

    chi249=np.sum(((dataHist-simHist_scaled249)**2)/dataHist)
    chi250=np.sum(((dataHist-simHist_scaled250)**2)/dataHist)
    chi251=np.sum(((dataHist-simHist_scaled251)**2)/dataHist)
    plt.title("Z-chi")
    y=np.array([chi249,chi250,chi251])
    #print(y)
    x=np.array([1,2,3])
    plt.xlabel("target width")
    plt.ylabel("Chi^2")
    plt.scatter(x,y)
    qfit = np.polyfit(x, y, deg=2)
    functionpfit = np.poly1d(qfit) 
    #print(qfit)
    xx = np.linspace((0.8), (3.2), 300)
    plt.plot(xx,functionpfit(xx))
    plt.savefig(f"transient/Z-chi width({'real' if err==False else 'dont-use'}).pdf")

    min=-qfit[1]/(2*qfit[0])
    f=qfit[0]*min**2+qfit[1]*min+qfit[2]
    if err== False:   # if it is calibration with 0 error to get the values of mass
        print(fitParam)
        print(y)
        print(qfit)
        print("min ch^2 is",f)
        print("mz width is", min)
        xerror_p=(-qfit[1]+np.sqrt(qfit[1]**2 -4*qfit[0]*(qfit[2]-1-f)))/(2*qfit[0])
        xerror_n=(-qfit[1]-np.sqrt(qfit[1]**2 -4*qfit[0]*(qfit[2]-1-f)))/(2*qfit[0])
    #print(xerror_p)
    #print(xerror_n)
    #print(min-xerror_p,min-xerror_n)
        print("the width error is",min-xerror_n)
    plt.figure()
    return min



#TrueZfit_allplots(tmass,simdatam,datam)
#Zdata_with_fit(tmass,simdatam,datam)
#mikasuggestions(tmass,simdatam,datam)
#Zstack_plot(tmass,simdatam,datam)


#functions to run the mass and width and create the calirbation errors
def mass(tmass,simdatam,datam,Cr,Cre,err1,err2):
    e2=plotmass(tmass,simdatam,datam,Cre,err2) 
    e1=plotmass(tmass,simdatam,datam,Cr,err1)
    print("the calibration uncertancy in mass is ",np.sqrt((e1-e2)**2))   #so true means it has error calibration in it

mass(tmass,simdatam,datam,C_rat,C_rat+C_err,False,True)

def width(tmass,simdatam,datam,Cr,Cre,err1,err2):
    e2=Widthchi(tmass,simdatam,datam,Cre,err2) 
    e1=Widthchi(tmass,simdatam,datam,Cr,err1)
    print("the calibration uncertancy in width is ",e1-e2)
width(tmass,simdatam,datam,C_rat,C_rat+C_err,False,True)
'''
output = {
    "Chi_sq":"placeholder",
    "Z_mass": "placeholder",
    "Z_err": "placeholder"
    }
with open(f'Z_result_calibration_{args.Calibration}.json',"w") as OutputFile:
    dump(output,OutputFile,indent=2)
''' #Delete ''' to use this to output result to a json file
