#!/usr/bin/env python
import numpy as np

fname1 = "../mphys-electroweak-project/smear_minus.txt"
Sm = np.loadtxt(fname1, dtype=float)
print(Sm)

fname2 = "../mphys-electroweak-project/smear_plus.txt"
Sp = np.loadtxt(fname2, dtype=float)
print(Sp)

smear_mass_error=(Sm[0]-Sp[0])/2
smear_width_error=(Sm[1]-Sp[1])/2

fname3 = "../mphys-electroweak-project/C_rat_minus.txt"
Sm = np.loadtxt(fname3, dtype=float)
print(Sm)


fname4 = "../mphys-electroweak-project/C_rat_plus.txt"
Sp = np.loadtxt(fname4, dtype=float)
print(Sp)

Calibration_ratio_mass_error=(Sm[0]-Sp[0])/2
Calibration_ratio_width_error=(Sm[1]-Sp[1])/2

print(smear_mass_error)
print(smear_width_error)
print(Calibration_ratio_mass_error)
print(Calibration_ratio_width_error)
with open("results_table.tex", "w") as f:
    f.write(r"\begin{table}[H]" + "\n")
    f.write(r"\centering" + "\n")
    f.write(r"\begin{tabular}{l c c}" + "\n")
    f.write(r"\hline" + "\n")
    f.write(r" & Mass & Width \\" + "\n")
    f.write(r"\hline" + "\n")
    f.write(f"Calibration error & {Calibration_ratio_mass_error:.5f} & {Calibration_ratio_width_error:.5f} \\\\\n")
    f.write(f"Smearing error & {smear_mass_error:.5f} & {smear_width_error:.5f} \\\\\n")
    f.write(r'$\Upsilon$(1S) mass error & 0.0001 & 0.0000 \\' + '\n')
    f.write(r"\hline" + "\n")
    f.write(r"\end{tabular}" + "\n")
    f.write(r"\caption{Z Systematic error on mass and width}" + "\n")
    f.write(r"\end{table}" + "\n")

fname5 = "../mphys-electroweak-project/mass-width_values_and_error.txt"
mass_width = np.loadtxt(fname5, dtype=float)
mass=mass_width[0]
width=mass_width[1]
mass_error=mass_width[2]
width_error=mass_width[3]
chi2=mass_width[4]
ndf=mass_width[5]
corelation=np.sqrt((mass_width[6])**2)

mass_syst_err=np.sqrt(Calibration_ratio_mass_error**2 +smear_mass_error**2)
width_syst_err=np.sqrt(Calibration_ratio_width_error**2 +smear_width_error**2)
with open("mass_width_results_table.tex", "w") as f:
    f.write(r"\begin{table}[H]" + "\n")
    f.write(r"\centering" + "\n")
    f.write(r"\begin{tabular}{l c c}" + "\n")
    f.write(r"\hline" + "\n")
    f.write(r" & Mass & Width \\" + "\n")
    f.write(r"\hline" + "\n")
    f.write(f"Values & {mass:.5f} & {width:.5f} \\\\\n")
    f.write(f"Statistical error & {mass_error:.5f} & {width_error:.5f} \\\\\n")
    f.write(r"\hline" + "\n")
    f.write(r"\end{tabular}" + "\n")
    f.write(r"\caption{Z Mass and width results with statistical error}" + "\n")
    f.write(r"\end{table}" + "\n")


with open("CHI2_ndf.tex", "w") as f:
    f.write(r"\begin{table}[H]" + "\n")
    f.write(r"\centering" + "\n")
    f.write(r"\begin{tabular}{cc}" + "\n")
    f.write(r"\hline" + "\n")
    f.write(r"$\chi^2$ & ndf \\" + "\n")
    f.write(r"\hline" + "\n")
    f.write(f"{chi2:.2f} & {ndf:.0f} \\\\\n")
    f.write(r"\hline" + "\n")
    f.write(r"\end{tabular}" + "\n")
    f.write(r"\caption{minimum $\chi^2$ and number of degrees of freedom}" + "\n")
    f.write(r"\end{table}" + "\n")

with open("constants.tex", "w") as f:
    f.write(r"\newcommand{\ZMass}{%.3f}" % mass + "\n")
    f.write(r"\newcommand{\ZWidth}{%.3f}" % width + "\n")
    f.write(r"\newcommand{\ZMasssysterror}{%.3f}" % mass_syst_err + "\n")
    f.write(r"\newcommand{\ZWidthsysterror}{%.3f}" % width_syst_err + "\n")
    f.write(r"\newcommand{\Zmassstaterror}{%.3f}" % mass_error + "\n")
    f.write(r"\newcommand{\Zwidthstaterror}{%.3f}" % width_error + "\n")
    f.write(r"\newcommand{\ZCorrelation}{%.3f}" % corelation + "\n")



fname7 = "../mphys-electroweak-project/pos_dipole"
mass_width = np.loadtxt(fname7, dtype=float)
massp=mass_width[0]
widthp=mass_width[1]
mass_errorp=mass_width[2]
width_errorp=mass_width[3]
chi2p=mass_width[4]
ndfp=mass_width[5]
corelationp=np.sqrt((mass_width[6])**2)

fname8 = "../mphys-electroweak-project/neg_dipole"
mass_width = np.loadtxt(fname8, dtype=float)
massn=mass_width[0]
widthn=mass_width[1]
mass_errorn=mass_width[2]
width_errorn=mass_width[3]
chi2n=mass_width[4]
ndfn=mass_width[5]
corelationn=np.sqrt((mass_width[6])**2)

differnce_mass = np.sqrt((massp-massn)**2)
differnce_width = np.sqrt((widthp-widthn)**2)

magnet_mass_significance1=differnce_mass/(np.sqrt(mass_errorp**2 + mass_errorn**2))
magnet_width_significance1=differnce_width/(np.sqrt(width_errorp**2 + width_errorn**2))

with open("dipole_table.tex", "w") as f:
    f.write(r"\begin{table}[H]" + "\n")
    f.write(r"\centering" + "\n")
    f.write(r"\begin{tabular}{lccc}" + "\n")
    f.write(r"\hline" + "\n")
    f.write(r"Measuremnt & Positive dipole & Negative dipole & Significance\\" + "\n")
    f.write(r"\hline" + "\n")

    f.write(f"Mass [GeV] & ${massp:.5f} \pm {mass_errorp:.2f}$ & "
            f"${massn:.5f} \pm {mass_errorn:.2f}$ & "
            f"${magnet_mass_significance1:.5f} $\\\\\n")

    f.write(f"Width [GeV] & ${widthp:.5f} \pm {width_errorp:.2f}$ & "
            f"${widthn:.5f} \pm {width_errorn:.2f}$ & "
            f"${magnet_width_significance1:.5f} $ \\\\\n")

    f.write(f"$\\chi^2$ & {chi2p:.2f} & {chi2n:.2f} & N/A\\\\\n")

    f.write(f"ndf & {ndfp:.0f} & {ndfn:.0f} & N/A\\\\\n")

    f.write(r"\hline" + "\n")
    f.write(r"\end{tabular}" + "\n")
    f.write(r"\caption{Reults for positive and negative dipoles}" + "\n")
    f.write(r"\end{table}" + "\n")


fname9 = "../mphys-electroweak-project/0__π"
mass_width = np.loadtxt(fname9, dtype=float)
massp=mass_width[0]
widthp=mass_width[1]
mass_errorp=mass_width[2]
width_errorp=mass_width[3]
chi2p=mass_width[4]
ndfp=mass_width[5]
corelationp=np.sqrt((mass_width[6])**2)

fname10 = "../mphys-electroweak-project/-0__π"
mass_width = np.loadtxt(fname10, dtype=float)
massn=mass_width[0]
widthn=mass_width[1]
mass_errorn=mass_width[2]
width_errorn=mass_width[3]
chi2n=mass_width[4]
ndfn=mass_width[5]
corelationn=np.sqrt((mass_width[6])**2)

differnce_mass = np.sqrt((massp-massn)**2)
differnce_width = np.sqrt((widthp-widthn)**2)

magnet_mass_significance2=differnce_mass/(np.sqrt(mass_errorp**2 + mass_errorn**2))
magnet_width_significance2=differnce_width/(np.sqrt(width_errorp**2 + width_errorn**2))
with open("Angle_table.tex", "w") as f:
    f.write(r"\begin{table}[H]" + "\n")
    f.write(r"\centering" + "\n")
    f.write(r"\begin{tabular}{lccc}" + "\n")
    f.write(r"\hline" + "\n")
    f.write(r"Measuremnt & angle range 0-$\pi$ & angle range-$\pi$-0 & Significance\\" + "\n")
    f.write(r"\hline" + "\n")

    f.write(f"Mass [GeV] & ${massp:.5f} \pm {mass_errorp:.2f}$ & "
            f"${massn:.5f} \pm {mass_errorn:.2f}$ & "
            f"${magnet_mass_significance2:.5f} $\\\\\n")

    f.write(f"Width [GeV] & ${widthp:.5f} \pm {width_errorp:.2f}$ & "
            f"${widthn:.5f} \pm {width_errorn:.2f}$ & "
            f"${magnet_width_significance2:.5f} $ \\\\\n")

    f.write(f"$\\chi^2$ & {chi2p:.2f} & {chi2n:.2f} & N/A\\\\\n")

    f.write(f"ndf & {ndfp:.0f} & {ndfn:.0f} & N/A\\\\\n")

    f.write(r"\hline" + "\n")
    f.write(r"\end{tabular}" + "\n")
    f.write(r"\caption{Reults for angles between 0 and $\pi$ and -$\pi$ to 0}" + "\n")
    f.write(r"\end{table}" + "\n")

fname11 = "../mphys-electroweak-project/Lower_Pt"
mass_width = np.loadtxt(fname11, dtype=float)
massp=mass_width[0]
widthp=mass_width[1]
mass_errorp=mass_width[2]
width_errorp=mass_width[3]
chi2p=mass_width[4]
ndfp=mass_width[5]
corelationp=np.sqrt((mass_width[6])**2)

fname12 = "../mphys-electroweak-project/Higher_Pt"
mass_width = np.loadtxt(fname12, dtype=float)
massn=mass_width[0]
widthn=mass_width[1]
mass_errorn=mass_width[2]
width_errorn=mass_width[3]
chi2n=mass_width[4]
ndfn=mass_width[5]
corelationn=np.sqrt((mass_width[6])**2)

differnce_mass = np.sqrt((massp-massn)**2)
differnce_width = np.sqrt((widthp-widthn)**2)

magnet_mass_significance3=differnce_mass/(np.sqrt(mass_errorp**2 + mass_errorn**2))
magnet_width_significance3=differnce_width/(np.sqrt(width_errorp**2 + width_errorn**2))
with open("Dimuon-transverse-P.tex", "w") as f:
    f.write(r"\begin{table}[H]" + "\n")
    f.write(r"\centering" + "\n")
    f.write(r"\begin{tabular}{lccc}" + "\n")
    f.write(r"\hline" + "\n")
    f.write(r"Measuremnt & 0-18.3 Gev & $>$ 18.3Gev  & Significance\\" + "\n")
    f.write(r"\hline" + "\n")

    f.write(f"Mass [GeV] & ${massp:.5f} \pm {mass_errorp:.2f}$ & "
            f"${massn:.5f} \pm {mass_errorn:.2f}$ & "
            f"${magnet_mass_significance3:.5f} $\\\\\n")

    f.write(f"Width [GeV] & ${widthp:.5f} \pm {width_errorp:.2f}$ & "
            f"${widthn:.5f} \pm {width_errorn:.2f}$ & "
            f"${magnet_width_significance3:.5f} $ \\\\\n")

    f.write(f"$\\chi^2$ & {chi2p:.2f} & {chi2n:.2f} & N/A\\\\\n")

    f.write(f"ndf & {ndfp:.0f} & {ndfn:.0f} & N/A\\\\\n")

    f.write(r"\hline" + "\n")
    f.write(r"\end{tabular}" + "\n")
    f.write(r"\caption{Reults for Dimuon transverse moemntum splitting up data at 18.3 Gev}" + "\n")
    f.write(r"\end{table}" + "\n")

fname13 = "../mphys-electroweak-project/30-bins"
mass_width = np.loadtxt(fname12, dtype=float)
massp=mass_width[0]
widthp=mass_width[1]
mass_errorp=mass_width[2]
width_errorp=mass_width[3]
chi2p=mass_width[4]
ndfp=mass_width[5]
corelationp=np.sqrt((mass_width[6])**2)


fname12 = "../mphys-electroweak-project/80-bins"
mass_width = np.loadtxt(fname12, dtype=float)
massn=mass_width[0]
widthn=mass_width[1]
mass_errorn=mass_width[2]
width_errorn=mass_width[3]
chi2n=mass_width[4]
ndfn=mass_width[5]
corelationn=np.sqrt((mass_width[6])**2)
with open("Bins_table.tex", "w") as f:
    f.write(r"\begin{table}[H]" + "\n")
    f.write(r"\centering" + "\n")
    f.write(r"\begin{tabular}{l c c}" + "\n")
    f.write(r"\hline" + "\n")
    f.write(r" & Mass & Width \\" + "\n")
    f.write(r"\hline" + "\n")
    f.write(f"30 bins & {massp:.5f} $\pm$ {mass_errorp:.5f} & {widthp:.5f} $\pm$ {width_errorp:.5f} \\\\\n")
    f.write(f"50 bins & {mass:.5f} $\pm$ {mass_error:.5f} & {width:.5f} $\pm$ {width_error:.5f} \\\\\n")
    f.write(f"80 bins & {massn:.5f} $\pm$ {mass_errorn:.5f} & {widthn:.5f} $\pm$ {width_errorn:.5f} \\\\\n")
    f.write(r"\hline" + "\n")
    f.write(r"\end{tabular}" + "\n")
    f.write(r"\caption{Z mass and width with statistical error for each bin size}" + "\n")
    f.write(r"\end{table}" + "\n")


fname13 = "../mphys-electroweak-project/seperate_Pt_selection"
mass_width = np.loadtxt(fname13, dtype=float)
mass1=mass_width[0]
width1=mass_width[1]
mass_error1=mass_width[2]
width_error1=mass_width[3]
chi21=mass_width[4]
ndf1=mass_width[5]
corelation1=np.sqrt((mass_width[6])**2)

fname14 = "../mphys-electroweak-project/psuedo"
mass_width = np.loadtxt(fname14, dtype=float)
mass2=mass_width[0]
width2=mass_width[1]
mass_error2=mass_width[2]
width_error2=mass_width[3]
chi22=mass_width[4]
ndf2=mass_width[5]
corelation2=np.sqrt((mass_width[6])**2)

fname15 = "../mphys-electroweak-project/seperate_Pt_psuedo_selection"
mass_width = np.loadtxt(fname15, dtype=float)
mass3=mass_width[0]
width3=mass_width[1]
mass_error3=mass_width[2]
width_error3=mass_width[3]
chi23=mass_width[4]
ndf3=mass_width[5]
corelation3=np.sqrt((mass_width[6])**2)



with open("selection_table.tex", "w") as f:
    f.write(r"\begin{table}[H]" + "\n")
    f.write(r"\centering" + "\n")
    f.write(r"\begin{tabular}{l c c c}" + "\n")  # 4 columns now
    f.write(r"\hline" + "\n")
    f.write(r" & Transverse momentum & pseudorapidity & Both applied simultaneously \\" + "\n")
    f.write(r"\hline" + "\n")
    f.write(f"mass & {mass1:.5f} $\pm$ {mass_error1:.5f} & {mass2:.5f} $\pm$ {mass_error2:.5f} & {mass3:.5f} $\pm$ {mass_error3:.5f} \\\\\n")
    f.write(f"width & {width1:.5f} $\pm$ {width_error1:.5f} & {width2:.5f} $\pm$ {width_error2:.5f} & {width3:.5f} $\pm$ {width_error3:.5f} \\\\\n")
    f.write(f"$\chi^2$ & {chi21:.3f}  & {chi22:.3f} & {chi23:.3f} \\\\\n")
    f.write(f"$\chi^2$ & {ndf1:.3f}  & {ndf2:.3f} & {ndf3:.3f} \\\\\n")
    f.write(f"corelation & {corelation1:.3f}  & {corelation2:.3f} & {corelation3:.3f} \\\\\n")
    f.write(r"\hline" + "\n")
    f.write(r"\end{tabular}" + "\n")
    f.write(r"\caption{Transver momentum selection of $>$ 20 Gev per in both muons and pseudorapidity selection 2.2$<\eta <$4.4}" + "\n")
    f.write(r"\end{table}" + "\n")


with open("significance_table.tex", "w") as f:
    f.write(r"\begin{table}[H]" + "\n")
    f.write(r"\centering" + "\n")
    f.write(r"\begin{tabular}{l c c c}" + "\n")  # 4 columns now
    f.write(r"\hline" + "\n")
    f.write(r" & magnetic dipole & Dimuon $P_T$ & azimuthal angle \\" + "\n")
    f.write(r"\hline" + "\n")
    f.write(f"Significance on the mass & {magnet_mass_significance1:.5f}  & {magnet_mass_significance3:.5f}  & {magnet_mass_significance2:.5f}  \\\\\n")
    f.write(f"Significance on the width & {magnet_width_significance1:.5f}  & {magnet_width_significance3:.5f}  & {magnet_mass_significance2:.5f}  \\\\\n")
    f.write(r"\hline" + "\n")
    f.write(r"\end{tabular}" + "\n")
    f.write(r"\caption{Significance values for three validation checks}" + "\n")
    f.write(r"\end{table}" + "\n")