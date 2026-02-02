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
    f.write(r"\begin{tabular}{l c c}" + "\n")
    f.write(r"\hline" + "\n")
    f.write(r" & Mass & Width \\" + "\n")
    f.write(r"\hline" + "\n")
    f.write(f"Calibration error & {Calibration_ratio_mass_error:f} & {Calibration_ratio_width_error:f} \\\\\n")
    f.write(f"Smearing error & {smear_mass_error:f} & {smear_width_error:f} \\\\\n")
    f.write(r"\hline" + "\n")
    f.write(r"\end{tabular}" + "\n")

fname5 = "../mphys-electroweak-project/mass-width_values_and_error.txt"
mass_width = np.loadtxt(fname5, dtype=float)
mass=mass_width[0]
width=mass_width[1]
mass_error=mass_width[2]
width_error=mass_width[3]
chi2_ndf=mass_width[4]
with open("mass_width_results_table.tex", "w") as f:
    f.write(r"\begin{tabular}{l c c}" + "\n")
    f.write(r"\hline" + "\n")
    f.write(r" & Mass & Width \\" + "\n")
    f.write(r"\hline" + "\n")
    f.write(f"values & {mass:f} & {width:f} \\\\\n")
    f.write(f"statistical error & {mass_error:f} & {width_error:f} \\\\\n")
    f.write(r"\hline" + "\n")
    f.write(r"\end{tabular}" + "\n")



with open("CHI2_ndf.tex", "w") as f: 
    f.write(r"\begin{tabular}{c}" + "\n")  
    f.write(r"\hline" + "\n") 
    f.write(r"CHI$^2$/ndf \\" + "\n") 
    f.write(r"\hline" + "\n") 
    f.write(f"{chi2_ndf:f} \\\\" + "\n")  
    f.write(r"\hline" + "\n") 
    f.write(r"\end{tabular}" + "\n")
