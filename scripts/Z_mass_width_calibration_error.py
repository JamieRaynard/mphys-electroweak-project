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

fname3 = "../mphys-electroweak-project/smear_minus.txt"
Sm = np.loadtxt(fname3, dtype=float)
print(Sm)


fname4 = "../mphys-electroweak-project/smear_plus.txt"
Sp = np.loadtxt(fname4, dtype=float)
print(Sp)

Calibration_ratio_mass_error=(Sm[0]-Sp[0])/2
Calibration_ratio_width_error=(Sm[1]-Sp[1])/2

print(smear_mass_error)
print(smear_width_error)
print(Calibration_ratio_mass_error)
print(Calibration_ratio_width_error)