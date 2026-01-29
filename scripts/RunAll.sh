#!/usr/bin/env bash

./run ./scripts/Upsilon.py --FullOutput="TRUE"
./run ./scripts/Upsilon.py --Compare="TRUE"
./run ./scripts/Z_Simaltaneous_fits.py --Calibration="TRUE" --Smear_error=PLUS 
./run ./scripts/Z_Simaltaneous_fits.py --Calibration="TRUE" --Smear_error=MINUS 
./run ./scripts/Z_Simaltaneous_fits.py --Calibration="TRUE" --C_ratio=PLUS 
./run ./scripts/Z_Simaltaneous_fits.py --Calibration="TRUE" --C_ratio=MINUS 

./run ./scripts/Z_Simaltaneous_fits.py --Calibration="TRUE" 
./run ./scripts/Z_mass_width_calibration_error.py
./run pdflatex Images.tex