#!/usr/bin/env bash

./scripts/Upsilon.py --FullOutput="TRUE"
./scripts/Upsilon.py --Compare="TRUE"

./scripts/Z_Simaltaneous_fits.py --Calibration="TRUE" --Smear_error=PLUS &
./scripts/Z_Simaltaneous_fits.py --Calibration="TRUE" --Smear_error=MINUS &
./scripts/Z_Simaltaneous_fits.py --Calibration="TRUE" --C_ratio=PLUS &
./scripts/Z_Simaltaneous_fits.py --Calibration="TRUE" --C_ratio=MINUS &
./scripts/Z_Simaltaneous_fits.py --Calibration="TRUE" &
wait

./scripts/Z_mass_width_calibration_error.py
pdflatex Images.tex