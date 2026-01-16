#!/usr/bin/env bash

./run ./scripts/Upsilon.py --FullOutput=TRUE
./run ./scripts/ZAttempt3.py
./run ./scripts/ZAttempt3.py --Calibration=TRUE
./run ./scripts/Z_Simaltaneous_fits.py --Calibration=TRUE