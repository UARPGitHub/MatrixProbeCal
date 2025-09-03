# 2D Matrix Probe Calibration

This repository contains an extract of the code used for performing calibration
of a 2D matrix ultrasound probe as used for IUS 2025 Paper ID #3559 

    "2D Ultrasound Transducer Sub-Array Misalignment Correction Using Chirp Transmission and Matched Filtering"

The main function is `MatrixProbeCal.m` which will need to be modified to suit
custom systems and lab equipment control.

The original code used for this paper is integrated into the UToolbox distribution
which forms the control libraries for the University of Leeds UARP system. 

The code in this repository has been split out for the benefit of users of other
ultrasound systems. While nominally correct, the code in this repository has not
been tested and should be used at your own risk. Please validate the calculated
coordinates for accuracy yourself before relying on them for experiments.

The code is released as open source under a GPL3 license. No warranty or fitness
for purpose is provided, implied or otherwise.
