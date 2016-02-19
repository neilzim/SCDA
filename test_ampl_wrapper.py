#!/usr/bin/env python3
"""
Test the functionaility of the core SCDA 

02/14/2016 -- created by NTZ
"""

import scda
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

if __name__ == "__main__":

    scda.configure_log("wrapper_test.log")

    local_test_dir = os.path.join(os.getcwd(), "test_scda_aplc") # nominal destination for new AMPL programs
    aux_dir = "../2d AMPL script - half pupil"
    irisao_telap_fname = os.path.join(aux_dir, "IRISAO_N=0200_center_half_spiders3=01_gapID=10_BW.dat")
    occspot_fname = os.path.join(aux_dir, "CircPupil_N=0050_obs=00_center_half.dat")
    irisao_lyotstop_fname = os.path.join(aux_dir, "IRISAO-0_N=0200_center_half_spiders3=02_ID=20_OD=098.dat")

    fileorg = { 'work dir': local_test_dir, 'TelAp fname': irisao_telap_fname, \
                'FPM dir': occspot_fname, 'LS fname': irisao_lyotstop_fname }

    # half-plane IRISAO example design parameters
    pupil_params = {'N': 200,'tel diam': 12.}
    fpm_params = {'rad': 9.898/2, 'M':50}
#    ls_params = {'id': 10, 'od': 0.9}
    ls_params = {}
    image_params = {'c': 10., 'iwa':3.5, 'owa':10., 'oca':25.}
    design_params = {'Pupil': pupil_params, 'FPM': fpm_params, 'LS': ls_params, 'Image': image_params}

    # Options for constraints and optimizer
#    solver_params = {'method': 'bar', 'presolve': False, 'Nthreads': 8}
    solver_params = {}
    
    irisao_coron = scda.HalfplaneAPLC( fileorg=fileorg, design=design_params, solver=solver_params )

    irisao_coron.write_ampl()
