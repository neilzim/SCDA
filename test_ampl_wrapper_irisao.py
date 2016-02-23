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

    test_dir = "test_scda_aplc" # nominal destination for new AMPL programs
    aux_dir = "~/SCDA/2d AMPL script - half pupil"

    fileorg = {'work dir': test_dir, 'TelAp dir': aux_dir, 'FPM dir': aux_dir, 'LS dir': aux_dir,
               'TelAp fname': "IRISAO_N=0150_center_half_spiders3=01_gapID=10_BW.dat",
               'FPM fname': "CircPupil_N=0050_obs=00_center_half.dat",
               'LS fname': "IRISAO-0_N=0150_center_half_spiders3=02_ID=20_OD=098.dat"}

    pupil_params = {'N': 150}
#    fpm_params = {'rad': 9.898/2, 'M':50}
    fpm_params = {'rad': 6.466/2, 'M':50}
#    ls_params = {'id': 10, 'od': 0.9}
    ls_params = {}
    image_params = {'c': 8., 'iwa':3., 'owa':10.}
    design_params = {'Pupil': pupil_params, 'FPM': fpm_params, 'LS': ls_params, 'Image': image_params}

#    solver_params = {'method': 'bar', 'presolve': False, 'Nthreads': 8}
    solver_params = {}
    
    irisao_coron = scda.HalfplaneAPLC(fileorg=fileorg, design=design_params, solver=solver_params, verbose=True)

    irisao_coron.write_ampl(ampl_src_fname="irisao_hpaplc_C08.mod", overwrite=True)
