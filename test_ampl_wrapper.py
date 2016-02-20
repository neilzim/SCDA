#!/usr/bin/env python3
"""
Test the functionaility of the core SCDA 

02/14/2016 -- created by NTZ
"""

import scda
import numpy as np
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":

    scda.configure_log("wrapper_test.log")

    test_dir = "test_aplc_wrapper" # nominal destination for new AMPL programs
    aux_dir = "~/STScI/SCDA/2d AMPL script - half pupil"

    fileorg = {'work dir': test_dir, 'TelAp dir': aux_dir, 'FPM dir': aux_dir, 'LS dir': aux_dir, \
               'TelAp fname': "IRISAO_N=0200_center_half_spiders3=01_gapID=10_BW.dat", \
               'FPM fname': "CircPupil_N=0050_obs=00_center_half.dat", \
               'LS fname': "IRISAO-0_N=0200_center_half_spiders3=02_ID=20_OD=098.dat"}

    # IrisAO design parameters
    pupil_params = {'N': 200,'tel diam': 12.}
    fpm_params = {'rad': 9.898/2, 'M':50}
    ls_params = {}
    image_params = {'c': 10., 'iwa':3.5, 'owa':10.}
    design_params = {'Pupil': pupil_params, 'FPM': fpm_params, 'LS': ls_params, 'Image': image_params}

    # Options for constraints and optimizer
    solver_params = {'method': 'bar', 'presolve': False, 'Nthreads': 8}

    irisao_coron = scda.HalfplaneAPLC( fileorg=fileorg, design=design_params, solver=solver_params )

    irisao_coron.write_ampl(ampl_src_fname="test_aplc.mod", overwrite=True)
