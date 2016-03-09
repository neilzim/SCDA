#!/usr/bin/env python3
"""
Test the functionaility of the core SCDA 

02/14/2016 -- created by NTZ
"""

import scda
import numpy as np
import os
import sys

if __name__ == "__main__":

    scda.configure_log("wrapper_test.log")

    test_dir = "test_scda_aplc" # nominal destination for new AMPL programs
    #aux_dir = "~/SCDA/2d AMPL script - quarter pupil"
    aux_dir = "../2d AMPL script - quarter pupil"

    fileorg = {'work dir': test_dir, 'TelAp dir': aux_dir, 'FPM dir': aux_dir, 'LS dir': aux_dir,
               'TelAp fname': "CircPupil_N=0300_obs=20_center_quarter_spiders3=01_gaps=01.dat",
               'FPM fname': "CircPupil_N=0050_obs=00_center_quarter.dat",
               'LS fname': "CircPupil_N=0300_obs=40_center_quarter_spiders3=02.dat"}

    pupil_params = {'N': 300}
#    fpm_params = {'rad': 9.898/2, 'M':50}
#    fpm_params = {'rad': 6.466/2, 'M':50}
    fpm_params = {'rad': 8./2, 'M':50}
#    ls_params = {'id': 10, 'od': 0.9}
    ls_params = {}
    image_params = {'c': 10., 'iwa':3.5, 'owa':7., 'bw':0., 'Nlam':1}
    design_params = {'Pupil': pupil_params, 'FPM': fpm_params, 'LS': ls_params, 'Image': image_params}

#    solver_params = {'method': 'bar', 'presolve': False, 'Nthreads': 8}
    solver_params = {}
    
    atlast_coron = scda.QuarterplaneAPLC(fileorg=fileorg, design=design_params, solver=solver_params, verbose=True)

    atlast_coron.write_ampl(ampl_src_fname="ref_qpaplc_dev2.mod", overwrite=True)
