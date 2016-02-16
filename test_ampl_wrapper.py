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

#    ampl_src_dir = os.path.join(top_path, "test_scda_aplc") # nominal destination for new AMPL programs
    ampl_src_dir = os.path.join(os.getcwd(), "test_scda_aplc_new") # nominal destination for new AMPL programs
    fileorg = {'ampl src dir': ampl_src_dir}

    # Design parameters, encoded in a nested dictionary
    pupil_params = {'N': 200,'tel diam': 12.}
    fpm_params = {'rad': 4.}
    ls_params = {'id': 10, 'od': 0.9}
    image_params = {'c': 10.1}
    design_params = {'Pupil': pupil_params, 'FPM': fpm_params, 'LS': ls_params, 'Image': image_params}
    # Options for constraints and optimizer
    solver_params = {'method': 'bar', 'presolve': False}
    
    test_coron = scda.HalfplaneAPLC(fileorg=fileorg, design=design_params, solver=solver_params)

    print("file organization: {0}".format(test_coron.fileorg))
    print("solver parameters: {0}".format(test_coron.solver))
    print("design parameters: {0}".format(test_coron.design))

    print
#    print("AMPL source file name: {0}".format(test_coron.ampl_src_fname))

    test_coron.write_ampl()
