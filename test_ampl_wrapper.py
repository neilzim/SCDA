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

    # IrisAO design parameters
    pupil_params = {'N': 200,'tel diam': 12.}
    fpm_params = {'rad': 9.898/2, 'M':50}
    ls_params = {'id': 10, 'od': 0.9}
    image_params = {'c': 10., 'iwa':3.5, 'owa':10., 'owa2':25.}
    design_params = {'Pupil': pupil_params, 'FPM': fpm_params, 'LS': ls_params, 'Image': image_params}
    # Options for constraints and optimizer
    solver_params = {'method': 'bar', 'presolve': False}
    
    irisao_coron = scda.HalfplaneAPLC(fileorg=fileorg, design=design_params, solver=solver_params)

    irisao_coron.logger.info("file organization: {0}".format(irisao_coron.fileorg))
    irisao_coron.logger.info("solver parameters: {0}".format(irisao_coron.solver))
    irisao_coron.logger.info("design parameters: {0}".format(irisao_coron.design))

    print
#    print("AMPL source file name: {0}".format(test_coron.ampl_src_fname))

    irisao_coron.write_ampl()
