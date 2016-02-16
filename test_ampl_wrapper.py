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

    top_path = os.getcwd() # During a "real" survey preparation, this will be an agreed
                           # location on central store
    assert os.path.exists(top_path)
    ampl_src_dir = os.path.join(top_path, "test_scda_aplc") # nominal destination for new AMPL programs
    if not os.path.exists(ampl_src_dir):
        os.mkdir(ampl_src_dir) # TODO: set context-dependent group ID and permissions
    # The file_org structure holds general locations of telescope apertures,
    # intermediate masks, and co-eval AMPL programs.
    # File names specific to a coronagraph design (AMPL program, apodizer solution, optimizer log, etc)
    # are set downstream as object attributes.
    file_org = dict([('tel ap dir', ampl_src_dir), ('FPM dir', ampl_src_dir), ('LS dir', ampl_src_dir),\
                     ('ampl src dir', ampl_src_dir)]) # place holders

    # Design parameters, encoded in a hierarchical dictionary
    pupil_params = dict([('N', 200),('tel diam', 12.)])
    fpm_params = dict([('rad',4), ('M', 50)])
    ls_params = dict([('id', 0.1), ('od', 0.9)])
    image_params = dict([('c', -10)])
    design_params = dict([('Pupil', pupil_params), ('FPM', fpm_params),\
                          ('LS', ls_params), ('Image', image_params)])
    # Options for constraints and optimizer
    solver_params = dict([('method', 'barrier')])
    
    test_coron = scda.HalfplaneAPLC(fileorg=file_org, design=design_params, solver=solver_params)

    print("file organization: {0}".format(test_coron.fileorg))
    print("solver parameters: {0}".format(test_coron.solver))
    print("design parameters: {0}".format(test_coron.design))

    test_coron.write_ampl()
