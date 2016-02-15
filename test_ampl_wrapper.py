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

    scda.config_logging("wrapper_test.log")

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
    file_org = dict([('aperture dir', ampl_src_dir), ('intermediate mask dir', ampl_src_dir),\
                     ('ampl src dir', ampl_src_dir)]) # place holders

    # Design parameters, encoded in a hierarchical dictionary
    pupil_params = dict([('array diam', 200),('telescope diam', 12.)])
    FPM_params = dict([('spot rad',4), ('spot array size', 50)])
    LS_params = dict([('inner diam', 0.1), ('outer diam', 0.9)])
    image_params = dict([('logc', -10)])
    design_params = dict([('pupil_params', pupil_params), ('FPM_params', FPM_params),\
                          ('LS_params', LS_params), ('image_params', image_params)])
    # Options for constraints and optimizer
    solver_params = dict([('type', 'barrier')])
    
    test_coron = scda.HalfplaneAPLC(file_org=file_org, design_params=design_params, solver_params=solver_params)

    print test_coron.design_params['pupil_params']

    test_coron.write_ampl()
