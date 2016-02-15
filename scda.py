"""
Core definitions for STScI's
Segmented Coronagraph Design & Analysis investigation

02/14/2016 -- created by NTZ
"""

import os
import sys
import logging
import datetime

def config_logging(log_fname=None):
    logger = logging.getLogger("scda.logger")
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout) # console dump
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)
    if log_fname is not None: # optional file log in addition to the console dump
        fh = logging.FileHandler(log_fname, mode="w")
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)

class LyotCoronagraph(object): # Lyot coronagraph superclass
    def __init__(self, *args, **kwargs):
        self.logger = logging.getLogger('scda.logger')
        self.file_org = kwargs.get('file_org')
        self.design_params = kwargs.get('design_params')
        self.solver_params = kwargs.get('solver_params')
        self.eval_params = kwargs.get('eval_params')
        self.ampl_fname = "lyot_coronagraph.mod"
    def write_ampl(self):
        logger.info("Writing the AMPL program")
    def read_solution(self):
        logger.info("Reading in the apodizer solution and parse the optimizer log")
    def create_eval_model(self):
        logger.info("Defining a python/Poppy model to evaluate the solution")
    def eval_solution(self):
        logger.info("Evaluating the design throughput, PSF FWHM, etc., and writing summary")

class NdiayeAPLC(LyotCoronagraph): # Image-constrained APLC following N'Diaye et al. (2015, 2016)
    def __init__(self, *args, **kwargs):
        super(NdiayeAPLC, self).__init__(self, *args, **kwargs)
        self.ampl_fname = "ndiaye_aplc.mod"

class HalfplaneAPLC(NdiayeAPLC): # Subclass for the half-plane symmetry case
    def __init__(self, *args, **kwargs):
        super(HalfplaneAPLC, self).__init__(self, *args, **kwargs)
        self.ampl_fname = "halfplane_aplc.mod"
    def write_ampl(self): 
        self.logger.info("Writing the AMPL program for this APLC with halfplane symmetry")
        self.abs_ampl_fname = os.path.join(self.file_org['ampl src dir'], self.ampl_fname)
        if os.path.exists(self.abs_ampl_fname):
            self.logger.warning("Overwriting the existing copy of %s"%self.abs_ampl_fname)
        mod_fobj = open(self.abs_ampl_fname, "w")
        mod_fobj.write("""# AMPL program to optimize an APLC with half-plane symmetry
# created with %s at %s
load amplgsl.dll;
param pi:= 4*atan(1);
""" % (os.path.basename(__file__), datetime.datetime.now()))
        mod_fobj.close()
        self.logger.info("Wrote %s"%self.abs_ampl_fname)
