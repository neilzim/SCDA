"""
Core definitions of python tools for the STScI
Segmented Coronagraph Design & Analysis investigation

02/14/2016 -- created by NTZ
"""

import os
import sys
import logging
import datetime
import textwrap

def configure_log(log_fname=None):
    logger = logging.getLogger("scda.logger")
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout) # console dump
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)
    if log_fname is not None: # optional file log in addition to the console dump
        fh = logging.FileHandler(log_fname, mode="w")
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)

class LyotCoronagraph(object): # Lyot coronagraph base class
    _key_fields = dict([('fileorg', ['tel ap dir', 'FPM dir', 'LS dir', 'ampl src dir', 'sol dir', 'eval dir']),\
                        ('solver', ['constr', 'method', 'presolve'])])

    def __init__(self, **kwargs):
        self.logger = logging.getLogger('scda.logger')
     
        # Only set fileorg and solver attributes here, since design and eval parameter checking is design-specific. 
        setattr(self, 'fileorg', {})
        if 'fileorg' in kwargs:
            for field, location in kwargs['fileorg'].items():
                if field in self._key_fields['fileorg']:
                    self.fileorg[field] = location
                else:
                    self.logger.warning("Unrecognized field {0} in fileorg argument".format(field))
        else:
            pass # TODO: default settings
        setattr(self, 'solver', {})
        if 'solver' in kwargs:
            for field, value in kwargs['solver'].items():
                if field in self._key_fields['solver']:
                    self.solver[field] = value
                else:
                    self.logger.warning("Unrecognized field {0} in solver argument".format(field))
        else:
            pass # TODO: default settings

class NdiayeAPLC(LyotCoronagraph): # Image-constrained APLC following N'Diaye et al. (2015, 2016)
    _design_fields = dict([('Pupil', ['N', 'pm', 'sm', 'ss', 'tel diam']),\
                           ('FPM', ['M', 'rad']),\
                           ('LS', ['id', 'od', 'ovszfac', 'altol']),\
                           ('Image', ['c', 'c1', 'c2', 'iwa', 'owa', 'owa2', 'bw', 'Nlam'])])
    _eval_fields =   dict([('Pupil', ['N', 'pm', 'sm', 'ss', 'tel diam']),\
                           ('FPM', ['M', 'rad']),\
                           ('LS', ['id', 'od', 'ovszfac', 'altol']),\
                           ('Image', ['c', 'c1', 'c2', 'iwa', 'owa', 'owa2', 'bw', 'Nlam']),\
                           ('Target', []), ('Aber', []), ('WFSC', [])])

    def __init__(self, **kwargs):
        super(NdiayeAPLC, self).__init__(**kwargs)

        setattr(self, 'design', {})
        for keycat in self._design_fields:
            self.design[keycat] = {}

        if 'design' in kwargs:
            design_dict = kwargs.get('design')
            for keycat, param_dict in design_dict.items():
                if keycat in self._design_fields:
                    for param, value in param_dict.items():
                        if param in self._design_fields[keycat]:
                            self.design[keycat][param] = value
                        else:
                            self.logger.warning("Unrecognized parameter \"{0}\" under category \"{1}\" in design initialization argument".format(param, keycat))
                else:
                    self.logger.warning("Unrecognized key category \"{0}\" in design initialization argument".format(keycat))
                    self.design[keycat] = None
        else:
            pass # TODO: default values
        
        self.ampl_fname = "ndiaye_aplc.mod"

    def write_ampl(self):
        self.logger.info("Writing the AMPL program")
    def read_solution(self):
        self.logger.info("Reading in the apodizer solution and parse the optimizer log")
    def create_eval_model(self):
        self.logger.info("Defining a python/Poppy model to evaluate the solution")
    def eval_solution(self):
        self.logger.info("Evaluating the design throughput, PSF FWHM, etc., and writing summary")

class HalfplaneAPLC(NdiayeAPLC): # N'Diaye APLC subclass for the half-plane symmetry case
    def __init__(self, **kwargs):
        super(HalfplaneAPLC, self).__init__(**kwargs)
        self.ampl_fname = "halfplane_aplc.mod"

    def write_ampl(self): 
        self.logger.info("Writing the AMPL program for this APLC with halfplane symmetry")
        self.abs_ampl_fname = os.path.join(self.fileorg['ampl src dir'], self.ampl_fname)
        if os.path.exists(self.abs_ampl_fname):
            self.logger.warning("Overwriting the existing copy of {0}".format(self.abs_ampl_fname))
        mod_fobj = open(self.abs_ampl_fname, "w")

        header = """\
        # AMPL program to optimize an APLC with half-plane symmetry
        # Created with {0} at {1}
        load amplgsl.dll;
        param pi:= 4*atan(1);
        """.format(os.path.basename(__file__), datetime.datetime.now())

        mod_fobj.write( textwrap.dedent(header) )

        mod_fobj.close()
        self.logger.info("Wrote %s"%self.abs_ampl_fname)
