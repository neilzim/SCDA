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
import numpy as np
import pdb

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
    _key_fields = dict([('fileorg', ['ampl src dir', 'tel ap dir', 'FPM dir', 'LS dir', 'sol dir', 'eval dir']),\
                        ('solver', ['constr', 'method', 'presolve'])])
    _solver_menu = dict([('constr',['lin', 'quad']), ('method', ['bar', 'barhom', 'dualsimp']), ('presolve',[True, False])])
    _aperture_menu = dict([('pm', ['hex1', 'hex2', 'hex3', 'key24', 'pie12', 'pie8', 'irisao']),\
                           ('ss', ['y60','y60off','x','cross','t','y90']),\
                           ('sst', ['025','100']),\
                           ('so', [True, False])])

    def __init__(self, **kwargs):
        # Only set fileorg and solver attributes in this constructor,
        # since design and eval parameter checking is design-specific. 
        self.logger = logging.getLogger('scda.logger')

        #////////////////////////////////////////////////////////////////////////////////////////////////////
        #   The fileorg attribute holds the locations of telescope apertures,
        #   intermediate masks, and co-eval AMPL programs.
        #   File names specific to a coronagraph design (AMPL program, apodizer solution, optimizer log, etc)
        #   are set downstream as object attributes.
        #////////////////////////////////////////////////////////////////////////////////////////////////////
        setattr(self, 'fileorg', {})
        if 'fileorg' in kwargs:
            for dirkey, location in kwargs['fileorg'].items():
                if dirkey in self._key_fields['fileorg']:
                    self.fileorg[dirkey] = location
                    if location is not None and not os.path.exists(location):
                        self.logger.warning("Warning: The specified location for \"{0}\", \"{1}\" does not yet exist".format(dirkey, location))
                else:
                    self.logger.warning("Warning: Unrecognized field {0} in fileorg argument".format(dirkey))
        # Handle missing values
        if 'ampl src dir' not in self.fileorg or self.fileorg['ampl src dir'] is None: # Set ampl src dir to cwd if unspecified
            self.fileorg['ampl src dir'] = os.getcwd()
        for dirkey in self._key_fields['fileorg']: # Revert missing locations to ampl src dir
            if dirkey not in self.fileorg or self.fileorg[dirkey] is None:
                self.fileorg[dirkey] = self.fileorg['ampl src dir']
                
        #////////////////////////////////////////////////////////////////////////////////////////////////////
        #   The solver attribute holds the options handed from AMPL to Gurobi, 
        #   and determines how the field constraints are mathematically expressed.
        #////////////////////////////////////////////////////////////////////////////////////////////////////
        setattr(self, 'solver', {})
        if 'solver' in kwargs:
            for field, value in kwargs['solver'].items():
                if field in self._key_fields['solver']:
                    if value in self._solver_menu[field]:
                        self.solver[field] = value
                    else:
                        self.logger.warning("Warning: Unrecognized solver option \"{0}\" in field \"{1}\", reverting to default".format(value, field))
                else:
                    self.logger.warning("Warning: Unrecognized field {0} in solver argument".format(field))
        # Handle missing values
        if 'constr' not in self.solver or self.solver['constr'] is None: self.solver['constr'] = 'lin'
        if 'method' not in self.solver or self.solver['method'] is None: self.solver['method'] = 'bar'
        if 'presolve' not in self.solver or self.solver['presolve'] is None: self.solver['presolve'] = True

class NdiayeAPLC(LyotCoronagraph): # Image-constrained APLC following N'Diaye et al. (2015, 2016)
    _design_fields = { 'Pupil': {'N':(int, 1000), 'pm':(str, 'hex1'), 'so':(bool, True), 'ss':(str, 'x'), 'sst':(str, '100'), 'tel diam':(float, 12.)}, \
                       'FPM':   {'M':(int, 50), 'rad':(float, 4.)}, \
                       'LS':    {'id':(float, 0.2), 'od':(float, 0.9), 'ovszfac':(float, None), 'altol':(float, None)}, \
                       'Image': {'c':(float, 10.), 'c1':(float, None), 'c2':(float, None), 'iwa':(float, 4.), \
                                 'owa':(float, 10.), 'owa2':(float, None), 'fpres':(int,2), 'bw':(float, 0.1), 'Nlam':(int, 3)} }
    _eval_fields =   { 'Pupil': _design_fields['Pupil'], 'FPM': _design_fields['FPM'], \
                       'LS': _design_fields['LS'], 'Image': _design_fields['Image'], \
                       'Target': {}, 'Aber': {}, 'WFSC': {} }

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
                            if value is not None:
                                if isinstance(value, self._design_fields[keycat][param][0]):
                                    self.design[keycat][param] = value
                                else:
                                    warnstr = ("Warning: Invalid {0} for parameter \"{1}\" under category \"{2}\" " + \
                                               "design initialization argument, expecting a {3}").format(type(value), param, keycat, self._design_fields[keycat][param][0]) 
                                    self.logger.warning(warnstr)
                        else:
                            self.logger.warning("Warning: Unrecognized parameter \"{0}\" under category \"{1}\" in design initialization argument".format(param, keycat))
                else:
                    self.logger.warning("Warning: Unrecognized key category \"{0}\" in design initialization argument".format(keycat))
                    self.design[keycat] = None
        for keycat in self._design_fields: # Fill in default values where appropriate
            for param in self._design_fields[keycat]:
                if param not in self.design[keycat] or (self.design[keycat][param] is None and \
                                                        self._design_fields[keycat][param][1] is not None):
                    self.design[keycat][param] = self._design_fields[keycat][param][1]
     
        self.amplname_coron = "APLC" 
        if self.design['Pupil']['so'] == True:
            self.amplname_pupil = "{}{}t{}so1D{:04}".format(self.design['Pupil']['pm'], self.design['Pupil']['ss'], self.design['Pupil']['sst'], self.design['Pupil']['N'])
        else:                                                                                                                                      
            self.amplname_pupil = "{}{}t{}so0D{:04}".format(self.design['Pupil']['pm'], self.design['Pupil']['ss'], self.design['Pupil']['sst'], self.design['Pupil']['N'])

        self.amplname_fpm = "FPM{:02}M{:03}".format(int(round(100*self.design['FPM']['rad'])), self.design['FPM']['M'])
        if self.design['LS']['altol'] is None and self.design['LS']['ovszfac'] is None:
            self.amplname_ls = "LS{:02}D{:02}ovsz{:02}tol{:02}".format(int(round(100*self.design['LS']['id'])), \
                               int(round(100*self.design['LS']['od'])), 0, 0)
        elif self.design['LS']['altol'] is None and self.design['LS']['ovszfac'] is not None:
            self.amplname_ls = "LS{:02}D{:02}ovsz{:02}tol{:02}".format(int(round(100*self.design['LS']['id'])), \
                               int(round(100*self.design['LS']['od'])), int(round(100*self.design['LS']['ovszfac'])), 0)
        elif self.design['LS']['altol'] is not None and self.design['LS']['ovszfac'] is None:
            self.amplname_ls = "LS{:02}D{:02}ovsz{:02}tol{:02}".format(int(round(100*self.design['LS']['id'])), \
                               int(round(100*self.design['LS']['od'])), 0, int(round(100*self.design['LS']['altol'])))
        else:
            self.amplname_ls = "LS{:02}D{:02}ovsz{:02}tol{:02}".format(int(round(100*self.design['LS']['id'])), \
                               int(round(100*self.design['LS']['od'])), int(round(100*self.design['LS']['ovszfac'])), \
                               int(round(100*self.design['LS']['altol'])))

        self.amplname_image = "Img{:03}C{:02}WA{:03}BW{:02}Nlam{:02}fpres{:1}".format(int(round(10*self.design['Image']['c'])), \
                               int(round(10*self.design['Image']['iwa'])), int(round(10*self.design['Image']['owa'])), \
                               int(round(100*self.design['Image']['bw'])), self.design['Image']['Nlam'], self.design['Image']['fpres'])
        if self.design['Image']['c1'] is not None:
            self.amplname_image += "C1{:03}".format(self.design['Image']['c1'])
        if self.design['Image']['c2'] is not None:
            self.amplname_image += "C2{:03}".format(self.design['Image']['c2'])
        if self.design['Image']['owa2'] is not None:
            self.amplname_image += "OWA2{:03}".format(self.design['Image']['owa2'])

        if self.solver['presolve']:
            self.amplname_solver = "{}{}pre1".format(self.solver['constr'], self.solver['method'])
        else:
            self.amplname_solver = "{}{}pre0".format(self.solver['constr'], self.solver['method'])

        self.ampl_src_fname = self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                              self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".mod" 

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
        self.amplname_coron = "hpAPLC" 
        self.ampl_src_fname = self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                              self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".mod" 

    def write_ampl(self): 
        self.logger.info("Writing the AMPL program for this APLC with halfplane symmetry")
        self.abs_ampl_src_fname = os.path.join(self.fileorg['ampl src dir'], self.ampl_src_fname)
        if not os.path.exists(self.fileorg['ampl src dir']):
           os.mkdir(self.fileorg['ampl src dir'])
        if os.path.exists(self.abs_ampl_src_fname):
            self.logger.warning("Warning: Overwriting the existing copy of {0}".format(self.abs_ampl_src_fname))
        mod_fobj = open(self.abs_ampl_src_fname, "w")

        header = """\
        # AMPL program to optimize an APLC with half-plane symmetry
        # Created with {0} at {1}
        load amplgsl.dll;
        param pi:= 4*atan(1);
        """.format(os.path.basename(__file__), datetime.datetime.now())

        mod_fobj.write( textwrap.dedent(header) )

        mod_fobj.close()
        self.logger.info("Wrote %s"%self.abs_ampl_src_fname)
