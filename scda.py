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
import getpass
import socket
from collections import defaultdict, OrderedDict
import itertools
import pprint

def configure_log(log_fname=None):
    logger = logging.getLogger("scda.logger")
    if not len(logger.handlers):
        logger.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout) # console dump
        ch.setLevel(logging.INFO)
        logger.addHandler(ch)
        if sys.version_info[0] >= 3:
            ch.setFormatter(WrappedFixedIndentingLog(indent=8, width=120))

        if log_fname is not None: # optional file log in addition to the console dump
            fh = logging.FileHandler(log_fname, mode="w")
            fh.setLevel(logging.DEBUG)
            logger.addHandler(fh)

class WrappedFixedIndentingLog(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None, style='%', width=70, indent=4):
        super(WrappedFixedIndentingLog, self).__init__(fmt=fmt, datefmt=datefmt)
        self.wrapper = textwrap.TextWrapper(width=width)
        #self.wrapper = textwrap.TextWrapper(width=width, subsequent_indent=' '*indent)
    def format(self, record):
        return self.wrapper.fill(super().format(record))

class DesignParamSurvey(object):
    def __init__(self, coron_class, survey_config, **kwargs):
        self.logger = logging.getLogger('scda.logger')
        _param_menu = coron_class._design_fields.copy()
        _file_fields = coron_class._file_fields.copy()
        setattr(self, 'survey_config', {})
        for keycat, param_dict in survey_config.items():
            self.survey_config[keycat] = {}
            if keycat in _param_menu:
                for param, values in param_dict.items():
                    if param in _param_menu[keycat]:
                        if values is not None:
                            if hasattr(values, '__iter__'): #check the type of all items
                                if all(isinstance(value, _param_menu[keycat][param][0]) for value in values):
                                    self.survey_config[keycat][param] = values
                                    #self.survey_config[keycat][param] = tuple(values)
                                else:
                                    warnstr = ("Warning: Invalid type found in survey set {0} for parameter {1} under category \"{2}\" " + \
                                               "design initialization argument, expecting {3}").format(value, param, keycat, _param_menu[keycat][param][0]) 
                                    self.logger.warning(warnstr)
                            else:
                                if isinstance(values, _param_menu[keycat][param][0]):
                                    self.survey_config[keycat][param] = values
                                else:
                                    warnstr = ("Warning: Invalid {0} for parameter \"{1}\" under category \"{2}\" " + \
                                               "design initialization argument, expecting a {3}").format(type(value), param, keycat, _param_menu[keycat][param][0]) 
                                    self.logger.warning(warnstr)
                    else:
                        self.logger.warning("Warning: Unrecognized parameter \"{0}\" under category \"{1}\" in design initialization argument".format(param, keycat))
            else:
                self.logger.warning("Warning: Unrecognized key category \"{0}\" in design initialization argument".format(keycat))
                self.survey_config[keycat] = None
        varied_param_flat = []
        varied_param_index = []
        fixed_param_flat = []
        fixed_param_index = []
        for keycat in _param_menu: # Fill in default values where appropriate
            if keycat not in self.survey_config:
                self.survey_config[keycat] = {}
            for param in _param_menu[keycat]:
                if param not in self.survey_config[keycat] or (self.survey_config[keycat][param] is None and \
                                                               _param_menu[keycat][param][1] is not None):
                    self.survey_config[keycat][param] = _param_menu[keycat][param][1] # default value
                elif param in self.survey_config[keycat] and self.survey_config[keycat][param] is not None and \
                not hasattr(self.survey_config[keycat][param], '__iter__'):
                    fixed_param_flat.append(self.survey_config[keycat][param])
                    fixed_param_index.append((keycat, param))
                elif hasattr(self.survey_config[keycat][param], '__iter__'):
                    varied_param_flat.append(self.survey_config[keycat][param])
                    varied_param_index.append((keycat, param))
     
        varied_param_combos = []
        for combo in itertools.product(*varied_param_flat):
            varied_param_combos.append(combo)
        self.varied_param_combos = tuple(varied_param_combos)
        self.varied_param_index = tuple(varied_param_index)
        self.fixed_param_vals = tuple(fixed_param_flat)
        self.fixed_param_index = tuple(fixed_param_index)

        #////////////////////////////////////////////////////////////////////////////////////////////////////
        #   The fileorg attribute holds the locations of telescope apertures,
        #   intermediate masks, co-eval AMPL programs, solutilons, logs, etc.
        #////////////////////////////////////////////////////////////////////////////////////////////////////
        setattr(self, 'fileorg', {})
        if 'fileorg' in kwargs:
            for namekey, location in kwargs['fileorg'].items():
                if namekey in _file_fields['fileorg']:
                    if location is not None:
                        if namekey.endswith('dir'):
                            self.fileorg[namekey] = os.path.abspath(os.path.expanduser(location)) # Convert all directory names to absolute paths
                            if not os.path.exists(self.fileorg[namekey]):
                                self.logger.warning("Warning: The specified location of '{0}', \"{1}\" does not exist".format(namekey, self.fileorg[namekey]))
                        else:
                            self.fileorg[namekey] = location
                    else:
                        self.fileorg[namekey] = None
                else:
                    self.logger.warning("Warning: Unrecognized field {0} in fileorg argument".format(dirkey))
        # Handle missing directory values
        if 'work dir' not in self.fileorg or self.fileorg['work dir'] is None:
            self.fileorg['work dir'] = os.getcwd()
        for namekey in _file_fields['fileorg']: # Set other missing directory locations to 'work dir'
            if namekey.endswith('dir') and ( namekey not in self.fileorg or self.fileorg[namekey] is None ):
                self.fileorg[namekey] = self.fileorg['work dir']
      
        # In most cases we don't expect to directly specify file names for apertures, FPM, or LS files
        # for the SCDA parameter survey. However, it easy enough to make this option available.
        # If the location of the optimizer input file is not known, 
        # look for it in the directory corresponding to its specific category
        if 'TelAp fname' in self.fileorg and self.fileorg['TelAp fname'] is not None and \
        not os.path.exists(self.fileorg['TelAp fname']) and os.path.exists(self.fileorg['TelAp dir']) and \
        os.path.dirname(self.fileorg['TelAp fname']) == '':
            try_fname = os.path.join(self.fileorg['TelAp dir'], self.fileorg['TelAp fname']) 
            if os.path.exists(try_fname):
                self.fileorg['TelAp fname'] = try_fname
            else:
                self.logger.warning("Warning: Could not find the specified telescope aperture file \"{0}\" in {1}".format(self.fileorg['TelAp fname'],
                                    self.fileorg['TelAp dir']))
        if 'FPM fname' in self.fileorg and self.fileorg['FPM fname'] is not None and \
        not os.path.exists(self.fileorg['FPM fname']) and os.path.exists(self.fileorg['FPM dir']) and \
        os.path.dirname(self.fileorg['FPM fname']) == '':
            try_fname = os.path.join(self.fileorg['FPM dir'], self.fileorg['FPM fname']) 
            if os.path.exists(try_fname):
                self.fileorg['FPM fname'] = try_fname
            else:
                self.logger.warning("Warning: Could not find the specified FPM file \"{0}\" in {1}".format(self.fileorg['FPM fname'],
                                    self.fileorg['FPM dir']))
        if 'LS fname' in self.fileorg and self.fileorg['LS fname'] is not None and \
        not os.path.exists(self.fileorg['LS fname']) and os.path.exists(self.fileorg['LS dir']) and \
        os.path.dirname(self.fileorg['LS fname']) == '':
            try_fname = os.path.join(self.fileorg['LS dir'], self.fileorg['LS fname']) 
            if os.path.exists(try_fname):
                self.fileorg['LS fname'] = try_fname
            else:
                self.logger.warning("Warning: Could not find the specified LS file \"{0}\" in {1}".format(self.fileorg['LS fname'],
                                    self.fileorg['LS dir']))
        #////////////////////////////////////////////////////////////////////////////////////////////////////
        #   The solver attribute holds the options handed from AMPL to Gurobi, 
        #   and determines how the field constraints are mathematically expressed.
        #////////////////////////////////////////////////////////////////////////////////////////////////////
        setattr(self, 'solver', {})
        if 'solver' in kwargs:
            for field, value in kwargs['solver'].items():
                if field in self._file_fields['solver']:
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
        if 'Nthreads' not in self.solver or self.solver['Nthreads'] is None: self.solver['Nthreads'] = None
         
        setattr(self, 'coron_list', [])
        design = {}
        for keycat in _param_menu:
            design[keycat] = {}
        for (fixed_keycat, fixed_parname), fixed_val in zip(self.fixed_param_index, self.fixed_param_vals):
            design[fixed_keycat][fixed_parname] = fixed_val
        self.coron_list = []
        for param_combo in self.varied_param_combos:
            for (varied_keycat, varied_parname), current_val in zip(self.varied_param_index, param_combo):
                design[varied_keycat][varied_parname] = current_val
            self.coron_list.append( coron_class(design=design, fileorg=self.fileorg, solver=self.solver) )

    def write_ampl(self, overwrite=False, override_infile_status=False):
        for coron in self.coron_list:
            coron.write_ampl(overwrite, override_infile_status)

class LyotCoronagraph(object): # Lyot coronagraph base class
    _file_fields = { 'fileorg': ['work dir', 'ampl src dir', 'TelAp dir', 'FPM dir', 'LS dir', 'sol dir', 'eval dir',
                                 'ampl src fname', 'TelAp fname', 'FPM fname', 'LS fname', 'sol fname'],
                     'solver': ['constr', 'method', 'presolve', 'Nthreads'] }

    _solver_menu = { 'constr': ['lin', 'quad'], 'solver': ['LOQO', 'gurobi', 'gurobix'], 
                     'method': ['bar', 'barhom', 'dualsimp'],
                     'presolve': [True, False], 'Nthreads': [None]+range(1,33) }

    _aperture_menu = { 'pm': ['hex1', 'hex2', 'hex3', 'hex4', 'key24', 'pie12', 'pie8', 'irisao'],
                       'ss': ['y60','y60off','x','cross','t','y90'],
                       'sst': ['025','100'],
                       'so': [True, False] }

    def __init__(self, verbose=False, **kwargs):
        # Only set fileorg and solver attributes in this constructor,
        # since design and eval parameter checking is design-specific. 
        self.logger = logging.getLogger('scda.logger')

        #////////////////////////////////////////////////////////////////////////////////////////////////////
        #   The fileorg attribute holds the locations of telescope apertures,
        #   intermediate masks, co-eval AMPL programs, solutilons, logs, etc.
        #////////////////////////////////////////////////////////////////////////////////////////////////////
        setattr(self, 'fileorg', {})
        if 'fileorg' in kwargs:
            for namekey, location in kwargs['fileorg'].items():
                if namekey in self._file_fields['fileorg']:
                    if location is not None:
                        if namekey.endswith('dir'):
                            self.fileorg[namekey] = os.path.abspath(os.path.expanduser(location)) # Convert all directory names to absolute paths
                            if not os.path.exists(self.fileorg[namekey]):
                                self.logger.warning("Warning: The specified location of '{0}', \"{1}\" does not exist".format(namekey, self.fileorg[namekey]))
                        else:
                            self.fileorg[namekey] = location
                    else:
                        self.fileorg[namekey] = None
                else:
                    self.logger.warning("Warning: Unrecognized field {0} in fileorg argument".format(dirkey))
        # Handle missing directory values
        if 'work dir' not in self.fileorg or self.fileorg['work dir'] is None:
            self.fileorg['work dir'] = os.getcwd()
        for namekey in self._file_fields['fileorg']: # Set other missing directory locations to 'work dir'
            if namekey.endswith('dir') and ( namekey not in self.fileorg or self.fileorg[namekey] is None ):
                self.fileorg[namekey] = self.fileorg['work dir']
       
        # If the location of the optimizer input file is not known, 
        # look for it in the directory corresponding to its specific category
        if 'TelAp fname' in self.fileorg and self.fileorg['TelAp fname'] is not None and \
        not os.path.exists(self.fileorg['TelAp fname']) and os.path.exists(self.fileorg['TelAp dir']) and \
        os.path.dirname(self.fileorg['TelAp fname']) == '':
            try_fname = os.path.join(self.fileorg['TelAp dir'], self.fileorg['TelAp fname']) 
            if os.path.exists(try_fname):
                self.fileorg['TelAp fname'] = try_fname
            else:
                self.logger.warning("Warning: Could not find the specified telescope aperture file \"{0}\" in {1}".format(self.fileorg['TelAp fname'], \
                                    self.fileorg['TelAp dir']))
        if 'FPM fname' in self.fileorg and self.fileorg['FPM fname'] is not None and \
        not os.path.exists(self.fileorg['FPM fname']) and os.path.exists(self.fileorg['FPM dir']) and \
        os.path.dirname(self.fileorg['FPM fname']) == '':
            try_fname = os.path.join(self.fileorg['FPM dir'], self.fileorg['FPM fname']) 
            if os.path.exists(try_fname):
                self.fileorg['FPM fname'] = try_fname
            else:
                self.logger.warning("Warning: Could not find the specified FPM file \"{0}\" in {1}".format(self.fileorg['FPM fname'], \
                                    self.fileorg['FPM dir']))
        if 'LS fname' in self.fileorg and self.fileorg['LS fname'] is not None and \
        not os.path.exists(self.fileorg['LS fname']) and os.path.exists(self.fileorg['LS dir']) and \
        os.path.dirname(self.fileorg['LS fname']) == '':
            try_fname = os.path.join(self.fileorg['LS dir'], self.fileorg['LS fname']) 
            if os.path.exists(try_fname):
                self.fileorg['LS fname'] = try_fname
            else:
                self.logger.warning("Warning: Could not find the specified LS file \"{0}\" in {1}".format(self.fileorg['LS fname'], \
                                    self.fileorg['LS dir']))

        # If the specified ampl source filename is a simple name with no directory, append it to the ampl source directory.
        if 'ampl src fname' in self.fileorg and self.fileorg['ampl src fname'] is not None and \
        not os.path.exists(self.fileorg['ampl src fname']) and os.path.dirname(self.fileorg['ampl src fname']) == '':
            self.fileorg['ampl src fname'] = os.path.join(self.fileorg['ampl src dir'], self.fileorg['ampl src fname'])
                
        #////////////////////////////////////////////////////////////////////////////////////////////////////
        #   The solver attribute holds the options handed from AMPL to Gurobi, 
        #   and determines how the field constraints are mathematically expressed.
        #////////////////////////////////////////////////////////////////////////////////////////////////////
        setattr(self, 'solver', {})
        if 'solver' in kwargs:
            for field, value in kwargs['solver'].items():
                if field in self._file_fields['solver']:
                    if value in self._solver_menu[field]:
                        self.solver[field] = value
                    else:
                        self.logger.warning("Warning: Unrecognized solver option \"{0}\" in field \"{1}\", reverting to default".format(value, field))
                else:
                    self.logger.warning("Warning: Unrecognized field {0} in solver argument".format(field))
        # Handle missing values
        if 'constr' not in self.solver or self.solver['constr'] is None: self.solver['constr'] = 'lin'
        if 'solver' not in self.solver or self.solver['solver'] is None: self.solver['solver'] = 'gurobi'
        if 'method' not in self.solver or self.solver['method'] is None: self.solver['method'] = 'bar'
        if 'presolve' not in self.solver or self.solver['presolve'] is None: self.solver['presolve'] = True
        if 'Nthreads' not in self.solver or self.solver['Nthreads'] is None: self.solver['Nthreads'] = None

        setattr(self, 'ampl_infile_status', None)
        if not issubclass(self.__class__, LyotCoronagraph):
            self.check_ampl_input_files()

    def check_ampl_input_files(self):
        status = True
        checklist = ['TelAp fname', 'FPM fname', 'LS fname']
        for fname in checklist:
            if not os.path.exists(self.fileorg[fname]):
                status = False
                break
        self.ampl_infile_status = status

class NdiayeAPLC(LyotCoronagraph): # Image-constrained APLC following N'Diaye et al. (2015, 2016)
    _design_fields = OrderedDict([ ( 'Pupil', OrderedDict([('N',(int, 1000)), ('pm',(str, 'hex1')), ('ss',(str, 'x')), 
                                                          ('sst',(str, '100')), ('so',(bool, True))]) ),
                                   ( 'FPM', OrderedDict([('rad',(float, 4.)), ('M',(int, 50))]) ),
                                   ( 'LS', OrderedDict([('id',(int, 20)), ('od',(int, 90)), ('ovsz',(int, 0)), ('altol',(int, None))]) ),
                                   ( 'Image', OrderedDict([('c',(float, 10.)), ('iwa',(float, 4.)), ('owa',(float, 10.)),
                                                         ('bw',(float, 0.1)), ('Nlam',(int, 3)), ('fpres',(int,2)),
                                                         ('oca',(float, 10.)), ('ci',(float, None)), ('co',(float, None))]) ) ])
    _eval_fields =   { 'Pupil': _design_fields['Pupil'], 'FPM': _design_fields['FPM'], \
                       'LS': _design_fields['LS'], 'Image': _design_fields['Image'], \
                       'Tel': {'TelAp diam':(float, 12.)}, 'Target': {}, 'Aber': {}, 'WFSC': {} }

    def __init__(self, verbose=False, **kwargs):
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
            if keycat not in self.design:
                self.design[keycat] = {}
            for param in self._design_fields[keycat]:
                if param not in self.design[keycat] or (self.design[keycat][param] is None and \
                                                        self._design_fields[keycat][param][1] is not None):
                    self.design[keycat][param] = self._design_fields[keycat][param][1]
        # Finally, set a private attribute for the number of image plane samples between the center and the OCA
        self.design['Image']['_Nimg'] = int( np.ceil( self.design['Image']['fpres']*self.design['Image']['oca']/(1. - self.design['Image']['bw']/2) ) )
        if verbose: # Print summary of the set parameters
            self.logger.info("Design parameters: {}".format(self.design))
            self.logger.info("Optimization and solver parameters: {}".format(self.solver))
            self.logger.info("File organization parameters: {}".format(self.fileorg))
     
        self.amplname_coron = "APLC_full"
        if self.design['Pupil']['so'] == True:
            self.amplname_pupil = "{0:s}{1:s}t{2:s}so1_N{3:04d}".format(self.design['Pupil']['pm'], self.design['Pupil']['ss'], self.design['Pupil']['sst'], self.design['Pupil']['N'])
        else:                                                                                                                                      
            self.amplname_pupil = "{0:s}{1:s}t{2:s}so0_N{3:04d}".format(self.design['Pupil']['pm'], self.design['Pupil']['ss'], self.design['Pupil']['sst'], self.design['Pupil']['N'])

        self.amplname_fpm = "FPM{:02}M{:03}".format(int(round(100*self.design['FPM']['rad'])), self.design['FPM']['M'])
        if self.design['LS']['altol'] is None and self.design['LS']['ovsz'] is None:
            self.amplname_ls = "LS{:02}D{:02}ovsz{:02}tol{:02}".format(self.design['LS']['id'], \
                               self.design['LS']['od'], 0, 0)
        elif self.design['LS']['altol'] is None and self.design['LS']['ovsz'] is not None:
            self.amplname_ls = "LS{:02}D{:02}ovsz{:02}tol{:02}".format(self.design['LS']['id'], \
                               self.design['LS']['od'], self.design['LS']['ovsz'], 0)
        elif self.design['LS']['altol'] is not None and self.design['LS']['ovsz'] is None:
            self.amplname_ls = "LS{:02}D{:02}ovsz{:02}tol{:02}".format(self.design['LS']['id'], \
                               self.design['LS']['od'], 0, self.design['LS']['altol'])
        else:
            self.amplname_ls = "LS{:02}D{:02}ovsz{:02}tol{:02}".format(self.design['LS']['id'], \
                               self.design['LS']['od'], self.design['LS']['ovsz'], self.design['LS']['altol'])

        self.amplname_image = "Img{:03}C_{:02}WA{:03}CA{:03}_BW{:02}Nlam{:02}fpres{:1}".format(int(round(10*self.design['Image']['c'])), \
                               int(round(10*self.design['Image']['iwa'])), int(round(10*self.design['Image']['owa'])), \
                               int(round(10*self.design['Image']['oca'])), \
                               int(round(100*self.design['Image']['bw'])), self.design['Image']['Nlam'], self.design['Image']['fpres'])
        if self.design['Image']['ci'] is not None:
            self.amplname_image += "Cin{:03}".format(self.design['Image']['ci'])
        if self.design['Image']['co'] is not None:
            self.amplname_image += "Cout{:03}".format(self.design['Image']['co'])

        if self.solver['presolve']:
            self.amplname_solver = "{}{}pre1".format(self.solver['constr'], self.solver['method'])
        else:
            self.amplname_solver = "{}{}pre0".format(self.solver['constr'], self.solver['method'])
        if self.solver['Nthreads'] is not None:
            self.amplname_solver += "thr{:02d}".format(self.solver['Nthreads'])

        if not issubclass(self.__class__, NdiayeAPLC): # Only set these file names if this is a full-plane APLC.
            if 'ampl src fname' not in self.fileorg or self.fileorg['ampl src fname'] is None:
                ampl_src_fname_tail = self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                                      self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".mod"
                self.fileorg['ampl src fname'] = os.path.join(self.fileorg['ampl src dir'], ampl_src_fname_tail)

            if 'sol fname' not in self.fileorg or self.fileorg['sol fname'] is None:
                self.fileorg['sol fname'] = self.fileorg['ampl src fname'][:-4] + "solAp.dat"

            if 'TelAp fname' not in self.fileorg or self.fileorg['TelAp fname'] is None:
                self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_full_" + self.amplname_pupil + ".dat") )

            if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
                self.fileorg['FPM fname'] = os.path.join( self.fileorg['FPM dir'], "FPM_occspot_M{:03}.dat".format(self.design['FPM']['M']) )

            if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_full_" + self.amplname_pupil + \
                                                         "_{0:02d}D{1:02d}ovsz{2:02d}.dat".format(self.design['LS']['id'], \
                                                         self.design['LS']['od'], self.design['LS']['ovsz'])) )
            self.check_ampl_input_files()
                                                              
    def write_ampl(self, overwrite=False):
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
        self.amplname_coron = "APLC_half"
        if 'ampl src fname' not in self.fileorg or self.fileorg['ampl src fname'] is None:
            ampl_src_fname_tail = self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                                  self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".mod"
            self.fileorg['ampl src fname'] = os.path.join(self.fileorg['ampl src dir'], ampl_src_fname_tail)

        if 'sol fname' not in self.fileorg or self.fileorg['sol fname'] is None:
            self.fileorg['sol fname'] = self.fileorg['ampl src fname'][:-4] + "solAp.dat"

        if 'TelAp fname' not in self.fileorg or self.fileorg['TelAp fname'] is None:
            self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_half_" + self.amplname_pupil + ".dat") )

        if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
            self.fileorg['FPM fname'] = os.path.join( self.fileorg['FPM dir'], "FPM_occspot_M{:03}.dat".format(self.design['FPM']['M']) )

        if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
            self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + self.amplname_pupil + \
                                                     "_{0:02d}D{1:02d}ovsz{2:02d}.dat".format(self.design['LS']['id'], \
                                                     self.design['LS']['od'], self.design['LS']['ovsz'])) )
        self.check_ampl_input_files()
    def write_ampl(self, overwrite=False, override_infile_status=False, ampl_src_fname=None):
        if self.ampl_infile_status is False and not override_infile_status:
            self.logger.error("Error: the most recent input file check for this design configuration failed.")
            self.logger.error("The override_infile_status switch is off, so write_ampl() will now abort.")
            self.logger.error("See previous warnings in the log to see what file was missing during the initialization")
            return
        if ampl_src_fname is not None:
            if os.path.dirname(ampl_src_fname) == '' and self.fileorg['ampl src dir'] is not None:
                self.fileorg['ampl src fname'] = os.path.join(self.fileorg['ampl src dir'], ampl_src_fname)
            else:
                self.fileorg['ampl src fname'] = os.path.abspath(ampl_src_fname)    
                self.fileorg['ampl src dir'] = os.path.dirname(self.fileorg['ampl src fname'])
        if os.path.exists(self.fileorg['ampl src fname']):
            if overwrite == True:
                self.logger.warning("Warning: Overwriting the existing copy of {0}".format(self.fileorg['ampl src fname']))
            else:
                self.logger.warning("Error: {0} already exists and overwrite switch is off, so write_ampl() will now abort".format(self.fileorg['ampl src fname']))
                return
        elif not os.path.exists(self.fileorg['ampl src dir']):
            os.mkdir(self.fileorg['ampl src dir'])
            self.logger.info("Created new AMPL source code directory, {0:s}".format(self.fileorg['ampl src dir']))
#        self.logger.info("Writing the AMPL program for the specified half-plane APLC")
        mod_fobj = open(self.fileorg['ampl src fname'], "w")

        header = """\
        # AMPL program to optimize a half-plane symmetric APLC
        # Created by {0:s} with {1:s} on {2:s} at {3:s}
        load amplgsl.dll;
        """.format(getpass.getuser(), os.path.basename(__file__), socket.gethostname(), datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))

        params = """
        #---------------------

        param pi:= 4*atan(1);

        #---------------------
        param c := {0:.2f};

        #---------------------
        param Rmask := {1:0.3f};
        param rho0 := {2:0.2f};
        param rho1 := {3:0.2f};
        
        #---------------------
        param N := {4:d};				# discretization parameter (pupil)
        param M := {5:d};				# discretization parameter (mask)
        
        param Nimg := {6:d};           # discretization parameter (image)
        param Fmax := {7:0.2f};        # NTZ: We paramaterize our image plane resolution by fpres = sampling rate at 
                                    #      the shortest wavelength. Then Nimg is an integer function of fpres, oca,
        #---------------------      #      and bw. This integer is not specified directly by the user, but computed "privately"
        param bw := {8:0.2f};           #      by the APLC class constructor.
        param lam0 := 1.;
        param dl := bw*lam0;
        param Nlam := {9:d};
        
        #---------------------
        param obs := 20;             # NTZ: We will eliminate this section of parameter definitions, 
        param spiders :=01;          #      since we determine the file names of the aperture and Lyot stop
        param lsobs := 20;           #      outside the program, as well as the file name of the apodizer solution.
        param lsspiders :=02;
        param gap :=01;
        param lsgap :=00;
        
        param OD :=0.98;
        
        #---------------------
        param Normterm := 1.00; 	# 0.380
        
        #---------------------
        param CoeffOverSizePup :=0.824*OD;
        """.format(self.design['Image']['c'], self.design['FPM']['rad'], self.design['Image']['iwa'], self.design['Image']['owa'], \
                   self.design['Pupil']['N'], self.design['FPM']['M'], self.design['Image']['_Nimg'], \
                   self.design['Image']['oca'], self.design['Image']['bw'], self.design['Image']['Nlam'])
        
        load_masks = """\
        #---------------------
        # Loading Pupil
        param PupilFile {{1..2*N,1..N}};
        
        read {{i in 1..2*N,j in 1..N}} PupilFile[i,j] < "{0:s}";
        close "{1:s}";
        
        # Loading FPM
        param MaskFile {{1..2*M,1..M}};
        
        read {{i in 1..2*M,j in 1..M}} MaskFile[i,j] < "{2:s}"; 
        close "{3:s}";
        
        # Loading Lyot stop
        param LyotFile {{1..2*N,1..N}};
        
        read {{i in 1..2*N,j in 1..N}} LyotFile[i,j] < "{4:s}";
        close "{5:s}";
        """.format(self.fileorg['TelAp fname'], self.fileorg['TelAp fname'], self.fileorg['FPM fname'], self.fileorg['FPM fname'], \
                   self.fileorg['LS fname'], self.fileorg['LS fname'])

        define_coords = """\
         
        #---------------------
        # steps in each plane
        param dx := 1/(2*N);
        
        param dy := dx;
        
        param dmx := 2*Rmask/(2*M);
        param dmy := dmx;
        
        param dxi := (Fmax/Nimg)*(1/CoeffOverSizePup);
        param deta := dxi;

        #---------------------
        # coordinate vectors in each plane
        set Xs := setof {i in -N+0.5..N-0.5 by 1} i*dx;
        set Ys := setof {j in 0.5..N-0.5 by 1} j*dy;
        
        set MXs := setof {i in -M+0.5..M-0.5 by 1} i*dmx;
        set MYs := setof {j in 0.5..M-0.5 by 1} j*dmy;
        
        set Ls := setof {l in 1..Nlam} lam0*(1+((l-1)/(Nlam-1)-0.5)*dl);
        """

        sets_and_arrays = """
        #---------------------
        set Pupil := setof {x in Xs, y in Ys: PupilFile[round(x/dx+0.5+N),round(y/dy+0.5)] != 0} (x,y);
        set Mask := setof {mx in MXs, my in MYs: MaskFile[round(mx/dmx+0.5+M),round(my/dmy+0.5)] != 0} (mx,my);
        set Lyot := setof {x in Xs, y in Ys: LyotFile[round(x/dx+0.5+N),round(y/dy+0.5)] != 0} (x,y);
        
        param TR := sum {i in 1..2*N, j in 1..N} PupilFile[i,j]*dx*dy; # Transmission of the Pupil. Used for calibration.
        param I00 := (sum {i in 1..2*N, j in 1..N} PupilFile[i,j]*LyotFile[i,j]*dx*dy)^2; # Peak intensity in the absence of coronagraph (Pupil and Lyot terms are squared but square exponiential is unnecessary because of the binary values).
        
        var A {x in Xs, y in Ys} >= 0, <= 1, := 0.5;
        
        #---------------------
        
        set Xis := setof {i in -Nimg..Nimg-1 by 1} i*dxi;
        set Etas := setof {j in 0..Nimg-1 by 1} j*deta;
        
        set DarkHole := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= rho0 && sqrt(xi^2+eta^2) <= rho1} (xi,eta); # Only for 360deg masks.
        #set PSFCore := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= 0 && sqrt(xi^2+eta^2) < rho0} (xi,eta); # Only for 360deg masks.
        #set InsideArea := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= 0 && sqrt(xi^2+eta^2) <= rho1} (xi,eta); # Only for 360deg masks.
        """

        field_propagation = """
        #---------------------
        var EBm_cx {x in Xs, my in MYs, lam in Ls} := 0.0;
        var EBm_real {mx in MXs, my in MYs, lam in Ls} := 0.0;
        var EBm_imag {mx in MXs, my in MYs, lam in Ls} := 0.0;
        
        subject to st_EBm_cx {x in Xs, my in MYs, lam in Ls}: EBm_cx[x,my,lam] = 2.*sum {y in Ys: (x,y) in Pupil} A[x,y]*PupilFile[round(x/dx+0.5+N),round(y/dy+0.5)]*cos(2*pi*y*my*(lam0/lam))*dy;
        subject to st_EBm_real {(mx,my) in Mask, lam in Ls}: EBm_real[mx,my,lam] = (lam0/lam)*sum {x in Xs} EBm_cx[x,my,lam]*cos(2*pi*x*mx*(lam0/lam))*dx;
        subject to st_EBm_imag {(mx,my) in Mask, lam in Ls}: EBm_imag[mx,my,lam] = (lam0/lam)*sum {x in Xs} -EBm_cx[x,my,lam]*sin(2*pi*x*mx*(lam0/lam))*dx;
        
        #---------------------
        var ECm1_Bmreal_cx {mx in MXs, y in Ys, lam in Ls} := 0.0;
        var ECm1_real {x in Xs, y in Ys, lam in Ls} := 0.0;
        var ECm1_imag {x in Xs, y in Ys, lam in Ls} := 0.0;
        
        subject to st_ECm1_Bmreal_cx {mx in MXs, y in Ys, lam in Ls}: ECm1_Bmreal_cx[mx,y,lam] = 2.*sum {my in MYs: (mx,my) in Mask} EBm_real[mx,my,lam]*cos(2*pi*y*my*(lam0/lam))*dmy;
        subject to st_ECm1_real {(x,y) in Lyot, lam in Ls}: ECm1_real[x,y,lam] = (lam0/lam)*sum {mx in MXs} ECm1_Bmreal_cx[mx,y,lam]*cos(2*pi*x*mx*(lam0/lam))*dmx;
        subject to st_ECm1_imag {(x,y) in Lyot, lam in Ls}: ECm1_imag[x,y,lam] = (lam0/lam)*sum {mx in MXs} ECm1_Bmreal_cx[mx,y,lam]*sin(2*pi*x*mx*(lam0/lam))*dmx;
        
        
        var ECm2_Bmimag_cx {mx in MXs, y in Ys, lam in Ls} := 0.0;
        var ECm2_real {x in Xs, y in Ys, lam in Ls} := 0.0;
        var ECm2_imag {x in Xs, y in Ys, lam in Ls} := 0.0;
        
        subject to st_ECm2_Bmimag_cx {mx in MXs, y in Ys, lam in Ls}: ECm2_Bmimag_cx[mx,y,lam] = 2.*sum {my in MYs: (mx,my) in Mask} EBm_imag[mx,my,lam]*cos(2*pi*y*my*(lam0/lam))*dmy;
        subject to st_ECm2_real {(x,y) in Lyot, lam in Ls}: ECm2_real[x,y,lam] = (lam0/lam)*sum {mx in MYs} ECm2_Bmimag_cx[mx,y,lam]*cos(2*pi*x*mx*(lam0/lam))*dmx;
        subject to st_ECm2_imag {(x,y) in Lyot, lam in Ls}: ECm2_imag[x,y,lam] = (lam0/lam)*sum {mx in MYs} ECm2_Bmimag_cx[mx,y,lam]*sin(2*pi*x*mx*(lam0/lam))*dmx;
        
        
        var ECm_real {x in Xs, y in Ys, lam in Ls} := 0.0;
        var ECm_imag {x in Xs, y in Ys, lam in Ls} := 0.0;
        subject to st_ECm_real {x in Xs, y in Ys, lam in Ls}: ECm_real[x,y,lam] = ECm1_real[x,y,lam]-ECm2_imag[x,y,lam];
        subject to st_ECm_imag {x in Xs, y in Ys, lam in Ls}: ECm_imag[x,y,lam] = ECm1_imag[x,y,lam]+ECm2_real[x,y,lam];
        
        #---------------------
        var ED1_ECmreal_cx {x in Xs, eta in Etas, lam in Ls} := 0.0;
        var ED1_real {xi in Xis, eta in Etas, lam in Ls} := 0.0;
        var ED1_imag {xi in Xis, eta in Etas, lam in Ls} := 0.0;
        
        subject to st_ED1_ECmreal_cx {x in Xs, eta in Etas, lam in Ls}: ED1_ECmreal_cx[x,eta,lam] = 2.*sum {y in Ys: (x,y) in Lyot} (A[x,y]*PupilFile[round(x/dx+0.5+N),round(y/dy+0.5)]-ECm_real[x,y,lam])*cos(2*pi*y*eta*(lam0/lam))*dy;
        subject to st_ED1_real {xi in Xis, eta in Etas, lam in Ls}: ED1_real[xi,eta,lam] = (lam0/lam)*sum {x in Xs} ED1_ECmreal_cx[x,eta,lam]*cos(2*pi*x*xi*(lam0/lam))*dx;
        subject to st_ED1_imag {xi in Xis, eta in Etas, lam in Ls}: ED1_imag[xi,eta,lam] = (lam0/lam)*sum {x in Xs} -ED1_ECmreal_cx[x,eta,lam]*sin(2*pi*x*xi*(lam0/lam))*dx;
        
        
        var ED2_ECmimag_cx {x in Xs, eta in Etas, lam in Ls} := 0.0;
        var ED2_real {xi in Xis, eta in Etas, lam in Ls} := 0.0;
        var ED2_imag {xi in Xis, eta in Etas, lam in Ls} := 0.0;
        
        subject to st_ED2_ECmimag_cx {x in Xs, eta in Etas, lam in Ls}: ED2_ECmimag_cx[x,eta,lam] = 2.*sum {y in Ys: (x,y) in Lyot} (-ECm_imag[x,y,lam])*cos(2*pi*y*eta*(lam0/lam))*dy;
        subject to st_ED2_real {xi in Xis, eta in Etas, lam in Ls}: ED2_real[xi,eta,lam] = (lam0/lam)*sum {x in Xs} ED2_ECmimag_cx[x,eta,lam]*cos(2*pi*x*xi*(lam0/lam))*dx;
        subject to st_ED2_imag {xi in Xis, eta in Etas, lam in Ls}: ED2_imag[xi,eta,lam] = (lam0/lam)*sum {x in Xs} -ED2_ECmimag_cx[x,eta,lam]*sin(2*pi*x*xi*(lam0/lam))*dx;
        
        
        var ED_real {xi in Xis, eta in Etas, lam in Ls} := 0.0;
        var ED_imag {xi in Xis, eta in Etas, lam in Ls} := 0.0;
        subject to st_ED_real {xi in Xis, eta in Etas, lam in Ls}: ED_real[xi,eta,lam] = ED1_real[xi,eta,lam]-ED2_imag[xi,eta,lam];
        subject to st_ED_imag {xi in Xis, eta in Etas, lam in Ls}: ED_imag[xi,eta,lam] = ED1_imag[xi,eta,lam]+ED2_real[xi,eta,lam];
        
        
        #---------------------
        var ED00_real := 0.0;
        subject to st_ED00_real: ED00_real = 2.*sum {x in Xs, y in Ys: (x,y) in Lyot} (A[x,y]*PupilFile[round(x/dx+0.5+N),round(y/dy+0.5)])*dx*dy;
        """

        constraints = """
        #---------------------
        maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;
        
        #subject to sidelobe {(xi,eta) in DarkHole, lam in Ls}: (lam/lam0)^4*(ED_real[xi,eta,lam]^2+ED_imag[xi,eta,lam]^2) <= 10^(-c)*Normterm*I00;
        
        
        subject to sidelobe_zero_real_pos {(xi,eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] <= 10^(-c/2)*ED00_real/sqrt(2.);
        subject to sidelobe_zero_real_neg {(xi,eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] >= -10^(-c/2)*ED00_real/sqrt(2.);
        subject to sidelobe_zero_imag_pos {(xi,eta) in DarkHole, lam in Ls}: ED_imag[xi,eta,lam] <= 10^(-c/2)*ED00_real/sqrt(2.);
        subject to sidelobe_zero_imag_neg {(xi,eta) in DarkHole, lam in Ls}: ED_imag[xi,eta,lam] >= -10^(-c/2)*ED00_real/sqrt(2.);
        """

        misc_options = """
        option times 1;
        option gentimes 1;
        option show_stats 1;
        """
        
        solver = """
        option solver gurobi;
        """

        solver_options = """
        option gurobi_options "outlev=1 lpmethod=2 crossover=0";
        """
 
        execute = """
        solve;

        display solve_result_num, solve_result;
        display ED00_real; 
        """

        store_results = """
        #---------------------

        printf {{x in Xs, y in Ys}}: "%15g %15g %15g \\n", x, y, A[x,y] > "{0:s}"
        """.format(self.fileorg['sol fname'])

        mod_fobj.write( textwrap.dedent(header) )
        mod_fobj.write( textwrap.dedent(params) )
        mod_fobj.write( textwrap.dedent(load_masks) )
        mod_fobj.write( textwrap.dedent(define_coords) )
        mod_fobj.write( textwrap.dedent(sets_and_arrays) )
        mod_fobj.write( textwrap.dedent(field_propagation) )
        mod_fobj.write( textwrap.dedent(constraints) )
        mod_fobj.write( textwrap.dedent(misc_options) )
        mod_fobj.write( textwrap.dedent(solver) )
        mod_fobj.write( textwrap.dedent(solver_options) )
        mod_fobj.write( textwrap.dedent(execute) )
        mod_fobj.write( textwrap.dedent(store_results) )

        mod_fobj.close()
        self.logger.info("Wrote %s"%self.fileorg['ampl src fname'])
