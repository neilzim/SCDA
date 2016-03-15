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
import csv
import numpy as np
import pdb
import getpass
import socket
from collections import defaultdict, OrderedDict
import itertools
import pprint
import pickle

def configure_log(log_fname=None):
#    logger = logging.getLogger("scda.logger")
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
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
    
def load_design_param_survey(pkl_fname):
    fobj = open(pkl_fname, 'rb')
    survey_obj = pickle.load(fobj)
    fobj.close() 
    return survey_obj

class WrappedFixedIndentingLog(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None, style='%', width=70, indent=4):
        super(WrappedFixedIndentingLog, self).__init__(fmt=fmt, datefmt=datefmt)
        self.wrapper = textwrap.TextWrapper(width=width)
        #self.wrapper = textwrap.TextWrapper(width=width, subsequent_indent=' '*indent)
    def format(self, record):
        return self.wrapper.fill(super().format(record))

class DesignParamSurvey(object):
    def __init__(self, coron_class, survey_config, **kwargs):
        #self.logger = logging.getLogger('scda.logger')
        setattr(self, 'coron_class', coron_class)
        self._param_menu = coron_class._design_fields.copy()
        self._file_fields = coron_class._file_fields.copy()
        self._file_fields['fileorg'].append('survey fname')
        setattr(self, 'survey_config', {})
        for keycat, param_dict in survey_config.items():
            self.survey_config[keycat] = {}
            if keycat in self._param_menu:
                for param, values in param_dict.items():
                    if param in self._param_menu[keycat]:
                        if values is not None:
                            if hasattr(values, '__iter__'): #check the type of all items
                                if all(isinstance(value, self._param_menu[keycat][param][0]) for value in values):
                                    self.survey_config[keycat][param] = values
                                    #self.survey_config[keycat][param] = tuple(values)
                                else:
                                    warnstr = ("Warning: Invalid type found in survey set {0} for parameter {1} under category \"{2}\" " + \
                                               "design initialization argument, expecting {3}").format(value, param, keycat, self._param_menu[keycat][param][0]) 
                                    logging.warning(warnstr)
                            else:
                                if isinstance(values, self._param_menu[keycat][param][0]):
                                    self.survey_config[keycat][param] = values
                                else:
                                    warnstr = ("Warning: Invalid {0} for parameter \"{1}\" under category \"{2}\" " + \
                                               "design initialization argument, expecting a {3}").format(type(value), param, keycat, self._param_menu[keycat][param][0]) 
                                    logging.warning(warnstr)
                    else:
                        logging.warning("Warning: Unrecognized parameter \"{0}\" under category \"{1}\" in design initialization argument".format(param, keycat))
            else:
                logging.warning("Warning: Unrecognized key category \"{0}\" in design initialization argument".format(keycat))
                self.survey_config[keycat] = None
        varied_param_flat = []
        varied_param_index = []
        fixed_param_flat = []
        fixed_param_index = []
        for keycat in self._param_menu: # Fill in default values where appropriate
            if keycat not in self.survey_config:
                self.survey_config[keycat] = {}
            for param in self._param_menu[keycat]:
                if param not in self.survey_config[keycat] or (self.survey_config[keycat][param] is None and \
                                                               self._param_menu[keycat][param][1] is not None):
                    self.survey_config[keycat][param] = self._param_menu[keycat][param][1] # default value
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
        self.N_combos = len(varied_param_combos)

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
                                logging.warning("Warning: The specified location of '{0}', \"{1}\" does not exist".format(namekey, self.fileorg[namekey]))
                        else:
                            self.fileorg[namekey] = location
                    else:
                        self.fileorg[namekey] = None
                else:
                    logging.warning("Warning: Unrecognized field {0} in fileorg argument".format(namekey))
        # Handle missing directory values
        if 'work dir' not in self.fileorg or self.fileorg['work dir'] is None:
            self.fileorg['work dir'] = os.getcwd()
        for namekey in self._file_fields['fileorg']: # Set other missing directory locations to 'work dir'
            if namekey.endswith('dir') and ( namekey not in self.fileorg or self.fileorg[namekey] is None ):
                self.fileorg[namekey] = self.fileorg['work dir']
      
        #////////////////////////////////////////////////////////////////////////////////////////////////////
        # In most cases we don't expect to directly specify file names for apertures, FPM, or LS files
        # for the SCDA parameter survey. However, it easy enough to make this option available.
        # If the location of the optimizer input file is not known, 
        # look for it in the directory corresponding to its specific category
        #////////////////////////////////////////////////////////////////////////////////////////////////////
        if 'TelAp fname' in self.fileorg and self.fileorg['TelAp fname'] is not None and \
        not os.path.exists(self.fileorg['TelAp fname']) and os.path.exists(self.fileorg['TelAp dir']) and \
        os.path.dirname(self.fileorg['TelAp fname']) == '':
            try_fname = os.path.join(self.fileorg['TelAp dir'], self.fileorg['TelAp fname']) 
            if os.path.exists(try_fname):
                self.fileorg['TelAp fname'] = try_fname
            else:
                logging.warning("Warning: Could not find the specified telescope aperture file \"{0}\" in {1}".format(self.fileorg['TelAp fname'],
                                    self.fileorg['TelAp dir']))
        if 'FPM fname' in self.fileorg and self.fileorg['FPM fname'] is not None and \
        not os.path.exists(self.fileorg['FPM fname']) and os.path.exists(self.fileorg['FPM dir']) and \
        os.path.dirname(self.fileorg['FPM fname']) == '':
            try_fname = os.path.join(self.fileorg['FPM dir'], self.fileorg['FPM fname']) 
            if os.path.exists(try_fname):
                self.fileorg['FPM fname'] = try_fname
            else:
                logging.warning("Warning: Could not find the specified FPM file \"{0}\" in {1}".format(self.fileorg['FPM fname'],
                                    self.fileorg['FPM dir']))
        if 'LS fname' in self.fileorg and self.fileorg['LS fname'] is not None and \
        not os.path.exists(self.fileorg['LS fname']) and os.path.exists(self.fileorg['LS dir']) and \
        os.path.dirname(self.fileorg['LS fname']) == '':
            try_fname = os.path.join(self.fileorg['LS dir'], self.fileorg['LS fname']) 
            if os.path.exists(try_fname):
                self.fileorg['LS fname'] = try_fname
            else:
                logging.warning("Warning: Could not find the specified LS file \"{0}\" in {1}".format(self.fileorg['LS fname'],
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
                        logging.warning("Warning: Unrecognized solver option \"{0}\" in field \"{1}\", reverting to default".format(value, field))
                else:
                    logging.warning("Warning: Unrecognized field {0} in solver argument".format(field))
        # Handle missing values
        if 'constr' not in self.solver or self.solver['constr'] is None: self.solver['constr'] = 'lin'
        if 'method' not in self.solver or self.solver['method'] is None: self.solver['method'] = 'bar'
        if 'presolve' not in self.solver or self.solver['presolve'] is None: self.solver['presolve'] = True
        if 'threads' not in self.solver or self.solver['threads'] is None: self.solver['threads'] = None
         
        setattr(self, 'coron_list', [])
        design = {}
        for keycat in self._param_menu:
            design[keycat] = {}
        for (fixed_keycat, fixed_parname), fixed_val in zip(self.fixed_param_index, self.fixed_param_vals):
            design[fixed_keycat][fixed_parname] = fixed_val
        self.coron_list = []
        for param_combo in self.varied_param_combos: # TODO: Switch the coronagraph type depending on the symmetry of the telescope aperture and support struts 
            for (varied_keycat, varied_parname), current_val in zip(self.varied_param_index, param_combo):
                design[varied_keycat][varied_parname] = current_val
            coron_fileorg = self.fileorg.copy()
            if 'survey fname' in coron_fileorg:
                coron_fileorg.pop('survey fname')
            self.coron_list.append( coron_class(design=design, fileorg=self.fileorg, solver=self.solver) )
 
        setattr(self, 'ampl_infile_status', None)
        self.check_ampl_input_files()
        setattr(self, 'ampl_src_status', None)
        setattr(self, 'solution_status', None)
        setattr(self, 'evaluation_status', None)

    def write_ampl_batch(self, overwrite=False, override_infile_status=False):
        for coron in self.coron_list:
            coron.write_ampl(overwrite, override_infile_status, verbose=False)
        logging.info("Wrote the batch of design survey AMPL programs into {:s}".format(self.fileorg['ampl src dir']))
 
    def check_ampl_input_files(self):
        survey_status = True
        for coron in self.coron_list: # Update all individual statuses
            coron_status = coron.check_ampl_input_files()
            if coron_status is False: # If one is missing input files, set the survey-wide input file status to False
                survey_status = False
        self.ampl_infile_status = survey_status
        return survey_status

    def check_ampl_src_files(self):
        status = True
        for coron in self.coron_list:
            if not os.path.exists(coron.fileorg['ampl src fname']):
                status = False
                break
        self.ampl_src_status = status
        return status

    def check_solution_files(self):
        status = True
        for coron in self.coron_list:
            if not os.path.exists(coron.fileorg['sol fname']):
                status = False
                break
        self.solution_status = status
        return status

    def check_eval_status(self):
        status = True
        for coron in self.coron_list: # Update all individual statuses
            if coron.eval_metrics['thrupt'] is None or coron.eval_metrics['psf area'] is None \
               or coron.eval_metrics['psf width'] is None:
                status = False
                break
        self.eval_status = status
        return status

    def write(self, fname=None):
        if fname is not None:
            if os.path.dirname(fname) is '': # if no path specified, assume work dir
                self.fileorg['survey fname'] = os.path.join(self.fileorg['work dir'], fname)
            else:
                self.fileorg['survey fname'] = fname
        else:
            if 'survey fname' not in self.fileorg or \
               ('survey fname' in self.fileorg and self.fileorg['survey fname'] is None): # set the filename based on the coronagraph type, user, and date
                fname_tail = "scda_{:s}_survey_{:s}_{:s}.pkl".format(self.coron_class.__name__, getpass.getuser(), datetime.datetime.now().strftime("%Y-%m-%d"))
                self.fileorg['survey fname'] = os.path.join(self.fileorg['work dir'], fname_tail)
        fobj = open(self.fileorg['survey fname'], 'wb')
        pickle.dump(self, fobj)
        fobj.close()
        logging.info("Wrote the design parameter survey object to {:s}".format(self.fileorg['survey fname']))
 
    def write_spreadsheet(self, overwrite=False, csv_fname=None):
        if csv_fname is not None:
            if os.path.dirname(csv_fname) is '': # if no path specified, assume work dir
                csv_fname = os.path.join(self.fileorg['work dir'], fname)
            else:
                csv_fname = fname
        else:
            if 'survey fname' not in self.fileorg or ('survey fname' in self.fileorg and self.fileorg['survey fname'] is None):
                csv_fname_tail = "scda_{:s}_survey_{:s}_{:s}.csv".format(self.coron_class.__name__, getpass.getuser(), datetime.datetime.now().strftime("%Y-%m-%d"))
                csv_fname = os.path.join(self.fileorg['work dir'], csv_fname_tail)
            else:
                csv_fname = self.fileorg['survey fname'][-4:] + ".csv"
        with open(csv_fname, 'wb') as survey_spreadsheet:
            self.check_ampl_src_files()
            self.check_ampl_input_files()
            self.check_solution_files()
            self.check_eval_status()
            surveywriter = csv.writer(survey_spreadsheet)
            #/////////////////////////////////////////////////
            #    Write a header for the spreadsheet
            #/////////////////////////////////////////////////
            surveywriter.writerow(["Created by {:s} on {:s} at {:s}".format(getpass.getuser(), socket.gethostname(), datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))])
            surveywriter.writerow(["FILE ORGANIZATION AND STATUS"])
            surveywriter.writerow(["Work dir", self.fileorg['work dir']])
            surveywriter.writerow(["AMPL source location", self.fileorg['ampl src dir']])
            surveywriter.writerow(["Solution location", self.fileorg['sol dir']])
            surveywriter.writerow(["Telescope aperture location", self.fileorg['TelAp dir']])
            surveywriter.writerow(["Focal plane mask location", self.fileorg['FPM dir']])
            surveywriter.writerow(["Lyot stop location", self.fileorg['LS dir']])
            if self.ampl_src_status is True:
                surveywriter.writerow(["All AMPL source files exist?", 'Y'])
            else:
                surveywriter.writerow(["All AMPL source files exist?", 'N'])
            if self.ampl_infile_status is True:
                surveywriter.writerow(["All input files exist?", 'Y'])
            else:
                surveywriter.writerow(["All input files exist?", 'N'])
            if self.solution_status is True:
                surveywriter.writerow(["All solution files exist?", 'Y'])
            else:
                surveywriter.writerow(["All solution files exist?", 'N'])
            if self.eval_status is True:
                surveywriter.writerow(["All evaluation metrics extracted?", 'Y'])
            else:
                surveywriter.writerow(["All evaluation metrics extracted?", 'N'])
            #/////////////////////////////////////////////////
            #    Write out the fixed design parameters
            #/////////////////////////////////////////////////
            surveywriter.writerow([""])
            surveywriter.writerow(["FIXED design parameters"])
            fixed_param_category_row = []
            fixed_param_subheading_row = []
            for cat in self._param_menu:
                fixed_param_category_row.extend([cat,''])
                fixed_param_subheading_row.extend(['param name', 'value'])
            surveywriter.writerow(fixed_param_category_row)
            surveywriter.writerow(fixed_param_subheading_row)
            self.fixed_param_table = []
            max_N_params = 0
            for cat in self._param_menu:
                param_col = []
                val_col = []
                for (param_cat, param) in self.fixed_param_index:
                    if param_cat is cat:
                        param_col.append(param)
                        val_col.append(self.survey_config[cat][param])
                N_params = len(param_col)
                if N_params > max_N_params:
                    max_N_params = N_params
                self.fixed_param_table.append(param_col)
                self.fixed_param_table.append(val_col)
            N_cols = len(self.fixed_param_table)
            for ci in range(N_cols):
                N_rows = len(self.fixed_param_table[ci])
                if N_rows < max_N_params:
                    self.fixed_param_table[ci].extend(['']*(max_N_params - N_rows))
            for ri in range(max_N_params):
                fixed_table_row = []
                for ci in range(N_cols):
                    fixed_table_row.append(self.fixed_param_table[ci][ri])
                surveywriter.writerow(fixed_table_row)
            #/////////////////////////////////////////////////
            #    Write out the varied design parameters
            #/////////////////////////////////////////////////
            surveywriter.writerow([""])
            surveywriter.writerow(["VARIED design parameters ({:d} total combinations)".format(self.N_combos)])
            varied_param_category_row = []
            varied_param_subheading_row = []
            for cat in self._param_menu:
                varied_param_category_row.extend([cat,'',''])
                varied_param_subheading_row.extend(['param name', 'value list', 'num'])
            surveywriter.writerow(varied_param_category_row)
            surveywriter.writerow(varied_param_subheading_row)
            self.varied_param_table = []
            max_N_params = 0
            for cat in self._param_menu:
                param_col = []
                vals_col = []
                num_col = []
                for (param_cat, param) in self.varied_param_index:
                    if param_cat is cat:
                        param_col.append(param)
                        vals_col.append(self.survey_config[cat][param])
                        num_col.append(len(self.survey_config[cat][param]))
                N_params = len(param_col)
                if N_params > max_N_params:
                    max_N_params = N_params
                self.varied_param_table.append(param_col)
                self.varied_param_table.append(vals_col)
                self.varied_param_table.append(num_col)
            N_cols = len(self.varied_param_table)
            for ci in range(N_cols):
                N_rows = len(self.varied_param_table[ci])
                if N_rows < max_N_params:
                    self.varied_param_table[ci].extend(['']*(max_N_params - N_rows))
            for ri in range(max_N_params):
                varied_table_row = []
                for ci in range(N_cols):
                    varied_table_row.append(self.varied_param_table[ci][ri])
                surveywriter.writerow(varied_table_row)
            #/////////////////////////////////////////////////////////
            #    Write out the survey design parameter combinations
            #/////////////////////////////////////////////////////////
            surveywriter.writerow([""])
            surveywriter.writerow(["SURVEY TABLE"])
            catrow = []
            paramrow = []
            if len(self.varied_param_index) > 0:
                for (cat, name) in self.varied_param_index:
                    catrow.append(cat)
                    paramrow.append(name)
                catrow.extend(['', 'AMPL source', '', '', 'Solution', '', 'Evaluation metrics', '', ''])
                paramrow.extend(['', 'filename', 'exists?', 'input files?', 'filename', 'exists?', 'Thrupt', 'PSF area', 'PSF width'])
                surveywriter.writerow(catrow)
                surveywriter.writerow(paramrow)
                for ii, param_combo in enumerate(self.varied_param_combos):
                    param_combo_row = list(param_combo)
                    param_combo_row.append('')
                    param_combo_row.append(os.path.basename(self.coron_list[ii].fileorg['ampl src fname'])) 
                    if os.path.exists(self.coron_list[ii].fileorg['ampl src fname']):
                        param_combo_row.append('Y')
                    else:
                        param_combo_row.append('N')
                    if self.coron_list[ii].ampl_infile_status is True:
                        param_combo_row.append('Y')
                    else:
                        param_combo_row.append('N')
                    param_combo_row.append(os.path.basename(self.coron_list[ii].fileorg['sol fname'])) 
                    if os.path.exists(self.coron_list[ii].fileorg['sol fname']):
                        param_combo_row.append('Y')
                    else:
                        param_combo_row.append('N')
                    surveywriter.writerow(param_combo_row)
                    
        survey_spreadsheet.close()
        logging.info("Wrote design survey spreadsheet to {:s}".format(csv_fname))

class LyotCoronagraph(object): # Lyot coronagraph base class
    _file_fields = { 'fileorg': ['work dir', 'ampl src dir', 'TelAp dir', 'FPM dir', 'LS dir', 'sol dir', 'eval dir',
                                 'ampl src fname', 'TelAp fname', 'FPM fname', 'LS fname', 'sol fname'],
                     'solver': ['constr', 'method', 'presolve', 'threads'] }

    _solver_menu = { 'constr': ['lin', 'quad'], 'solver': ['LOQO', 'gurobi', 'gurobix'], 
                     'method': ['bar', 'barhom', 'dualsimp'],
                     'presolve': [True, False], 'threads': [None]+range(1,33) }

    _aperture_menu = { 'pm': ['hex1', 'hex2', 'hex3', 'key24', 'pie12', 'pie08', 'irisao', 'atlast'],
                       'ss': ['Y60d','Yoff60d','X','Cross','T','Y90d'],
                       'sst': ['025','100'],
                       'sm': [True, False] }

    def __init__(self, verbose=False, **kwargs):
        # Only set fileorg and solver attributes in this constructor,
        # since design and eval parameter checking is design-specific. 

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
                                logging.warning("Warning: The specified location of '{0}', \"{1}\" does not exist".format(namekey, self.fileorg[namekey]))
                        else:
                            self.fileorg[namekey] = location
                    else:
                        self.fileorg[namekey] = None
                else:
                    logging.warning("Warning: Unrecognized field {0} in fileorg argument".format(namekey))
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
                logging.warning("Warning: Could not find the specified telescope aperture file \"{0}\" in {1}".format(self.fileorg['TelAp fname'], \
                                    self.fileorg['TelAp dir']))
        if 'FPM fname' in self.fileorg and self.fileorg['FPM fname'] is not None and \
        not os.path.exists(self.fileorg['FPM fname']) and os.path.exists(self.fileorg['FPM dir']) and \
        os.path.dirname(self.fileorg['FPM fname']) == '':
            try_fname = os.path.join(self.fileorg['FPM dir'], self.fileorg['FPM fname']) 
            if os.path.exists(try_fname):
                self.fileorg['FPM fname'] = try_fname
            else:
                logging.warning("Warning: Could not find the specified FPM file \"{0}\" in {1}".format(self.fileorg['FPM fname'], \
                                    self.fileorg['FPM dir']))
        if 'LS fname' in self.fileorg and self.fileorg['LS fname'] is not None and \
        not os.path.exists(self.fileorg['LS fname']) and os.path.exists(self.fileorg['LS dir']) and \
        os.path.dirname(self.fileorg['LS fname']) == '':
            try_fname = os.path.join(self.fileorg['LS dir'], self.fileorg['LS fname']) 
            if os.path.exists(try_fname):
                self.fileorg['LS fname'] = try_fname
            else:
                logging.warning("Warning: Could not find the specified LS file \"{0}\" in {1}".format(self.fileorg['LS fname'], \
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
                        logging.warning("Warning: Unrecognized solver option \"{0}\" in field \"{1}\", reverting to default".format(value, field))
                else:
                    logging.warning("Warning: Unrecognized field {0} in solver argument".format(field))
        # Handle missing values
        if 'constr' not in self.solver or self.solver['constr'] is None: self.solver['constr'] = 'lin'
        if 'solver' not in self.solver or self.solver['solver'] is None: self.solver['solver'] = 'gurobi'
        if 'method' not in self.solver or self.solver['method'] is None: self.solver['method'] = 'bar'
        if 'presolve' not in self.solver or self.solver['presolve'] is None: self.solver['presolve'] = True
        if 'threads' not in self.solver or self.solver['threads'] is None: self.solver['threads'] = None

        setattr(self, 'ampl_infile_status', None)
        if not issubclass(self.__class__, LyotCoronagraph):
            self.check_ampl_input_files()

        setattr(self, 'eval_metrics', {})
        self.eval_metrics['thrupt'] = None
        self.eval_metrics['psf area'] = None
        self.eval_metrics['psf width'] = None

    def check_ampl_input_files(self):
        status = True
        checklist = ['TelAp fname', 'FPM fname', 'LS fname']
        for fname in checklist:
            if not os.path.exists(self.fileorg[fname]):
                status = False
                break
        self.ampl_infile_status = status
        return status

class NdiayeAPLC(LyotCoronagraph): # Image-constrained APLC following N'Diaye et al. (2015, 2016)
    _design_fields = OrderedDict([ ( 'Pupil', OrderedDict([('N',(int, 1000)), ('pm',(str, 'hex1')), ('ss',(str, 'x')), 
                                                          ('sst',(str, '100')), ('sm',(bool, True))]) ),
                                   ( 'FPM', OrderedDict([('rad',(float, 4.)), ('M',(int, 50))]) ),
                                   ( 'LS', OrderedDict([('shape',(str, 'ann')), ('id',(int, 20)), ('od',(int, 90)), ('obscure',(int, 0)),
                                                        ('spad',(int, 0)), ('ppad',(int, 0)), ('altol',(int, None))]) ),
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
                                    logging.warning(warnstr)
                        else:
                            logging.warning("Warning: Unrecognized parameter \"{0}\" under category \"{1}\" in design initialization argument".format(param, keycat))
                else:
                    logging.warning("Warning: Unrecognized key category \"{0}\" in design initialization argument".format(keycat))
                    self.design[keycat] = None
        for keycat in self._design_fields: # Fill in default values where appropriate
            if keycat not in self.design:
                self.design[keycat] = {}
            for param in self._design_fields[keycat]:
                if param not in self.design[keycat] or (self.design[keycat][param] is None and \
                                                        self._design_fields[keycat][param][1] is not None):
                    self.design[keycat][param] = self._design_fields[keycat][param][1]
        if self.design['Image']['co'] is None: # An OCA different from OWA is meaningless without a corresponding outer contrast constraint
            self.design['Image']['oca'] = self.design['Image']['owa']
        # Finally, set a private attribute for the number of image plane samples between the center and the OCA
        self.design['Image']['_Nimg'] = int( np.ceil( self.design['Image']['fpres']*self.design['Image']['oca']/(1. - self.design['Image']['bw']/2) ) )
        if verbose: # Print summary of the set parameters
            logging.info("Design parameters: {}".format(self.design))
            logging.info("Optimization and solver parameters: {}".format(self.solver))
            logging.info("File organization parameters: {}".format(self.fileorg))
     
        self.amplname_coron = "APLC_full"
        self.amplname_pupil = "{0:s}{1:s}{2:s}sm{3:d}_N{4:04d}".format(self.design['Pupil']['pm'], self.design['Pupil']['ss'], self.design['Pupil']['sst'], \
                                                                       int(self.design['Pupil']['sm']), self.design['Pupil']['N'])

        self.amplname_fpm = "FPM{:02}M{:03}".format(int(round(100*self.design['FPM']['rad'])), self.design['FPM']['M'])
        if self.design['LS']['obscure'] == 2: # LS includes primary and secondary aperture features
            self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}{3:s}Pad{4:02d}{5:s}{6:s}sm{7:d}Pad{8:02d}".format(self.design['LS']['shape'], self.design['LS']['id'], \
                               self.design['LS']['od'], self.design['Pupil']['pm'], self.design['LS']['ppad'], self.design['Pupil']['ss'], self.design['Pupil']['sst'], \
                               int(self.design['Pupil']['sm']), self.design['LS']['spad'])
        elif self.design['LS']['obscure'] == 1: # LS includes secondary aperture features
            self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}{3:s}{4:s}sm{5:d}Pad{6:02d}".format(self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'], \
                               self.design['Pupil']['ss'], self.design['Pupil']['sst'], int(self.design['Pupil']['sm']), self.design['LS']['spad'])
        else: # LS aperture is unobscured
            self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}clear".format(self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'])
        if self.design['LS']['altol'] is not None:
            self.amplname_ls += "tol{0:02d}".format(self.design['LS']['altol'])

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
        if self.solver['threads'] is not None:
            self.amplname_solver += "thr{:02d}".format(self.solver['threads'])

        if not issubclass(self.__class__, NdiayeAPLC): # Only set these file names if this is a full-plane APLC.
            if 'ampl src fname' not in self.fileorg or self.fileorg['ampl src fname'] is None:
                ampl_src_fname_tail = self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                                      self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".mod"
                self.fileorg['ampl src fname'] = os.path.join(self.fileorg['ampl src dir'], ampl_src_fname_tail)

            if 'sol fname' not in self.fileorg or self.fileorg['sol fname'] is None:
                self.fileorg['sol fname'] = self.fileorg['ampl src fname'][:-4] + "_ApodSol.dat"

            if 'TelAp fname' not in self.fileorg or self.fileorg['TelAp fname'] is None:
                self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_full_" + self.amplname_pupil + ".dat") )

            if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
                self.fileorg['FPM fname'] = os.path.join( self.fileorg['FPM dir'], "FPM_full_occspot_M{:03}.dat".format(self.design['FPM']['M']) )

            if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_full_" + self.amplname_pupil + \
                                                         "_{0:02d}D{1:02d}ovsz{2:02d}.dat".format(self.design['LS']['id'], \
                                                         self.design['LS']['od'], self.design['LS']['ovsz'])) )
            self.check_ampl_input_files()
                                                              
    def write_ampl(self, overwrite=False):
        logging.info("Writing the AMPL program")
    def read_solution(self):
        logging.info("Reading in the apodizer solution and parse the optimizer log")
    def create_eval_model(self):
        logging.info("Defining a python/Poppy model to evaluate the solution")
    def eval_solution(self):
        logging.info("Evaluating the design throughput, PSF FWHM, etc., and writing summary")

class HalfplaneAPLC(NdiayeAPLC): # N'Diaye APLC subclass for the half-plane symmetry case
    def __init__(self, **kwargs):
        super(HalfplaneAPLC, self).__init__(**kwargs)
        self.amplname_coron = "APLC_half"
        if 'ampl src fname' not in self.fileorg or self.fileorg['ampl src fname'] is None:
            ampl_src_fname_tail = self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                                  self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".mod"
            self.fileorg['ampl src fname'] = os.path.join(self.fileorg['ampl src dir'], ampl_src_fname_tail)

        if 'sol fname' not in self.fileorg or self.fileorg['sol fname'] is None:
            self.fileorg['sol fname'] = self.fileorg['ampl src fname'][:-4] + "_ApodSol.dat"

        if 'TelAp fname' not in self.fileorg or self.fileorg['TelAp fname'] is None:
            self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_half_" + self.amplname_pupil + ".dat") )

        if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
            self.fileorg['FPM fname'] = os.path.join( self.fileorg['FPM dir'], "FPM_half_occspot_M{:03}.dat".format(self.design['FPM']['M']) )

        if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}Pad{4:02d}{5:s}{6:s}sm{7:d}Pad{8:02d}_N{9:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['pm'], self.design['LS']['ppad'], self.design['Pupil']['ss'],
                                                         self.design['Pupil']['sst'], int(self.design['Pupil']['sm']), self.design['LS']['spad'],
                                                         self.design['Pupil']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}sm{5:d}Pad{6:02d}_N{7:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['ss'], self.design['Pupil']['sst'], int(self.design['Pupil']['sm']),
                                                         self.design['LS']['spad'], self.design['Pupil']['N'])) )
            else:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_clear_N{3:04d}.dat".format(self.design['LS']['shape'],
                                                         self.design['LS']['id'], self.design['LS']['od'], self.design['Pupil']['N'])) )
        self.check_ampl_input_files()
    def write_ampl(self, overwrite=False, override_infile_status=False, ampl_src_fname=None, verbose=True):
        if self.ampl_infile_status is False and not override_infile_status:
            logging.error("Error: the most recent input file check for this design configuration failed.")
            logging.error("The override_infile_status switch is off, so write_ampl() will now abort.")
            logging.error("See previous warnings in the log to see what file was missing during the initialization")
            return
        if ampl_src_fname is not None:
            if os.path.dirname(ampl_src_fname) == '' and self.fileorg['ampl src dir'] is not None:
                self.fileorg['ampl src fname'] = os.path.join(self.fileorg['ampl src dir'], ampl_src_fname)
            else:
                self.fileorg['ampl src fname'] = os.path.abspath(ampl_src_fname)    
                self.fileorg['ampl src dir'] = os.path.dirname(self.fileorg['ampl src fname'])
        if os.path.exists(self.fileorg['ampl src fname']):
            if overwrite == True:
                if verbose:
                    logging.warning("Warning: Overwriting the existing copy of {0}".format(self.fileorg['ampl src fname']))
            else:
                logging.warning("Error: {0} already exists and overwrite switch is off, so write_ampl() will now abort".format(self.fileorg['ampl src fname']))
                return
        elif not os.path.exists(self.fileorg['ampl src dir']):
            os.mkdir(self.fileorg['ampl src dir'])
            if verbose:
                logging.info("Created new AMPL source code directory, {0:s}".format(self.fileorg['ampl src dir']))
        mod_fobj = open(self.fileorg['ampl src fname'], "w")

        header = """\
        # AMPL program to optimize a half-plane symmetric APLC
        # Created by {0:s} with {1:s} on {2:s} at {3:s}
        load amplgsl.dll;


        #---------------------

        param pi:= 4*atan(1);

        #---------------------
        param c := {4:.2f};

        #---------------------
        param Rmask := {5:0.3f};
        param rho0 := {6:0.2f};
        param rho1 := {7:0.2f};
        
        #---------------------
        param N := {8:d};				# discretization parameter (pupil)
        param M := {9:d};				# discretization parameter (mask)
        
        param Nimg := {10:d};           # discretization parameter (image)
        param rho2 := {11:0.2f};        # NTZ: We paramaterize our image plane resolution by fpres = sampling rate at 
                                    #      the shortest wavelength. Then Nimg is an integer function of fpres, oca,
        #---------------------      #      and bw. This integer is not specified directly by the user, but computed "privately"
        param bw := {12:0.2f};           #      by the APLC class constructor.
        param lam0 := 1.;
        param dl := bw*lam0;
        param Nlam := {13:d};
        
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
        """.format(getpass.getuser(), os.path.basename(__file__), socket.gethostname(), datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), \
                   self.design['Image']['c'], self.design['FPM']['rad'], self.design['Image']['iwa'], self.design['Image']['owa'], \
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
        
        param dxi := (rho2/Nimg)*(1/CoeffOverSizePup);
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
        set PSFCore := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= 0 && sqrt(xi^2+eta^2) < rho0} (xi,eta); # Only for 360deg masks.
        set InsideArea := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= 0 && sqrt(xi^2+eta^2) <= rho1} (xi,eta); # Only for 360deg masks.
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
        
        #subject to st_ED1_real {xi in Xis, eta in Etas, lam in Ls}: ED1_real[xi,eta,lam] = (lam0/lam)*sum {x in Xs} ED1_ECmreal_cx[x,eta,lam]*cos(2*pi*x*xi*(lam0/lam))*dx;
        #subject to st_ED1_imag {xi in Xis, eta in Etas, lam in Ls}: ED1_imag[xi,eta,lam] = (lam0/lam)*sum {x in Xs} -ED1_ECmreal_cx[x,eta,lam]*sin(2*pi*x*xi*(lam0/lam))*dx;
        
        subject to st_ED1_real {(xi, eta) in InsideArea, lam in Ls}: ED1_real[xi,eta,lam] = (lam0/lam)*sum {x in Xs} ED1_ECmreal_cx[x,eta,lam]*cos(2*pi*x*xi*(lam0/lam))*dx;
        subject to st_ED1_imag {(xi, eta) in InsideArea, lam in Ls}: ED1_imag[xi,eta,lam] = (lam0/lam)*sum {x in Xs} -ED1_ECmreal_cx[x,eta,lam]*sin(2*pi*x*xi*(lam0/lam))*dx;
        
        
        var ED2_ECmimag_cx {x in Xs, eta in Etas, lam in Ls} := 0.0;
        var ED2_real {xi in Xis, eta in Etas, lam in Ls} := 0.0;
        var ED2_imag {xi in Xis, eta in Etas, lam in Ls} := 0.0;
        
        subject to st_ED2_ECmimag_cx {x in Xs, eta in Etas, lam in Ls}: ED2_ECmimag_cx[x,eta,lam] = 2.*sum {y in Ys: (x,y) in Lyot} (-ECm_imag[x,y,lam])*cos(2*pi*y*eta*(lam0/lam))*dy;
        
        #subject to st_ED2_real {xi in Xis, eta in Etas, lam in Ls}: ED2_real[xi,eta,lam] = (lam0/lam)*sum {x in Xs} ED2_ECmimag_cx[x,eta,lam]*cos(2*pi*x*xi*(lam0/lam))*dx;
        #subject to st_ED2_imag {xi in Xis, eta in Etas, lam in Ls}: ED2_imag[xi,eta,lam] = (lam0/lam)*sum {x in Xs} -ED2_ECmimag_cx[x,eta,lam]*sin(2*pi*x*xi*(lam0/lam))*dx;
        
        subject to st_ED2_real {(xi, eta) in InsideArea, lam in Ls}: ED2_real[xi,eta,lam] = (lam0/lam)*sum {x in Xs} ED2_ECmimag_cx[x,eta,lam]*cos(2*pi*x*xi*(lam0/lam))*dx;
        subject to st_ED2_imag {(xi, eta) in InsideArea, lam in Ls}: ED2_imag[xi,eta,lam] = (lam0/lam)*sum {x in Xs} -ED2_ECmimag_cx[x,eta,lam]*sin(2*pi*x*xi*(lam0/lam))*dx;
        
        
        var ED_real {xi in Xis, eta in Etas, lam in Ls} := 0.0;
        var ED_imag {xi in Xis, eta in Etas, lam in Ls} := 0.0;
        #subject to st_ED_real {xi in Xis, eta in Etas, lam in Ls}: ED_real[xi,eta,lam] = ED1_real[xi,eta,lam]-ED2_imag[xi,eta,lam];
        #subject to st_ED_imag {xi in Xis, eta in Etas, lam in Ls}: ED_imag[xi,eta,lam] = ED1_imag[xi,eta,lam]+ED2_real[xi,eta,lam];
        
        subject to st_ED_real {(xi, eta) in InsideArea, lam in Ls}: ED_real[xi,eta,lam] = ED1_real[xi,eta,lam]-ED2_imag[xi,eta,lam];
        subject to st_ED_imag {(xi, eta) in InsideArea, lam in Ls}: ED_imag[xi,eta,lam] = ED1_imag[xi,eta,lam]+ED2_real[xi,eta,lam];
        
        
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
        
        subject to sidelobe_zero_real_pos_core {(xi,eta) in PSFCore, lam in Ls}: ED_real[xi,eta,lam] <= 10^((-c+2)/2)*ED00_real/sqrt(2.);
        subject to sidelobe_zero_real_neg_core {(xi,eta) in PSFCore, lam in Ls}: ED_real[xi,eta,lam] >= -10^((-c+2)/2)*ED00_real/sqrt(2.);
        subject to sidelobe_zero_imag_pos_core {(xi,eta) in PSFCore, lam in Ls}: ED_imag[xi,eta,lam] <= 10^((-c+2)/2)*ED00_real/sqrt(2.);
        subject to sidelobe_zero_imag_neg_core {(xi,eta) in PSFCore, lam in Ls}: ED_imag[xi,eta,lam] >= -10^((-c+2)/2)*ED00_real/sqrt(2.);
        """
        
        solver = """
        option solver gurobi;
        """

        solver_options = """
        option gurobi_options "outlev=1 lpmethod=2 crossover=0";
        """

        misc_options = """
        option times 1;
        option gentimes 1;
        option show_stats 1;
        """
 
        execute = """
        solve;

        display solve_result_num, solve_result;
        display ED00_real; 
        """

        store_results = """
        #---------------------

        printf {{x in Xs, y in Ys}}: "%15g %15g %15g \\n", x, y, A[x,y] > "{0:s};"
        """.format(self.fileorg['sol fname'])

        mod_fobj.write( textwrap.dedent(header) )
        mod_fobj.write( textwrap.dedent(load_masks) )
        mod_fobj.write( textwrap.dedent(define_coords) )
        mod_fobj.write( textwrap.dedent(sets_and_arrays) )
        mod_fobj.write( textwrap.dedent(field_propagation) )
        mod_fobj.write( textwrap.dedent(constraints) )
        mod_fobj.write( textwrap.dedent(solver) )
        mod_fobj.write( textwrap.dedent(solver_options) )
        mod_fobj.write( textwrap.dedent(misc_options) )
        mod_fobj.write( textwrap.dedent(execute) )
        mod_fobj.write( textwrap.dedent(store_results) )

        mod_fobj.close()
        if verbose:
            logging.info("Wrote %s"%self.fileorg['ampl src fname'])

class QuarterplaneAPLC(NdiayeAPLC): # N'Diaye APLC subclass for the quarter-plane symmetry case
     def __init__(self, **kwargs):
         super(QuarterplaneAPLC, self).__init__(**kwargs)
         self.amplname_coron = "APLC_quart"
         if 'ampl src fname' not in self.fileorg or self.fileorg['ampl src fname'] is None:
             ampl_src_fname_tail = self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                                   self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".mod"
             self.fileorg['ampl src fname'] = os.path.join(self.fileorg['ampl src dir'], ampl_src_fname_tail)
 
         if 'sol fname' not in self.fileorg or self.fileorg['sol fname'] is None:
             sol_fname_tail = "ApodSol_" + self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                              self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".dat"
             self.fileorg['sol fname'] = os.path.join(self.fileorg['ampl src dir'], sol_fname_tail)
 
         if 'TelAp fname' not in self.fileorg or self.fileorg['TelAp fname'] is None:
             self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_quart_" + self.amplname_pupil + ".dat") )
 
         if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
             self.fileorg['FPM fname'] = os.path.join( self.fileorg['FPM dir'], "FPM_quart_occspot_M{:03}.dat".format(self.design['FPM']['M']) )

         if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
             if self.design['LS']['obscure'] == 2:
                 self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_quart_" + \
                                                          "{0:s}{1:02d}D{2:02d}_{3:s}Pad{4:02d}{5:s}{6:s}sm{7:d}Pad{8:02d}_N{9:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['Pupil']['pm'], self.design['LS']['ppad'], self.design['Pupil']['ss'],
                                                          self.design['Pupil']['sst'], int(self.design['Pupil']['sm']), self.design['LS']['spad'],
                                                          self.design['Pupil']['N'])) )
             elif self.design['LS']['obscure'] == 1:
                 self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_quart_" + \
                                                          "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}sm{5:d}Pad{6:02d}_N{7:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['Pupil']['ss'], self.design['Pupil']['sst'], int(self.design['Pupil']['sm']),
                                                          self.design['LS']['spad'], self.design['Pupil']['N'])) )
             else:
                 self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_quart_" + \
                                                          "{0:s}{1:02d}D{2:02d}_clear_N{3:04d}.dat".format(self.design['LS']['shape'],
                                                          self.design['LS']['id'], self.design['LS']['od'], self.design['Pupil']['N'])) )
         self.check_ampl_input_files()
     def write_ampl(self, overwrite=False, override_infile_status=False, ampl_src_fname=None, verbose=True):
         if self.ampl_infile_status is False and not override_infile_status:
             logging.error("Error: the most recent input file check for this design configuration failed.")
             logging.error("The override_infile_status switch is off, so write_ampl() will now abort.")
             logging.error("See previous warnings in the log to see what file was missing during the initialization")
             return
         if ampl_src_fname is not None:
             if os.path.dirname(ampl_src_fname) == '' and self.fileorg['ampl src dir'] is not None:
                 self.fileorg['ampl src fname'] = os.path.join(self.fileorg['ampl src dir'], ampl_src_fname)
             else:
                 self.fileorg['ampl src fname'] = os.path.abspath(ampl_src_fname)    
                 self.fileorg['ampl src dir'] = os.path.dirname(self.fileorg['ampl src fname'])
         if os.path.exists(self.fileorg['ampl src fname']):
             if overwrite == True:
                 if verbose:
                     logging.warning("Warning: Overwriting the existing copy of {0}".format(self.fileorg['ampl src fname']))
             else:
                 logging.warning("Error: {0} already exists and overwrite switch is off, so write_ampl() will now abort".format(self.fileorg['ampl src fname']))
                 return
         elif not os.path.exists(self.fileorg['ampl src dir']):
             os.mkdir(self.fileorg['ampl src dir'])
             if verbose:
                 logging.info("Created new AMPL source code directory, {0:s}".format(self.fileorg['ampl src dir']))
         mod_fobj = open(self.fileorg['ampl src fname'], "w")
 
         header = """\
         # AMPL program to optimize a quarter-plane symmetric APLC
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
         param rho2 := {7:0.2f};        # NTZ: We paramaterize our image plane resolution by fpres = sampling rate at 
                                     #      the shortest wavelength. Then Nimg is an integer function of fpres, oca,
         #---------------------      #      and bw. This integer is not specified directly by the user, but computed "privately"
         param bw := {8:0.2f};           #      by the APLC class constructor.
         param lam0 := 1.;
         param dl := bw*lam0;
         param Nlam := {9:d};
         
         #---------------------
         """.format(self.design['Image']['c'], self.design['FPM']['rad'], self.design['Image']['iwa'], self.design['Image']['owa'], \
                    self.design['Pupil']['N'], self.design['FPM']['M'], self.design['Image']['_Nimg'], \
                    self.design['Image']['oca'], self.design['Image']['bw'], self.design['Image']['Nlam'])

         define_coords = """
         #---------------------
         # steps in each plane
         param dx := 1/(2*N);
         param dy := dx;
         
         param dmx := Rmask/M;
         param dmy := dmx;
         
         param dxi := rho2/Nimg;
         param deta := dxi;
 
         #---------------------
         # coordinate vectors in each plane
         set Xs := setof {i in 0.5..N-0.5 by 1} i*dx;
         set Ys := setof {j in 0.5..N-0.5 by 1} j*dy;
         
         set MXs := setof {i in 0.5..M-0.5 by 1} i*dmx;
         set MYs := setof {j in 0.5..M-0.5 by 1} j*dmy;

         set Xis := setof {i in 0..Nimg-1 by 1} i*dxi;
         set Etas := setof {j in 0..Nimg-1 by 1} j*deta;
         """
         
         load_masks = """\
         #---------------------
         # Load telescope aperture
         param TelAp {{x in Xs, y in Ys}};
         
         read {{y in Ys, x in Xs}} TelAp[x,y] < "{0:s}";
         close "{1:s}";
         
         # Load FPM
         param FPM {{mx in MXs, my in MYs}};
         
         read {{my in MYs, mx in MXs}} FPM[mx,my] < "{2:s}"; 
         close "{3:s}";
         
         # Load Lyot stop
         param LS {{x in Xs, y in Ys}};
         
         read {{y in Ys,x in Xs}} LS[x,y] < "{4:s}";
         close "{5:s}";
         """.format(self.fileorg['TelAp fname'], self.fileorg['TelAp fname'], self.fileorg['FPM fname'], self.fileorg['FPM fname'], \
                    self.fileorg['LS fname'], self.fileorg['LS fname'])
 

         if self.design['Image']['Nlam'] > 1 and self.design['Image']['bw'] > 0:
             define_wavelengths = """
             set Ls := setof {l in 1..Nlam} lam0*(1 - bw/2 + (l-1)*bw/(Nlam-1));
             """
         else:
             define_wavelengths = """
             set Ls := setof {l in 1..1} lam0*l;
             """

         sets_and_arrays = """
         #---------------------

         set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] != 0.} (x,y);
         set Mask := setof {mx in MXs, my in MYs: FPM[mx,my] != 0.} (mx,my);
         set Lyot := setof {x in Xs, y in Ys: LS[x,y] != 0.} (x,y);

         param TR := sum {(x,y) in Pupil} TelAp[x,y]*dx*dy; # Transmission of the Pupil. Used for calibration.
         param I00 := (sum {(x,y) in Pupil} TelAp[x,y]*LS[x,y]*dx*dy)^2; # Peak intensity in the absence of coronagraph
         param A_all {{x in Xs, y in Ys}};
         
         var A {(x,y) in Pupil} >= 0, <= 1, := 0.5;
         
         #---------------------
         
         set DarkHole := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= rho0 && sqrt(xi^2+eta^2) <= rho1} (xi,eta);
         """
 
         field_propagation = """
         #---------------------
         var EBm_real_X {mx in MXs, y in Ys, lam in Ls};
         var EBm_real {mx in MXs, my in MYs, lam in Ls};
         
         subject to st_EBm_real_X {mx in MXs, y in Ys, lam in Ls}: EBm_real_X[mx,y,lam] = 2.*sum {x in Xs: (x,y) in Pupil} A[x,y]*TelAp[x,y]*cos(2.*pi*x*mx*(lam0/lam))*dx;
         subject to st_EBm_real {(mx, my) in Mask, lam in Ls}: EBm_real[mx,my,lam] = 2.*(lam0/lam)*sum {y in Ys} EBm_real_X[mx,y,lam]*cos(2.*pi*y*my*(lam0/lam))*dy;
         
         #---------------------
         var ECm_real_X {x in Xs, my in MYs, lam in Ls};
         var ECm_real {x in Xs, y in Ys, lam in Ls};
         
         subject to st_ECm_real_X {x in Xs, my in MYs, lam in Ls}: ECm_real_X[x,my,lam] = 2.*sum {mx in MXs: (mx,my) in Mask} EBm_real[mx,my,lam]*cos(2.*pi*x*mx*(lam0/lam))*dmx;
         subject to st_ECm_real {(x,y) in Lyot, lam in Ls}: ECm_real[x,y,lam] = 2.*(lam0/lam)*sum {my in MYs} ECm_real_X[x,my,lam]*cos(2.*pi*y*my*(lam0/lam))*dmy;
         
         #---------------------
         var ED_real_X {xi in Xis, y in Ys, lam in Ls};
         var ED_real {xi in Xis, eta in Etas, lam in Ls};
         
         subject to st_ED_real_X {xi in Xis, y in Ys, lam in Ls}: ED_real_X[xi,y,lam] = 2.*sum {x in Xs: (x,y) in Lyot} (A[x,y]*TelAp[x,y]-ECm_real[x,y,lam])*cos(2.*pi*x*xi*(lam0/lam))*dx;
         subject to st_ED_real {(xi, eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] = 2.*(lam0/lam)*sum {y in Ys} ED_real_X[xi,y,lam]*cos(2.*pi*y*eta*(lam0/lam))*dy;
         
         #---------------------
         
         var ED00_real := 0.0;
         subject to st_ED00_real: ED00_real = 4.*sum {x in Xs, y in Ys: (x,y) in Lyot} (A[x,y]*TelAp[x,y])*dx*dy;
         """
 
         constraints = """
         #---------------------
         maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;
         
         subject to sidelobe_zero_real_pos {(xi,eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] <= 10^(-c/2)*ED00_real/sqrt(2.); 
         subject to sidelobe_zero_real_neg {(xi,eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] >= -10^(-c/2)*ED00_real/sqrt(2.);
         """
 
         misc_options = """
         option times 1;
         option gentimes 1;
         option show_stats 1;
         """
         
         solver = """
         option solver gurobi;
         """

         if self.solver['method'] is 'barhom':
             if self.solver['presolve'] is True:
                 solver_options = """
                 option gurobi_options "outlev=1 lpmethod=2 barhomogeneous=1 crossover=0";
                 """
             else:
                 solver_options = """
                 option gurobi_options "outlev=1 lpmethod=2 barhomogeneous=1 crossover=0 presolve=0";
                 """
         elif self.solver['method'] is 'bar':
             if self.solver['presolve'] is True:
                 solver_options = """
                 option gurobi_options "outlev=1 lpmethod=2 crossover=0";
                 """
             else:
                 solver_options = """
                 option gurobi_options "outlev=1 lpmethod=2 crossover=0 presolve=0";
                 """
         else: # assume dual simplex
             if self.solver['presolve'] is True:
                solver_options = """
                option gurobi_options "outlev=1 lpmethod=1";
                """
             else:
                solver_options = """
                option gurobi_options "outlev=1 lpmethod=1 presolve=0";
                """
  
         execute = """
         solve;
 
         display solve_result_num, solve_result;
         display ED00_real; 
         """
 
         store_results = """
         #---------------------

         param A_fin {{x in Xs, y in Ys}};
         let {{x in Xs, y in Ys}} A_fin[x,y] := 0;
         let {{(x,y) in Pupil}} A_fin[x,y] := A[x,y];
 
         printf {{x in Xs, y in Ys}}: "%15g %15g %15g \\n", x, y, A_fin[x,y] > "{0:s}";
         """.format(self.fileorg['sol fname'])
 
         mod_fobj.write( textwrap.dedent(header) )
         mod_fobj.write( textwrap.dedent(params) )
         mod_fobj.write( textwrap.dedent(define_coords) )
         mod_fobj.write( textwrap.dedent(load_masks) )
         mod_fobj.write( textwrap.dedent(define_wavelengths) )
         mod_fobj.write( textwrap.dedent(sets_and_arrays) )
         mod_fobj.write( textwrap.dedent(field_propagation) )
         mod_fobj.write( textwrap.dedent(constraints) )
         #mod_fobj.write( textwrap.dedent(misc_options) )
         mod_fobj.write( textwrap.dedent(solver) )
         mod_fobj.write( textwrap.dedent(solver_options) )
         mod_fobj.write( textwrap.dedent(execute) )
         mod_fobj.write( textwrap.dedent(store_results) )
 
         mod_fobj.close()
         if verbose:
             logging.info("Wrote %s"%self.fileorg['ampl src fname'])
