"""
Core definitions of python tools for the STScI
Segmented Coronagraph Design & Analysis investigation

02/14/2016 -- created by NTZ
"""

import os
import shutil
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
from astropy.io import ascii

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

class WrappedFixedIndentingLog(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None, style='%', width=70, indent=4):
        super(WrappedFixedIndentingLog, self).__init__(fmt=fmt, datefmt=datefmt)
        self.wrapper = textwrap.TextWrapper(width=width)
        #self.wrapper = textwrap.TextWrapper(width=width, subsequent_indent=' '*indent)
    def format(self, record):
        return self.wrapper.fill(super().format(record))

def make_ampl_bundle(coron_list, bundled_dir, queue_spec='12h'):
    bundled_coron_list = []
    if not os.path.exists(bundled_dir):
        os.makedirs(bundled_dir)
        
    cwd = os.getcwd()
    os.chdir(bundled_dir)
    bundled_fileorg = {'work dir': "."}
    
    bash_fname = "run_" + os.path.basename(os.path.normpath(bundled_dir)) + ".sh"
    bash_fobj = open(bash_fname, "w")
    bash_fobj.write("#! /bin/bash -x\n")
    
    for coron in coron_list:
        shutil.copy2(coron.fileorg['TelAp fname'], ".")
        shutil.copy2(coron.fileorg['FPM fname'], ".")
        shutil.copy2(coron.fileorg['LS fname'], ".")
        if 'LDZ fname' in coron.fileorg and coron.fileorg['LDZ fname'] is not None:
            shutil.copy2(coron.fileorg['LDZ fname'], ".")
        design_params = coron.design.copy()
        design_params['Image'].pop('Nimg',None)
        design_params['LS'].pop('s',None)
        bundled_coron = coron.__class__(design=coron.design, fileorg=bundled_fileorg,
                                        solver=coron.solver)
        bundled_coron_list.append(bundled_coron)
        if bundled_coron.check_ampl_input_files() is True:
            bundled_coron.write_ampl(overwrite=True)
            bundled_coron.write_exec_script(queue_spec, overwrite=True, verbose=False)
        else:
            scda.logging.warning("Input file configuration check failed; AMPL source file not written")
            scda.logging.warning("Bundled file organization: {0}".format(bundled_coron.fileorg))
        bash_fobj.write("ampl {0:s}\n".format(bundled_coron.fileorg['ampl src fname']))
    bash_fobj.close()
    os.chmod(bash_fname, 0775)
    os.chdir(cwd)
    return bundled_coron_list
    
def load_design_param_survey(pkl_fname):
    fobj = open(pkl_fname, 'rb')
    survey_obj = pickle.load(fobj)
    fobj.close() 
    return survey_obj

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
                                               "design initialization argument, expecting {3}").format(values, param, keycat, self._param_menu[keycat][param][0]) 
                                    logging.warning(warnstr)
                            else:
                                if isinstance(values, self._param_menu[keycat][param][0]):
                                    self.survey_config[keycat][param] = values
                                else:
                                    warnstr = ("Warning: Invalid {0} for parameter \"{1}\" under category \"{2}\" " + \
                                               "design initialization argument, expecting a {3}").format(type(values), param, keycat, self._param_menu[keycat][param][0]) 
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
                            self.fileorg[namekey] = os.path.expanduser(location)
                            if not os.path.exists(self.fileorg[namekey]):
                                logging.warning("Warning: The specified location of '{0}', \"{1}\" does not exist".format(namekey, self.fileorg[namekey]))
                        else:
                            self.fileorg[namekey] = location
                    else:
                        self.fileorg[namekey] = None
                else:
                    logging.warning("Warning: Unrecognized field {0} in fileorg argument".format(namekey))
        # Handle missing directory values, and create the directories if they don't exist
        if 'work dir' not in self.fileorg or self.fileorg['work dir'] is None:
            self.fileorg['work dir'] = os.getcwd()
        for namekey in self._file_fields['fileorg']: # Set other missing directory locations to 'work dir'
            if namekey.endswith('dir'):
                if namekey not in self.fileorg or self.fileorg[namekey] is None:
                    self.fileorg[namekey] = self.fileorg['work dir']
                if not os.path.exists(self.fileorg[namekey]):
                    os.mkdir(self.fileorg[namekey])

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

    def write_exec_script_batch(self, queue_spec='auto', overwrite=False, override_infile_status=False):
        for coron in self.coron_list:
            coron.write_exec_script(queue_spec, overwrite, verbose=False)
        logging.info("Wrote the batch of execution scripts into {:s}".format(self.fileorg['exec script dir']))
 
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
            if coron.eval_metrics['airy thrupt'] is None or coron.eval_metrics['fwhm area'] is None \
               or coron.eval_metrics['apod nb res ratio'] is None:
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
                paramrow.extend(['', 'filename', 'exists?', 'input files?', 'filename', 'exists?', 'Thrupt', 'PSF area'])
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
    _file_fields = { 'fileorg': ['work dir', 'ampl src dir', 'TelAp dir', 'FPM dir', 'LS dir',
                                 'sol dir', 'log dir', 'eval dir', 'exec script dir',
                                 'ampl src fname', 'exec script fname', 'log fname', 'job name', 
                                 'TelAp fname', 'FPM fname', 'LS fname', 'LDZ fname', 'sol fname'],
                     'solver': ['constr', 'method', 'presolve', 'threads', 'solver'] }

    _solver_menu = { 'constr': ['lin', 'quad'], 'solver': ['LOQO', 'gurobi', 'gurobix'], 
                     'method': ['bar', 'barhom', 'dualsimp'],
                     'presolve': [True, False], 'threads': [None]+range(1,33) }

    _aperture_menu = { 'prim': ['hex1', 'hex2', 'hex3', 'key24', 'pie12', 'pie08', 'irisao', 'atlast'],
                       'secobs': ['Y60d','Yoff60d','X','Cross','T','Y90d'],
                       'thick': ['025','100'],
                       'centobs': [True, False] }

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
                            self.fileorg[namekey] = os.path.expanduser(location)
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

        # Make directories for optimization solutions and logs if they don't exist 
        if not os.path.exists(self.fileorg['sol dir']):
            os.mkdir(self.fileorg['sol dir'])
        if not os.path.exists(self.fileorg['log dir']):
            os.mkdir(self.fileorg['log dir'])
       
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
        if 'LDZ fname' in self.fileorg and self.fileorg['LDZ fname'] is not None and \
        not os.path.exists(self.fileorg['LDZ fname']) and os.path.exists(self.fileorg['LS dir']) and \
        os.path.dirname(self.fileorg['LDZ fname']) == '':
            try_fname = os.path.join(self.fileorg['LS dir'], self.fileorg['LDZ fname']) 
            if os.path.exists(try_fname):
                self.fileorg['LDZ fname'] = try_fname
            else:
                logging.warning("Warning: Could not find the specified LDZ file \"{0}\" in {1}".format(self.fileorg['LDZ fname'], \
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
                elif field is 'convtol':
                    if 8 <= value < 10:
                        self.solver['convtol'] = value
                    else:
                        self.solver['convtol'] = None
                else:
                    logging.warning("Warning: Unrecognized field {0} in solver argument".format(field))
        # Handle missing values
        if 'constr' not in self.solver or self.solver['constr'] is None: self.solver['constr'] = 'lin'
        if 'solver' not in self.solver or self.solver['solver'] is None: self.solver['solver'] = 'gurobi'
        if 'method' not in self.solver or self.solver['method'] is None: self.solver['method'] = 'bar'
        if 'convtol' not in self.solver: self.solver['convtol'] = None
        if 'threads' not in self.solver: self.solver['threads'] = None
        if 'presolve' not in self.solver or self.solver['presolve'] is None: self.solver['presolve'] = True

        setattr(self, 'ampl_infile_status', None)
        if not issubclass(self.__class__, LyotCoronagraph):
            self.check_ampl_input_files()

        setattr(self, 'eval_metrics', {})
        self.eval_metrics['airy thrupt'] = None
        self.eval_metrics['fwhm area'] = None
        self.eval_metrics['apod nb res ratio'] = None

    def check_ampl_input_files(self):
        status = True
        if self.design['LS']['aligntol'] is not None:
            checklist = ['TelAp fname', 'FPM fname', 'LS fname', 'LDZ fname']
        else:
            checklist = ['TelAp fname', 'FPM fname', 'LS fname']
        for fname in checklist:
            if not os.path.exists(self.fileorg[fname]):
                status = False
                break
        self.ampl_infile_status = status
        return status

class NdiayeAPLC(LyotCoronagraph): # Image-constrained APLC following N'Diaye et al. (2015, 2016)
    _design_fields = OrderedDict([ ( 'Pupil', OrderedDict([('N',(int, 250)), ('prim',(str, 'hex1')), ('secobs',(str, 'x')), 
                                                           ('thick',(str, '025')), ('centobs',(bool, True))]) ),
                                   ( 'FPM', OrderedDict([('rad',(float, 4.)), ('M',(int, 50))]) ),
                                   ( 'LS', OrderedDict([('shape',(str, 'ann')), ('id',(int, 20)), ('od',(int, 90)), ('obscure',(int, 0)),
                                                        ('spad',(int, 0)), ('ppad',(int, 0)), ('aligntol',(int, None)), ('aligntolcon',(float, 7.))]) ),
                                   ( 'Image', OrderedDict([('c',(float, 10.)), ('ida',(float, -0.5)), ('oda',(float, 10.)),
                                                           ('bw',(float, 0.1)), ('Nlam',(int, 1)), ('fpres',(int,2)),
                                                           ('wingang',(float, None)), ('incon',(float, None)), ('wingcon',(float, None))]) ) ])
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
        # Set number of wavelength samples based on bandwidth
        if self.design['Image']['Nlam'] is None:
            self.design['Image']['Nlam'] = int(np.round(self.design['Image']['bw']/(0.10/3)))
        # Finally, set a private attribute for the number of image plane samples between the center and the outer constraint angle
        if self.design['Image']['wingang'] is not None:
            self.design['Image']['Nimg'] = int( np.ceil( self.design['Image']['fpres']*self.design['Image']['wingang']/(1. - self.design['Image']['bw']/2) ) )
        else:
            self.design['Image']['Nimg'] = int( np.ceil( self.design['Image']['fpres']*self.design['Image']['oda']/(1. - self.design['Image']['bw']/2) ) )
        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon']:
            # The Lyot dark zone field suppression factor decreases with the square of pupil array size. The units of input parameter are arbitrarily normalized to N=125.
            self.design['LS']['s'] = self.design['LS']['aligntolcon'] - 2*np.log10(self.design['Pupil']['N']/125.)
        if verbose: # Print summary of the set parameters
            logging.info("Design parameters: {}".format(self.design))
            logging.info("Optimization and solver parameters: {}".format(self.solver))
            logging.info("File organization parameters: {}".format(self.fileorg))
     
        self.amplname_coron = "APLC_full"
        self.amplname_pupil = "{0:s}{1:s}{2:s}cobs{3:d}_N{4:04d}".format(self.design['Pupil']['prim'], self.design['Pupil']['secobs'], self.design['Pupil']['thick'], \
                                                                         int(self.design['Pupil']['centobs']), self.design['Pupil']['N'])

        self.amplname_fpm = "FPM{:02}M{:03}".format(int(round(100*self.design['FPM']['rad'])), self.design['FPM']['M'])
        if self.design['LS']['obscure'] == 2: # LS includes primary and secondary aperture features
            self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}{3:s}Pad{4:02d}{5:s}{6:s}cobs{7:d}Pad{8:02d}".format(self.design['LS']['shape'], self.design['LS']['id'], \
                               self.design['LS']['od'], self.design['Pupil']['prim'], self.design['LS']['ppad'], self.design['Pupil']['secobs'], self.design['Pupil']['thick'], \
                               int(self.design['Pupil']['centobs']), self.design['LS']['spad'])
        elif self.design['LS']['obscure'] == 1: # LS includes secondary aperture features
            self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}{3:s}{4:s}cobs{5:d}Pad{6:02d}".format(self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'], \
                               self.design['Pupil']['secobs'], self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['spad'])
        else: # LS aperture is unobscured
            self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}clear".format(self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'])
        if self.design['LS']['aligntol'] is not None:
            self.amplname_ls += "Tol{0:02d}s{1:02d}".format(self.design['LS']['aligntol'], int(round(10*self.design['LS']['aligntolcon'])))

        if self.design['Image']['wingang'] is not None:
            self.amplname_image = "Img{:03}C_{:02}DA{:03}W{:03}_BW{:02}Nlam{:02}fpres{:1}".format(int(round(10*self.design['Image']['c'])), \
                                   int(round(10*(self.design['FPM']['rad']+self.design['Image']['ida']))), int(round(10*self.design['Image']['oda'])), \
                                   int(round(10*self.design['Image']['wingang'])), int(round(100*self.design['Image']['bw'])), \
                                   self.design['Image']['Nlam'], self.design['Image']['fpres'])
        else:
            self.amplname_image = "Img{:03}C_{:02}DA{:03}_BW{:02}Nlam{:02}fpres{:1}".format(int(round(10*self.design['Image']['c'])), \
                                   int(round(10*(self.design['FPM']['rad']+self.design['Image']['ida']))), int(round(10*self.design['Image']['oda'])), \
                                   int(round(100*self.design['Image']['bw'])), self.design['Image']['Nlam'], self.design['Image']['fpres'])
        if self.design['Image']['incon'] is not None:
            self.amplname_image += "Cin{:03}".format(self.design['Image']['incon'])
        if self.design['Image']['wingcon'] is not None:
            self.amplname_image += "Cw{:03}".format(self.design['Image']['wingcon'])

        if self.solver['presolve']:
            self.amplname_solver = "{}{}pre1".format(self.solver['constr'], self.solver['method'])
        else:
            self.amplname_solver = "{}{}pre0".format(self.solver['constr'], self.solver['method'])
        if self.solver['convtol'] is not None:
            self.amplname_solver += "convtol{0:2d}".format(int(round(10*self.solver['convtol'])))
        if self.solver['threads'] is not None:
            self.amplname_solver += "thr{:02d}".format(self.solver['threads'])

        if not issubclass(self.__class__, NdiayeAPLC): # Only set these file names if this is a full-plane APLC.
            if 'ampl src fname' not in self.fileorg or self.fileorg['ampl src fname'] is None:
                ampl_src_fname_tail = self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                                      self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".mod"
                self.fileorg['ampl src fname'] = os.path.join(self.fileorg['ampl src dir'], ampl_src_fname_tail)

            if 'sol fname' not in self.fileorg or self.fileorg['sol fname'] is None:
                sol_fname_tail = "ApodSol_" + self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                                 self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".dat"
                self.fileorg['sol fname'] = os.path.join(self.fileorg['sol dir'], sol_fname_tail)

            if 'TelAp fname' not in self.fileorg or self.fileorg['TelAp fname'] is None:
                self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_full_" + self.amplname_pupil + ".dat") )

            if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
                self.fileorg['FPM fname'] = os.path.join( self.fileorg['FPM dir'], "FPM_full_occspot_M{:03}.dat".format(self.design['FPM']['M']) )

            if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
                if self.design['LS']['obscure'] == 2:
                    self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_{3:s}Pad{4:02d}{5:s}{6:s}cobs{7:d}Pad{8:02d}_N{9:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['Pupil']['prim'], self.design['LS']['ppad'], self.design['Pupil']['secobs'],
                                                             self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['spad'],
                                                             self.design['Pupil']['N'])) )
                elif self.design['LS']['obscure'] == 1:
                    self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}cobs{5:d}Pad{6:02d}_N{7:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['Pupil']['secobs'], self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']),
                                                             self.design['LS']['spad'], self.design['Pupil']['N'])) )
                else:
                    self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_clear_N{3:04d}.dat".format(self.design['LS']['shape'],
                                                             self.design['LS']['id'], self.design['LS']['od'], self.design['Pupil']['N'])) )

            if self.design['LS']['aligntol'] is not None and ('LDZ fname' not in self.fileorg or self.fileorg['LDZ fname'] is None):
                if self.design['LS']['obscure'] == 2:
                    self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_{3:s}Pad{4:02d}{5:s}{6:s}c{7:d}Pad{8:02d}_Tol{9:02d}_N{10:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['Pupil']['prim'], self.design['LS']['ppad'], self.design['Pupil']['secobs'],
                                                             self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['spad'],
                                                             self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
                elif self.design['LS']['obscure'] == 1:
                    self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}c{5:d}Pad{6:02d}_Tol{7:02d}_N{8:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['Pupil']['secobs'], self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']),
                                                             self.design['LS']['spad'], self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
                else:
                    self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_clear_Tol{3:02d}_N{4:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['LS']['aligntol'], self.design['Pupil']['N'])) )

            self.check_ampl_input_files()
                                                              
    def write_ampl(self, overwrite=False):
        logging.info("Writing the AMPL program")
    def eval_onax_psf(self, fp2res=8, rho_inc=0.25, rho_out=None, Nlam=None):
        TelAp_qp = np.loadtxt(self.fileorg['TelAp fname'])
        TelAp = np.concatenate((np.concatenate((TelAp_qp[::-1,::-1], TelAp_qp[:,::-1]),axis=0),
                                np.concatenate((TelAp_qp[::-1,:], TelAp_qp),axis=0)), axis=1)

        FPM_qp = np.loadtxt(self.fileorg['FPM fname'])
        FPM = np.concatenate((np.concatenate((FPM_qp[::-1,::-1], FPM_qp[:,::-1]),axis=0),
                              np.concatenate((FPM_qp[::-1,:], FPM_qp),axis=0)), axis=1)
        
        LS_qp = np.loadtxt(self.fileorg['LS fname'])
        LS = np.concatenate((np.concatenate((LS_qp[::-1,::-1], LS_qp[:,::-1]),axis=0),
                             np.concatenate((LS_qp[::-1,:], LS_qp),axis=0)), axis=1)
        
        A_col = np.loadtxt(self.fileorg['sol fname'])[:,-1]
        A_qp = A_col.reshape(TelAp_qp.shape)
        A = np.concatenate((np.concatenate((A_qp[::-1,::-1], A_qp[:,::-1]),axis=0),
                            np.concatenate((A_qp[::-1,:], A_qp),axis=0)), axis=1)

        D = 1.
        N = self.design['Pupil']['N']
        bw = self.design['Image']['bw']
        M_fp1 = self.design['FPM']['M']
        fpm_rad = self.design['FPM']['rad']
        if rho_out is None:
            rho_out = self.design['Image']['oda'] + 1.
        if Nlam is None:
            Nlam = self.design['Image']['Nlam']
        M_fp2 = int(np.ceil(rho_out*fp2res))
        
        # pupil plane
        dx = (D/2)/N
        dy = dx
        xs = np.matrix(np.linspace(-N+0.5,N-0.5,2*N)*dx)
        ys = xs
        
        # FPM
        dmx = fpm_rad/M_fp1
        dmy = dmx
        mxs = np.matrix(np.linspace(-M_fp1+0.5,M_fp1-0.5,2*M_fp1)*dmx)
        mys = mxs
        
        # FP2
        dxi = 1./fp2res
        xis = np.matrix(np.linspace(-M_fp2+0.5,M_fp2-0.5,2*M_fp2)*dxi)
        etas = xis
        
        # wavelength ratios
        wrs = np.linspace(1.-bw/2, 1.+bw/2, Nlam)

        intens_polychrom = np.zeros((Nlam, 2*M_fp2, 2*M_fp2))
        for wi, wr in enumerate(wrs):
            Psi_B = dx*dx/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(mxs.T, xs)), TelAp*A ),
                                           np.exp(-1j*2*np.pi/wr*np.dot(xs.T, mxs)))
            Psi_B_stop = np.multiply(Psi_B, FPM)
            Psi_C = A*TelAp - dmx*dmx/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(xs.T, mxs)), Psi_B_stop),
                                                       np.exp(-1j*2*np.pi/wr*np.dot(mxs.T, xs)))
            Psi_C_stop = np.multiply(Psi_C, LS)
            Psi_D = dx*dx/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(xis.T, xs)), Psi_C_stop),
                                           np.exp(-1j*2*np.pi/wr*np.dot(xs.T, xis)))
            Psi_D_0_peak = np.sum(A*TelAp*LS)*dx*dx/wr
            intens_polychrom[wi,:,:] = np.power(np.absolute(Psi_D)/Psi_D_0_peak, 2)
             
        seps = np.arange(self.design['FPM']['rad']+self.design['Image']['ida'], rho_out, rho_inc)
        radial_intens_polychrom = np.zeros((len(wrs), len(seps)))
        XXs = np.asarray(np.dot(np.matrix(np.ones(xis.shape)).T, xis))
        YYs = np.asarray(np.dot(etas.T, np.matrix(np.ones(etas.shape))))
        RRs = np.sqrt(XXs**2 + YYs**2)

        for si, sep in enumerate(seps):
            r_in = np.max([seps[0], sep-0.5])
            r_out = np.min([seps[-1], sep+0.5])
            meas_ann_mask = np.logical_and(np.greater_equal(RRs, r_in),
                                           np.less_equal(RRs, r_out))
            meas_ann_ind = np.nonzero(np.logical_and(np.greater_equal(RRs, r_in).ravel(),
                                                     np.less_equal(RRs, r_out).ravel()))[0]
            for wi, wr in enumerate(wrs):
                radial_intens_polychrom[wi, si] = np.mean(np.ravel(intens_polychrom[wi,:,:])[meas_ann_ind])

        return intens_polychrom, seps, radial_intens_polychrom
    
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
            sol_fname_tail = "ApodSol_" + self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                             self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".dat"
            self.fileorg['sol fname'] = os.path.join(self.fileorg['sol dir'], sol_fname_tail)

        if 'TelAp fname' not in self.fileorg or self.fileorg['TelAp fname'] is None:
            self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_half_" + self.amplname_pupil + ".dat") )

        if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
            self.fileorg['FPM fname'] = os.path.join( self.fileorg['FPM dir'], "FPM_half_occspot_M{:03}.dat".format(self.design['FPM']['M']) )

        if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}Pad{4:02d}{5:s}{6:s}cobs{7:d}Pad{8:02d}_N{9:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['prim'], self.design['LS']['ppad'], self.design['Pupil']['secobs'],
                                                         self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['spad'],
                                                         self.design['Pupil']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}cobs{5:d}Pad{6:02d}_N{7:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['secobs'], self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']),
                                                         self.design['LS']['spad'], self.design['Pupil']['N'])) )
            else:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_clear_N{3:04d}.dat".format(self.design['LS']['shape'],
                                                         self.design['LS']['id'], self.design['LS']['od'], self.design['Pupil']['N'])) )

        if self.design['LS']['aligntol'] is not None and ('LDZ fname' not in self.fileorg or self.fileorg['LDZ fname'] is None):
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}Pad{4:02d}{5:s}{6:s}c{7:d}Pad{8:02d}_Tol{9:02d}_N{10:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['prim'], self.design['LS']['ppad'], self.design['Pupil']['secobs'],
                                                         self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['spad'],
                                                         self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}c{5:d}Pad{6:02d}_Tol{7:02d}_N{8:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['secobs'], self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']),
                                                         self.design['LS']['spad'], self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
            else:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_clear_Tol{3:02d}_N{4:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['LS']['aligntol'], self.design['Pupil']['N'])) )

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
                self.fileorg['ampl src fname'] = ampl_src_fname
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
        
        #---------------------
        param bw := {12:0.2f};
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
                   self.design['Image']['c'], self.design['FPM']['rad'], self.design['FPM']['rad']+self.design['Image']['ida'],
                   self.design['Image']['oda'], self.design['Pupil']['N'], self.design['FPM']['M'], self.design['Image']['Nimg'], \
                   self.design['Image']['bw'], self.design['Image']['Nlam'])
        
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
        
        param dxi := (rho1/Nimg)*(1/CoeffOverSizePup);
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

        if 'job name' not in self.fileorg or self.fileorg['job name'] is None:
            self.fileorg['job name'] = os.path.basename(self.fileorg['ampl src fname'])[:-4]
 
        if 'sol fname' not in self.fileorg or self.fileorg['sol fname'] is None:
            sol_fname_tail = "ApodSol_" + self.fileorg['job name'] + ".dat"
            self.fileorg['sol fname'] = os.path.join(self.fileorg['sol dir'], sol_fname_tail)

        if 'exex script fname' not in self.fileorg or self.fileorg['exec script fname'] is None:
            exec_script_fname_tail = self.fileorg['job name'] + ".sh"
            self.fileorg['exec script fname'] = os.path.join(self.fileorg['exec script dir'], exec_script_fname_tail)

        if 'log fname' not in self.fileorg or self.fileorg['log fname'] is None:
            log_fname_tail = self.fileorg['job name'] + ".log"
            self.fileorg['log fname'] = os.path.join(self.fileorg['log dir'], log_fname_tail)
 
        if 'TelAp fname' not in self.fileorg or self.fileorg['TelAp fname'] is None:
            self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_quart_" + self.amplname_pupil + ".dat") )
 
        if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
            self.fileorg['FPM fname'] = os.path.join( self.fileorg['FPM dir'], "FPM_quart_occspot_M{:03}.dat".format(self.design['FPM']['M']) )

        if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_quart_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}Pad{4:02d}{5:s}{6:s}cobs{7:d}Pad{8:02d}_N{9:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['prim'], self.design['LS']['ppad'], self.design['Pupil']['secobs'],
                                                         self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['spad'],
                                                         self.design['Pupil']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_quart_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}cobs{5:d}Pad{6:02d}_N{7:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['secobs'], self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']),
                                                         self.design['LS']['spad'], self.design['Pupil']['N'])) )
            else:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_quart_" + \
                                                         "{0:s}{1:02d}D{2:02d}_clear_N{3:04d}.dat".format(self.design['LS']['shape'],
                                                         self.design['LS']['id'], self.design['LS']['od'], self.design['Pupil']['N'])) )

        if self.design['LS']['aligntol'] is not None and ('LDZ fname' not in self.fileorg or self.fileorg['LDZ fname'] is None):
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_quart_" + \
                                                          "{0:s}{1:02d}D{2:02d}_{3:s}Pad{4:02d}{5:s}{6:s}cobs{7:d}Pad{8:02d}_Tol{9:02d}_N{10:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['Pupil']['prim'], self.design['LS']['ppad'], self.design['Pupil']['secobs'],
                                                          self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['spad'],
                                                          self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_quart_" + \
                                                          "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}cobs{5:d}Pad{6:02d}_Tol{7:02d}_N{8:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['Pupil']['secobs'], self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']),
                                                          self.design['LS']['spad'], self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
            else:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_quart_" + \
                                                          "{0:s}{1:02d}D{2:02d}_clear_Tol{3:02d}_N{4:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
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
                self.fileorg['ampl src fname'] = ampl_src_fname
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
        # load amplgsl.dll;
        """.format(getpass.getuser(), os.path.basename(__file__), socket.gethostname(), datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
 
        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon'] is not None: 
            params = """
            #---------------------
 
            param pi:= 4*atan(1);
 
            #---------------------
            param c := {0:.2f};
            param s := {1:.2f};
 
            #---------------------
            param Rmask := {2:0.3f};
            param rho0 := {3:0.2f};
            param rho1 := {4:0.2f};
            
            #---------------------
            param N := {5:d};				# discretization parameter (pupil)
            param M := {6:d};				# discretization parameter (mask)
            param Nimg := {7:d};           # discretization parameter (image)
                                  
            #---------------------
            param bw := {8:0.2f};
            param Nlam := {9:d};
            
            #---------------------
            """.format(self.design['Image']['c'], self.design['LS']['s'], self.design['FPM']['rad'], \
                       self.design['FPM']['rad']+self.design['Image']['ida'], self.design['Image']['oda'], self.design['Pupil']['N'], \
                       self.design['FPM']['M'], self.design['Image']['Nimg'], self.design['Image']['bw'], self.design['Image']['Nlam'])
        else:
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
                                  
            #---------------------
            param bw := {7:0.2f};
            param Nlam := {8:d};
            
            #---------------------
            """.format(self.design['Image']['c'], self.design['FPM']['rad'], self.design['FPM']['rad']+self.design['Image']['ida'],
                       self.design['Image']['oda'], self.design['Pupil']['N'], self.design['FPM']['M'], self.design['Image']['Nimg'], \
                       self.design['Image']['bw'], self.design['Image']['Nlam'])

        define_coords = """
        #---------------------
        # steps in each plane
        param dx := 1/(2*N);
        param dy := dx;
        
        param dmx := Rmask/M;
        param dmy := dmx;
        
        param dxi := rho1/Nimg;
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
       
        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon'] is not None: 
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
            
            # Load Lyot dark zone
            param LDZ {{x in Xs, y in Ys}};
            read {{y in Ys,x in Xs}} LDZ[x,y] < "{6:s}";
            close "{7:s}";
            """.format(self.fileorg['TelAp fname'], self.fileorg['TelAp fname'], self.fileorg['FPM fname'], self.fileorg['FPM fname'], \
                       self.fileorg['LS fname'], self.fileorg['LS fname'], self.fileorg['LDZ fname'], self.fileorg['LDZ fname'])
        else:
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
            set Ls := setof {l in 1..Nlam} 1 - bw/2 + (l-1)*bw/(Nlam-1);
            """
        else:
            define_wavelengths = """
            set Ls := setof {l in 1..1} 1;
            """

        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon'] is not None: 
            sets_and_arrays = """
            #---------------------

            set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] > 0} (x,y);
            set Mask := setof {mx in MXs, my in MYs: FPM[mx,my] > 0} (mx,my);
            set Lyot := setof {x in Xs, y in Ys: LS[x,y] > 0} (x,y);
            set LyotDarkZone := setof {x in Xs, y in Ys: LDZ[x,y] > 0} (x,y);

            param TR := sum {(x,y) in Pupil} TelAp[x,y]*dx*dy; # Transmission of the Pupil. Used for calibration.
            
            var A {x in Xs, y in Ys} >= 0, <= 1, := 0.5;
            
            #---------------------
            set DarkHole := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= rho0 && sqrt(xi^2+eta^2) <= rho1} (xi,eta);
            """
        else:
            sets_and_arrays = """
            #---------------------

            set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] > 0} (x,y);
            set Mask := setof {mx in MXs, my in MYs: FPM[mx,my] > 0} (mx,my);
            set Lyot := setof {x in Xs, y in Ys: LS[x,y] > 0} (x,y);

            param TR := sum {(x,y) in Pupil} TelAp[x,y]*dx*dy; # Transmission of the Pupil. Used for calibration.
            
            var A {x in Xs, y in Ys} >= 0, <= 1, := 0.5;
            
            #---------------------
            set DarkHole := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= rho0 && sqrt(xi^2+eta^2) <= rho1} (xi,eta);
            """
 
        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon'] is not None: 
            field_propagation = """
            #---------------------
            var EBm_real_X {mx in MXs, y in Ys, lam in Ls};
            var EBm_real {mx in MXs, my in MYs, lam in Ls};
            
            subject to st_EBm_real_X {mx in MXs, y in Ys, lam in Ls}: EBm_real_X[mx,y,lam] = 2*sum {x in Xs: (x,y) in Pupil} TelAp[x,y]*A[x,y]*cos(2*pi*x*mx/lam)*dx;
            subject to st_EBm_real {(mx, my) in Mask, lam in Ls}: EBm_real[mx,my,lam] = 2/lam*sum {y in Ys} EBm_real_X[mx,y,lam]*cos(2*pi*y*my/lam)*dy;
            
            #---------------------
            var ECm_real_X {x in Xs, my in MYs, lam in Ls};
            var EC_real {x in Xs, y in Ys, lam in Ls};
            
            subject to st_ECm_real_X {x in Xs, my in MYs, lam in Ls}: ECm_real_X[x,my,lam] = 2*sum {mx in MXs: (mx,my) in Mask} EBm_real[mx,my,lam]*cos(2*pi*x*mx/lam))*dmx;
            subject to st_EC_real {(x,y) in Lyot union LyotDarkZone, lam in Ls}: EC_real[x,y,lam] = TelAp[x,y]*A[x,y] - 2/lam*sum {my in MYs} ECm_real_X[x,my,lam]*cos(2*pi*y*my/lam)*dmy;
            
            #---------------------
            var ED_real_X {xi in Xis, y in Ys, lam in Ls};
            var ED_real {xi in Xis, eta in Etas, lam in Ls};
            
            subject to st_ED_real_X {xi in Xis, y in Ys, lam in Ls}: ED_real_X[xi,y,lam] = 2*sum {x in Xs: (x,y) in Lyot} EC_real[x,y,lam]*cos(2*pi*x*xi/lam)*dx;
            subject to st_ED_real {(xi, eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] = 2/lam*sum {y in Ys} ED_real_X[xi,y,lam]*cos(2*pi*y*eta/lam)*dy;
            
            #---------------------
            var ED00_real := 0.0;
            subject to st_ED00_real: ED00_real = 4.*sum {x in Xs, y in Ys: (x,y) in Lyot} (A[x,y]*TelAp[x,y])*dx*dy;
            """
        else:
            field_propagation = """
            #---------------------
            var EBm_real_X {mx in MXs, y in Ys, lam in Ls};
            var EBm_real {mx in MXs, my in MYs, lam in Ls};
            
            subject to st_EBm_real_X {mx in MXs, y in Ys, lam in Ls}: EBm_real_X[mx,y,lam] = 2*sum {x in Xs: (x,y) in Pupil} TelAp[x,y]*A[x,y]*cos(2*pi*x*mx/lam)*dx;
            subject to st_EBm_real {(mx, my) in Mask, lam in Ls}: EBm_real[mx,my,lam] = 2/lam*sum {y in Ys} EBm_real_X[mx,y,lam]*cos(2*pi*y*my/lam)*dy;
            
            #---------------------
            var ECm_real_X {x in Xs, my in MYs, lam in Ls};
            var ECm_real {x in Xs, y in Ys, lam in Ls};
            
            subject to st_ECm_real_X {x in Xs, my in MYs, lam in Ls}: ECm_real_X[x,my,lam] = 2*sum {mx in MXs: (mx,my) in Mask} EBm_real[mx,my,lam]*cos(2*pi*x*mx/lam)*dmx;
            subject to st_ECm_real {(x,y) in Lyot, lam in Ls}: ECm_real[x,y,lam] = 2/lam*sum {my in MYs} ECm_real_X[x,my,lam]*cos(2*pi*y*my/lam)*dmy;
            
            #---------------------
            var ED_real_X {xi in Xis, y in Ys, lam in Ls};
            var ED_real {xi in Xis, eta in Etas, lam in Ls};
            
            subject to st_ED_real_X {xi in Xis, y in Ys, lam in Ls}: ED_real_X[xi,y,lam] = 2*sum {x in Xs: (x,y) in Lyot} (TelAp[x,y]*A[x,y]-ECm_real[x,y,lam])*cos(2*pi*x*xi/lam)*dx;
            subject to st_ED_real {(xi, eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] = 2/lam*sum {y in Ys} ED_real_X[xi,y,lam]*cos(2*pi*y*eta/lam)*dy;
            
            #---------------------
            var ED00_real := 0.0;
            subject to st_ED00_real: ED00_real = 4*sum {x in Xs, y in Ys: (x,y) in Lyot} (A[x,y]*TelAp[x,y])*dx*dy;
            """

        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon'] is not None:
            constraints = """
            #---------------------
            maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;
           
            subject to Lyot_aligntol_constr_pos {(x,y) in LyotDarkZone, lam in Ls}: EC_real[x,y,lam] <= 10^-s;
            subject to Lyot_aligntol_constr_neg {(x,y) in LyotDarkZone, lam in Ls}: EC_real[x,y,lam] >= -10^-s;
            subject to sidelobe_zero_real_pos {(xi,eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] <= 10^(-c/2)*ED00_real/lam/sqrt(2.);
            subject to sidelobe_zero_real_neg {(xi,eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] >= -10^(-c/2)*ED00_real/lam/sqrt(2.);
            """
        else:
            constraints = """
            #---------------------
            maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;
            
            subject to sidelobe_zero_real_pos {(xi,eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] <= 10^(-c/2)*ED00_real/lam/sqrt(2.);
            subject to sidelobe_zero_real_neg {(xi,eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] >= -10^(-c/2)*ED00_real/lam/sqrt(2.);
            """
 
        misc_options = """
        option times 1;
        option gentimes 1;
        option show_stats 1;
        """
        
        solver = """
        option solver gurobi;
        """

        gurobi_opt_str = "outlev=1"
        if self.solver['presolve'] is False:
            gurobi_opt_str += " presolve=0"
        if self.solver['method'] is 'bar' or self.solver['method'] is 'barhom':
            gurobi_opt_str += " lpmethod=2 crossover=0"
            if self.solver['convtol'] is not None:
                gurobi_opt_str += " barconvtol={0:.1e}".format(np.power(10,-self.solver['convtol']))
            if self.solver['method'] is 'barhom':
                gurobi_opt_str += " barhomogeneous=1"
        else: # assume dual simplex
            gurobi_opt_str += " lpmethod=1"

        solver_options = """
        option gurobi_options "{0:s}";
        """.format(gurobi_opt_str)

        execute = """
        solve;
 
        display solve_result_num, solve_result;
        display ED00_real; 
        """
 
        store_results = """
        #---------------------

        param A_fin {{y in Ys, x in Xs}};
        let {{y in Ys, x in Xs}} A_fin[x,y] := 0;
        let {{(x,y) in Pupil}} A_fin[x,y] := A[x,y];
 
        printf {{y in Ys, x in Xs}}: "%15g %15g %15g \\n", x, y, A_fin[x,y] > "{0:s}";
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

    def write_exec_script(self, queue_spec='auto', overwrite=False, verbose=True):
        if os.path.exists(self.fileorg['exec script fname']):
            if overwrite == True:
                if verbose:
                    logging.warning("Warning: Overwriting the existing copy of {0}".format(self.fileorg['exec script fname']))
            else:
                logging.warning("Error: {0} already exists and overwrite switch is off, so write_exec_script() will now abort".format(self.fileorg['exec script fname']))
                return
        elif not os.path.exists(self.fileorg['exec script dir']):
            os.mkdir(self.fileorg['exec script dir'])
            if verbose:
                logging.info("Created new execution script directory, {0:s}".format(self.fileorg['exec script dir']))

        bash_fobj = open(self.fileorg['exec script fname'], "w") 

        begin_script = """\
        #!/bin/bash

        #PBS -V
        #PBS -m e -M {0:s}@stsci.edu
        #SBATCH --job-name={1:s}
        #SBATCH -o {2:s}
        #SBATCH --account=s1649
        
        #SBATCH --constraint=hasw
        #SBATCH --ntasks=1 --nodes=1
        """.format(getpass.getuser(), self.fileorg['job name'], self.fileorg['log fname'])

        if queue_spec is 'auto':
            if self.design['LS']['aligntol'] is None:
                time_est_hrs = int(np.ceil(1*(self.design['Pupil']['N']/125.)**2*(self.design['Image']['Nlam']/3.)**3))
            else:
                time_est_hrs = int(np.ceil(3*(self.design['Pupil']['N']/125.)**2*(self.design['Image']['Nlam']/3.)**3))
            if time_est_hrs > 12:
                set_queue = """
                #SBATCH --qos=long
                #SBATCH --time={0:02d}:00:00
                """.format(time_est_hrs)
            else:
                set_queue = """
                #SBATCH --qos=allnccs
                #SBATCH --time={0:02d}:00:00
                """.format(time_est_hrs)
        elif queue_spec is '24h':
            set_queue = """
            #SBATCH --qos=long
            #SBATCH --time=24:00:00
            """
        elif queue_spec is '12h':
            set_queue = """
            #SBATCH --qos=allnccs
            #SBATCH --time=12:00:00
            """

        intel_module = """
        . /usr/share/modules/init/bash
        module purge
        module load comp/intel-10.1.017
        ulimit -s unlimited
        """

        monitor_mem = """
        #Optional: monitor the memory usage...
        mkdir -p ${NOBACKUP}/policeme
        /usr/local/other/policeme/policeme.exe -d ${NOBACKUP}/policeme
        """

        call_ampl = """
        ampl {0:s}
        
        exit 0""".format(self.fileorg['ampl src fname'])

        bash_fobj.write( textwrap.dedent(begin_script) )
        bash_fobj.write( textwrap.dedent(set_queue) )
        bash_fobj.write( textwrap.dedent(intel_module) )
        bash_fobj.write( textwrap.dedent(monitor_mem) )
        bash_fobj.write( textwrap.dedent(call_ampl) )

        bash_fobj.close()
        if verbose:
            logging.info("Wrote %s"%self.fileorg['exec script fname'])

    def get_metrics(self, fp2res=16, verbose=True):
        TelAp = np.loadtxt(self.fileorg['TelAp fname'])
        FPM = np.loadtxt(self.fileorg['FPM fname'])
        LS = np.loadtxt(self.fileorg['LS fname'])
        N = TelAp.shape[1]
        A_col = np.loadtxt(self.fileorg['sol fname'])[:,-1]
        A = A_col.reshape(TelAp.shape)
        self.eval_metrics['apod nb res ratio'] = np.sum(np.abs(A - np.round(A)))/np.sum(TelAp)
        dx = (1./2)/N
        dy = dx
        xs = np.matrix(np.linspace(0.5,N-0.5,N)*dx)
        ys = xs.copy()
        rho1 = 3.
        M_fp2 = int(np.ceil(rho1*fp2res))
        dxi = 1./fp2res
        xis = np.matrix(np.linspace(0.5,M_fp2-0.5,M_fp2)*dxi)
        etas = xis.copy()
        wrs = np.linspace(1.-self.design['Image']['bw']/2, 1.+self.design['Image']['bw']/2, self.design['Image']['Nlam'])
        airy_thrupt_polychrom = []
        fwhm_area_polychrom = []
        for wr in wrs:
            Psi_D_0 = 4*dx*dy/wr*np.dot(np.cos(2*np.pi/wr*np.dot(etas.T, ys)), A*TelAp*LS).dot(np.cos(2*np.pi/wr*np.dot(xs.T, xis)))
            Intens_D_0 = np.power(np.absolute(Psi_D_0), 2)
            Intens_D_0_peak = (4*np.sum(TelAp*A*LS)*dx*dy/wr)**2
            fwhm_ind_APLC = np.greater_equal(Intens_D_0, Intens_D_0_peak/2)
            Psi_TelAp = 4*dx*dy/wr*np.dot(np.cos(2*np.pi/wr*np.dot(etas.T, ys)), TelAp).dot(np.cos(2*np.pi/wr*np.dot(xs.T, xis)))
            Intens_TelAp = np.power(np.absolute(Psi_TelAp), 2)
            Intens_TelAp_peak = (4*np.sum(TelAp)*dx*dy/wr)**2
            fwhm_ind_TelAp = np.greater_equal(Intens_TelAp, Intens_TelAp_peak/2)
            fwhm_sum_TelAp = np.sum(Intens_TelAp[fwhm_ind_TelAp])
            fwhm_sum_APLC = np.sum(Intens_D_0[fwhm_ind_APLC])
            fwhm_area_polychrom.append(4.*np.sum(fwhm_ind_APLC)/(fp2res**2))
            airy_thrupt_polychrom.append(fwhm_sum_APLC/fwhm_sum_TelAp)
        self.eval_metrics['airy thrupt'] = np.mean(airy_thrupt_polychrom)
        self.eval_metrics['fwhm area'] = np.mean(fwhm_area_polychrom)
        if verbose:
            print("Non-binary residuals, as a percentage of clear telescope aperture area: {:.2f}%".format(100*self.eval_metrics['apod nb res ratio']))
            print("Band-averaged Airy throughput: {:.2f}%".format(100*self.eval_metrics['airy thrupt']))
            print("Band-averaged FWHM PSF area / (lambda0/D)^2: {:.2f}".format(self.eval_metrics['fwhm area']))
