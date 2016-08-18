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
import itertools
import pprint
import pickle
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

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

def make_ampl_bundle(coron_list, bundled_dir, queue_spec='auto', email=None, arch=None):
    bundled_coron_list = []
    if not os.path.exists(bundled_dir):
        os.makedirs(bundled_dir)
        
    cwd = os.getcwd()
    os.chdir(bundled_dir)
    
    serial_bash_fname = "run_" + os.path.basename(os.path.normpath(bundled_dir)) + "_serial.sh"
    serial_bash_fobj = open(serial_bash_fname, "w")
    serial_bash_fobj.write("#! /bin/bash -x\n")
    sbatch_bash_fname = "run_" + os.path.basename(os.path.normpath(bundled_dir)) + "_sbatch.sh"
    sbatch_bash_fobj = open(sbatch_bash_fname, "w")
    sbatch_bash_fobj.write("#! /bin/bash -x\n")
    
    for coron in coron_list:
        bundled_fileorg = {'work dir': ".", 'ampl src fname': os.path.basename(coron.fileorg['ampl src fname']),
                           'TelAp fname': os.path.basename(coron.fileorg['TelAp fname']), 
                           'FPM fname': os.path.basename(coron.fileorg['FPM fname']), 
                           'LS fname': os.path.basename(coron.fileorg['LS fname'])} 
        if not os.path.exists(os.path.basename(coron.fileorg['TelAp fname'])):
            shutil.copy2(coron.fileorg['TelAp fname'], ".")
        if not os.path.exists(os.path.basename(coron.fileorg['FPM fname'])):
            shutil.copy2(coron.fileorg['FPM fname'], ".")
        if not os.path.exists(os.path.basename(coron.fileorg['LS fname'])):
            shutil.copy2(coron.fileorg['LS fname'], ".")
        if 'LDZ fname' in coron.fileorg and coron.fileorg['LDZ fname'] is not None:
            if not os.path.exists(os.path.basename(coron.fileorg['LDZ fname'])):
                shutil.copy2(coron.fileorg['LDZ fname'], ".")
            bundled_fileorg['LDZ fname'] = os.path.basename(coron.fileorg['LDZ fname'])
        design_params = coron.design.copy()
        #design_params['FPM'].pop('M',None)
        design_params['LS'].pop('s',None)
        design_params['Image'].pop('Nimg',None)
        design_params['Image'].pop('bw+',None)
        bundled_coron = coron.__class__(design=coron.design, fileorg=bundled_fileorg,
                                        solver=coron.solver)
        bundled_coron_list.append(bundled_coron)
        if bundled_coron.check_ampl_input_files() is True:
            bundled_coron.write_ampl(overwrite=True)
            bundled_coron.write_slurm_script(queue_spec=queue_spec, email=email, arch=arch, 
                                            overwrite=True, verbose=False)
        else:
            logging.warning("Input file configuration check failed; AMPL source file not written")
            logging.warning("Bundled file organization: {0}".format(bundled_coron.fileorg))
        serial_bash_fobj.write("ampl {0:s}\n".format(bundled_coron.fileorg['ampl src fname']))
        sbatch_bash_fobj.write("sbatch {0:s}\n".format(bundled_coron.fileorg['slurm fname']))
    serial_bash_fobj.close()
    sbatch_bash_fobj.close()
    os.chmod(serial_bash_fname, 0775)
    os.chmod(sbatch_bash_fname, 0775)
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
        self._solver_menu = coron_class._solver_menu.copy()
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
        if 'crossover' not in self.solver: self.solver['crossover'] = None
         
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
 
        setattr(self, 'ampl_infile_status', False)
        self.check_ampl_input_files()
        setattr(self, 'ampl_src_status', False)
        setattr(self, 'ampl_submission_status', False)
        setattr(self, 'solution_status', False)
        setattr(self, 'eval_status', False)

    def write_ampl_batch(self, overwrite=False, override_infile_status=False):
        write_count = 0
        overwrite_deny_count = 0
        infile_deny_count = 0
        for coron in self.coron_list:
            status = coron.write_ampl(overwrite, override_infile_status, verbose=False)
            if status == 2:
                infine_deny_count += 1
            elif status == 1:
                overwrite_deny_count += 1
            else:
                write_count += 1
        if write_count == self.N_combos:
            logging.info("Wrote all {0:d} of {1:d} design survey AMPL programs into {2:s}".format(write_count, self.N_combos, self.fileorg['ampl src dir']))
        else:
            logging.warning("Wrote {0:d} of {1:d} design survey AMPL programs into {2:s}. {3:d} already existed and were denied overwriting. {4:d} were denied writing because of a failed input file configuration status.".format(write_count, self.N_combos, self.fileorg['ampl src dir'], overwrite_deny_count, infile_deny_count))

    def write_slurm_batch(self, queue_spec='auto', account='s1649', email=None, arch=None,
                          overwrite=False, override_infile_status=False):
        write_count = 0
        overwrite_deny_count = 0
        for coron in self.coron_list:
            status = coron.write_slurm_script(queue_spec=queue_spec, account=account, email=email, arch=arch,
                                              overwrite=overwrite, verbose=False)
            if status == 1:
                overwrite_deny_count += 1
            else:
                write_count += 1
        if write_count == self.N_combos:
            logging.info("Wrote all {0:d} of {1:d} design survey slurm scripts into {2:s}".format(write_count, self.N_combos, self.fileorg['slurm dir']))
        else:
            logging.warning("Wrote {0:d} of {1:d} design survey AMPL programs into {2:s}. {3:d} already existed and were denied overwriting.".format(write_count, self.N_combos, self.fileorg['slurm dir'], overwrite_deny_count))

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
            if coron.eval_metrics['fwhm area'] is None or coron.eval_metrics['apod nb res ratio'] is None:
                status = False
                break
        self.eval_status = status
        return status

    def get_metrics(self, fp2res=16, verbose=False):
        for coron in self.coron_list:
            if os.path.exists(coron.fileorg['sol fname']) and \
                (coron.eval_metrics['fwhm area'] is None \
                 or coron.eval_metrics['apod nb res ratio'] is None):
                coron.get_metrics(verbose=verbose)
                coron.eval_status = True

    def write(self, fname=None):
        if fname is not None:
            if os.path.dirname(fname) is '': # if no path specified, assume work dir
                self.fileorg['survey fname'] = os.path.join(self.fileorg['work dir'], fname)
            else:
                self.fileorg['survey fname'] = fname
        else:
            if 'survey fname' not in self.fileorg or \
               ('survey fname' in self.fileorg and self.fileorg['survey fname'] is None): # set the filename based on the coronagraph type, user, and date
                fname_tail = "{0:s}_{1:s}_{2:s}.pkl".format(os.path.basename(os.path.abspath(self.fileorg['work dir'])), getpass.getuser(), datetime.datetime.now().strftime("%Y-%m-%d"))
                self.fileorg['survey fname'] = os.path.join(self.fileorg['work dir'], fname_tail)
        fobj = open(self.fileorg['survey fname'], 'wb')
        pickle.dump(self, fobj)
        fobj.close()
        os.chmod(self.fileorg['survey fname'], 0644)
        logging.info("Wrote the design parameter survey object to {:s}".format(self.fileorg['survey fname']))
 
    def write_spreadsheet(self, overwrite=False, csv_fname=None):
        if csv_fname is not None:
            if os.path.dirname(csv_fname) is '': # if no path specified, assume work dir
                csv_fname = os.path.join(self.fileorg['work dir'], csv_fname)
            else:
                csv_fname = csv_fname
        else:
            if 'survey fname' not in self.fileorg or ('survey fname' in self.fileorg and self.fileorg['survey fname'] is None):
                #csv_fname_tail = "scda_{:s}_survey_{:s}_{:s}.csv".format(self.coron_class.__name__, getpass.getuser(), datetime.datetime.now().strftime("%Y-%m-%d"))
                csv_fname_tail = "{0:s}_{1:s}_{2:s}.csv".format(os.path.basename(os.path.abspath(self.fileorg['work dir'])), getpass.getuser(), datetime.datetime.now().strftime("%Y-%m-%d"))
                csv_fname = os.path.join(self.fileorg['work dir'], csv_fname_tail)
            else:
                csv_fname = self.fileorg['survey fname'][:-4] + ".csv"
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
            if self.ampl_submission_status is True:
                surveywriter.writerow(["All AMPL jobs submitted?", 'Y'])
            else:
                surveywriter.writerow(["All AMPL jobs submitted?", 'N'])
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
                catrow.extend(['', 'AMPL program', '', '', '', 'Solution', '', 'Evaluation metrics', '', ''])
                paramrow.extend(['', 'filename', 'exists?', 'input files?', 'submitted?', 'filename', 'exists?',
                                 'inc. energy', 'apodizer non-binarity', 'Tot thrupt', 'half-max thrupt', 'r=0.7 thrupt', 'rel. half-max thrupt', 'rel. r=0.7 thrupt', 'PSF area'])
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
                    if self.coron_list[ii].ampl_submission_status is True:
                        param_combo_row.append('Y')
                    else:
                        param_combo_row.append('N')
                    param_combo_row.append(os.path.basename(self.coron_list[ii].fileorg['sol fname'])) 
                    if os.path.exists(self.coron_list[ii].fileorg['sol fname']):
                        param_combo_row.append('Y')
                    else:
                        param_combo_row.append('N')
                    if self.coron_list[ii].eval_metrics['inc energy'] is not None:
                        param_combo_row.append(self.coron_list[ii].eval_metrics['inc energy'])
                    else:
                        param_combo_row.append('')
                    if self.coron_list[ii].eval_metrics['apod nb res ratio'] is not None:
                        param_combo_row.append(self.coron_list[ii].eval_metrics['apod nb res ratio'])
                    else:
                        param_combo_row.append('')
                    if self.coron_list[ii].eval_metrics['tot thrupt'] is not None:
                        param_combo_row.append(self.coron_list[ii].eval_metrics['tot thrupt'])
                    else:
                        param_combo_row.append('')
                    if self.coron_list[ii].eval_metrics['fwhm thrupt'] is not None:
                        param_combo_row.append(self.coron_list[ii].eval_metrics['fwhm thrupt'])
                    else:
                        param_combo_row.append('')
                    if self.coron_list[ii].eval_metrics['p7ap thrupt'] is not None:
                        param_combo_row.append(self.coron_list[ii].eval_metrics['p7ap thrupt'])
                    else:
                        param_combo_row.append('')
                    if self.coron_list[ii].eval_metrics['rel fwhm thrupt'] is not None:
                        param_combo_row.append(self.coron_list[ii].eval_metrics['rel fwhm thrupt'])
                    else:
                        param_combo_row.append('')
                    if self.coron_list[ii].eval_metrics['rel p7ap thrupt'] is not None:
                        param_combo_row.append(self.coron_list[ii].eval_metrics['rel p7ap thrupt'])
                    else:
                        param_combo_row.append('')
                    if self.coron_list[ii].eval_metrics['fwhm area'] is not None:
                        param_combo_row.append(self.coron_list[ii].eval_metrics['fwhm area'])
                    else:
                        param_combo_row.append('')
                    surveywriter.writerow(param_combo_row)
                    
        survey_spreadsheet.close()
        os.chmod(csv_fname, 0644)
        logging.info("Wrote design survey spreadsheet to {:s}".format(csv_fname))

class LyotCoronagraph(object): # Lyot coronagraph base class
    _file_fields = { 'fileorg': ['work dir', 'ampl src dir', 'TelAp dir', 'FPM dir', 'LS dir',
                                 'sol dir', 'log dir', 'eval dir', 'slurm dir',
                                 'ampl src fname', 'slurm fname', 'log fname', 'job name', 
                                 'TelAp fname', 'FPM fname', 'LS fname', 'LDZ fname', 'sol fname'],
                     'solver': ['constr', 'method', 'presolve', 'threads', 'solver', 'crossover'] }

    _solver_menu = { 'constr': ['lin', 'quad'], 'solver': ['LOQO', 'gurobi', 'gurobix'], 
                     'method': ['bar', 'barhom', 'dualsimp'],
                     'presolve': [True, False], 'threads': [None]+range(1,33), 'crossover': [None]+[True, False] }

    _aperture_menu = { 'prim': ['hex1', 'hex2', 'hex3', 'hex4', 'key24', 'pie12', 'pie08', 'irisao', 'atlast'],
                       'secobs': ['Y60d','Yoff60d','X','Cross','T','Y90d'],
                       'thick': ['025','100'],
                       'centobs': [True, False],
                       'edge': ['gray', 'round', 'floor'] }

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
        if 'crossover' not in self.solver: self.solver['crossover'] = None

        setattr(self, 'ampl_infile_status', None)
        if not issubclass(self.__class__, LyotCoronagraph):
            self.check_ampl_input_files()

        setattr(self, 'ampl_submission_status', None) # Only changed by the queue filler program

        setattr(self, 'eval_metrics', {})
        self.eval_metrics['inc energy'] = None
        self.eval_metrics['tot thrupt'] = None
        self.eval_metrics['fwhm thrupt'] = None
        self.eval_metrics['p7ap thrupt'] = None
        self.eval_metrics['rel fwhm thrupt'] = None
        self.eval_metrics['rel p7ap thrupt'] = None
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
                logging.warning("Missing {:s}".format(self.fileorg[fname]))
                break
        self.ampl_infile_status = status
        return status

class SPLC(LyotCoronagraph): # SPLC following Zimmerman et al. (2016), uses diaphragm FPM 
    _design_fields = OrderedDict([ ( 'Pupil', OrderedDict([('N',(int, 125)), ('prim',(str, 'hex3')), ('secobs',(str, 'X')), 
                                                           ('thick',(str, '025')), ('centobs',(bool, True)),
                                                           ('gap',(int, 1)), ('edge',(str, 'gray'))]) ),
                                   ( 'FPM', OrderedDict([('R0',(float, 4.)), ('R1',(float, 10.)), ('openang',(int, 180)),
                                                         ('orient',(str, 'H')), ('fpmres',(int, 10))]) ),
                                   ( 'LS', OrderedDict([('N',(int, 125)), ('shape',(str, 'ann')), ('id',(int, 25)), ('od',(int, 75)),
                                                        ('obscure',(int, 0)), ('pad',(int, 0)),
                                                        ('aligntol',(int, None)), ('aligntolcon',(float, 3.))]) ),
                                   ( 'Image', OrderedDict([('c',(float, 10.)), ('bw',(float, 0.10)), ('Nlam',(int, 5)),
                                                           ('dR',(float, -0.5)), ('fpres',(int,2))]) ) ])
    _eval_fields =   { 'Pupil': _design_fields['Pupil'], 'FPM': _design_fields['FPM'], \
                       'LS': _design_fields['LS'], 'Image': _design_fields['Image'], \
                       'Tel': {'TelAp diam':(float, 12.)}, 'Target': {}, 'Aber': {}, 'WFSC': {} }

    def __init__(self, verbose=False, **kwargs):
        super(SPLC, self).__init__(**kwargs)

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
        self.design['FPM']['M'] = int(np.ceil(self.design['FPM']['fpmres']*self.design['FPM']['R1']))
        # NOTE: temporarily removing bandpass augmentation in favor of direct image constraints down to R0+dR
        #if False and self.design['Image']['dR'] < 0: # expand optimization bandpass according to the focal plane stellar diam./pointing tolerance parameter, 'dR'
        if self.design['Image']['dR'] < 0: # expand optimization bandpass according to the focal plane stellar diam./pointing tolerance parameter, 'dR'
            self.design['Image']['bw+'] = self.design['Image']['bw']*self.design['FPM']['R0']/(self.design['FPM']['R0'] + self.design['Image']['dR'])
        else:
            self.design['Image']['bw+'] = self.design['Image']['bw']
        # Unless Nlam is explicitly specified, set the number of wavelength samples according to the bandwidth
        if self.design['Image']['Nlam'] == 1 and self.design['Image']['bw+'] > 0:
            self.design['Image']['Nlam'] = int(np.ceil(self.design['Image']['bw+']/(0.10/3)))
        # Set a private attribute for the number of image plane samples between the center and the outer constraint angle
        self.design['Image']['Nimg'] = int( np.ceil( self.design['Image']['fpres']*self.design['FPM']['R1']/(1. - self.design['Image']['bw+']/2) ) )
        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon']:
            # The Lyot dark zone field suppression factor decreases with the square of pupil array size. The units of input parameter are arbitrarily normalized to N=250.
            self.design['LS']['s'] = self.design['LS']['aligntolcon'] - 2*np.log10(self.design['Pupil']['N']/250.)

        if verbose: # Print summary of the set parameters
            logging.info("Design parameters: {}".format(self.design))
            logging.info("Optimization and solver parameters: {}".format(self.solver))
            logging.info("File organization parameters: {}".format(self.fileorg))
     
        self.amplname_coron = "SPLC_full"
        self.telap_descrip = "{0:s}{1:s}{2:s}cobs{3:d}gap{4:d}_N{5:04d}".format(self.design['Pupil']['prim'], self.design['Pupil']['secobs'], self.design['Pupil']['thick'], \
                                                                                int(self.design['Pupil']['centobs']), self.design['Pupil']['gap'], self.design['Pupil']['N'])
        self.amplname_pupil = "{0:s}{1:s}".format(self.telap_descrip, self.design['Pupil']['edge'][0])

        self.amplname_fpm = "FPM{0:02d}R{1:03d}{2:s}{3:03d}res{4:02d}".format(int(round(10*self.design['FPM']['R0'])), int(round(10*self.design['FPM']['R1'])),
                                                                              self.design['FPM']['orient'], self.design['FPM']['openang'], self.design['FPM']['fpmres'])
        if self.design['LS']['obscure'] == 2: # LS includes primary and secondary aperture features
            self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}".format(self.design['LS']['shape'], self.design['LS']['id'], \
                               self.design['LS']['od'], self.design['Pupil']['prim'], self.design['Pupil']['secobs'], self.design['Pupil']['thick'], \
                               int(self.design['Pupil']['centobs']), self.design['LS']['pad'])
        elif self.design['LS']['obscure'] == 1: # LS includes secondary aperture features
            self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}{3:s}{4:s}cobs{5:d}Pad{6:02d}".format(self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'], \
                               self.design['Pupil']['secobs'], self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'])
        else: # LS aperture is unobscured
            self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}clear".format(self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'])
        if self.design['LS']['aligntol'] is not None:
            self.amplname_ls += "Tol{0:02d}s{1:02d}".format(self.design['LS']['aligntol'], int(round(10*self.design['LS']['aligntolcon'])))
        self.amplname_ls += "N{0:04d}".format(self.design['LS']['N'])

        self.amplname_image = "ImgC{0:03d}BW{1:02d}Nlam{2:02d}dR{3:1d}res{4:1d}".format(int(round(10*self.design['Image']['c'])), \
                               int(round(100*self.design['Image']['bw'])), self.design['Image']['Nlam'], \
                               int(round(-10*self.design['Image']['dR'])), self.design['Image']['fpres'])
        if self.solver['presolve']:
            self.amplname_solver = "{}{}pre1".format(self.solver['constr'], self.solver['method'])
        else:
            self.amplname_solver = "{}{}pre0".format(self.solver['constr'], self.solver['method'])
        if self.solver['convtol'] is not None:
            self.amplname_solver += "convtol{0:2d}".format(int(round(10*self.solver['convtol'])))
        if self.solver['crossover'] is not None:
            self.amplname_solver += "cross"
        if self.solver['threads'] is not None:
            self.amplname_solver += "thr{:02d}".format(self.solver['threads'])

        if not issubclass(self.__class__, SPLC): # Only set these file names if this is a full-plane SPLC.
            if 'ampl src fname' not in self.fileorg or self.fileorg['ampl src fname'] is None:
                ampl_src_fname_tail = self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                                      self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".mod"
                self.fileorg['ampl src fname'] = os.path.join(self.fileorg['ampl src dir'], ampl_src_fname_tail)

            if 'job name' not in self.fileorg or self.fileorg['job name'] is None:
                self.fileorg['job name'] = os.path.basename(self.fileorg['ampl src fname'])[:-4]

            if 'sol fname' not in self.fileorg or self.fileorg['sol fname'] is None:
                sol_fname_tail = "ApodSol_" + self.fileorg['job name'] + ".dat"
                self.fileorg['sol fname'] = os.path.join(self.fileorg['sol dir'], sol_fname_tail)

            if 'slurm fname' not in self.fileorg or self.fileorg['slurm fname'] is None:
                exec_script_fname_tail = self.fileorg['job name'] + ".sh"
                self.fileorg['slurm fname'] = os.path.join(self.fileorg['slurm dir'], exec_script_fname_tail)

            if 'log fname' not in self.fileorg or self.fileorg['log fname'] is None:
                log_fname_tail = self.fileorg['job name'] + ".log"
                self.fileorg['log fname'] = os.path.join(self.fileorg['log dir'], log_fname_tail)

            if 'TelAp fname' not in self.fileorg or self.fileorg['TelAp fname'] is None:
                self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_full_" + self.telap_descrip + ".dat") )

            if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
                self.fileorg['FPM fname'] = \
                  os.path.join( self.fileorg['FPM dir'], "FPM_full_diaphragm_{0:03d}M{1:03d}_{2:s}{3:03d}deg.dat".format(
                                                          int(round(self.design['FPM']['fpmres']*self.design['FPM']['R0'])),
                                                          int(np.ceil(self.design['FPM']['fpmres']*self.design['FPM']['R1'])),
                                                          self.design['FPM']['orient'], self.design['FPM']['openang']) )

            if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
                if self.design['LS']['obscure'] == 2:
                    self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}_N{8:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['Pupil']['prim'], self.design['Pupil']['secobs'],
                                                             self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'],
                                                             self.design['LS']['N'])) )
                elif self.design['LS']['obscure'] == 1:
                    self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}Pad{5:02d}_N{6:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['Pupil']['secobs'], self.design['Pupil']['thick'],
                                                             self.design['LS']['pad'], self.design['LS']['N'])) )
                else:
                    self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_clear_N{3:04d}.dat".format(self.design['LS']['shape'],
                                                             self.design['LS']['id'], self.design['LS']['od'], self.design['LS']['N'])) )

            if self.design['LS']['aligntol'] is not None and ('LDZ fname' not in self.fileorg or self.fileorg['LDZ fname'] is None):
                if self.design['LS']['obscure'] == 2:
                    self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}_Tol{8:02d}_N{9:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['Pupil']['prim'], self.design['Pupil']['secobs'],
                                                             self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'],
                                                             self.design['LS']['aligntol'], self.design['LS']['N'])) )
                elif self.design['LS']['obscure'] == 1:
                    self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}Pad{5:02d}_Tol{6:02d}_N{7:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['Pupil']['secobs'], self.design['Pupil']['thick'],
                                                             self.design['LS']['pad'], self.design['LS']['aligntol'], self.design['LS']['N'])) )
                else:
                    self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_clear_Tol{3:02d}_N{4:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['LS']['aligntol'], self.design['LS']['N'])) )
            self.check_ampl_input_files()

    def write_ampl(self, overwrite=False):
        logging.info("Writing the AMPL program") # Not yet written for full-plane SPLC

    def get_onax_psf(self, fp2res=8, rho_inc=0.25, rho_out=None, Nlam=None): # for SPLC
        if self.design['Pupil']['edge'] is 'floor': # floor to binary
            TelAp_p = np.floor(np.loadtxt(self.fileorg['TelAp fname'])).astype(int)
        elif self.design['Pupil']['edge'] is 'round': # round to binary
            TelAp_p = np.round(np.loadtxt(self.fileorg['TelAp fname'])).astype(int)
        else: # keey it gray
            TelAp_p = np.loadtxt(self.fileorg['TelAp fname'])
        A_col = np.loadtxt(self.fileorg['sol fname'])[:,-1]
        FPM_p = np.loadtxt(self.fileorg['FPM fname'])
        LS_p = np.loadtxt(self.fileorg['LS fname'])
        A_p = A_col.reshape(TelAp_p.shape)
        if isinstance(self, QuarterplaneSPLC):
            TelAp = np.concatenate((np.concatenate((TelAp_p[::-1,::-1], TelAp_p[:,::-1]),axis=0),
                                    np.concatenate((TelAp_p[::-1,:], TelAp_p),axis=0)), axis=1)
            A = np.concatenate((np.concatenate((A_p[::-1,::-1], A_p[:,::-1]),axis=0),
                                np.concatenate((A_p[::-1,:], A_p),axis=0)), axis=1)
            FPM = np.concatenate((np.concatenate((FPM_p[::-1,::-1], FPM_p[:,::-1]),axis=0),
                                  np.concatenate((FPM_p[::-1,:], FPM_p),axis=0)), axis=1)
            LS = np.concatenate((np.concatenate((LS_p[::-1,::-1], LS_p[:,::-1]),axis=0),
                                 np.concatenate((LS_p[::-1,:], LS_p),axis=0)), axis=1)
        elif isinstance(self, HalfplaneSPLC):
            TelAp = np.concatenate((TelAp_p[:,::-1], TelAp_p), axis=1)
            A = np.concatenate((A_p[:,::-1], A_p), axis=1)
            FPM = np.concatenate((np.concatenate((FPM_p[::-1,::-1], FPM_p[:,::-1]),axis=0),
                                  np.concatenate((FPM_p[::-1,:], FPM_p),axis=0)), axis=1)
            LS = np.concatenate((LS_p[:,::-1], LS_p), axis=1)
        else:
            TelAp = TelAp_p
            A = A_p
            FPM = FPM_p
            LS = LS_p
        
        D = 1.
        N_A = self.design['Pupil']['N']
        N_L = self.design['LS']['N']
        bw = self.design['Image']['bw']
        fpmres = self.design['FPM']['fpmres']
        M_fp1 = FPM_p.shape[0]
        if rho_out is None:
            rho_out = self.design['FPM']['R1'] + 1.
        if Nlam is None:
            Nlam = self.design['Image']['Nlam']
        M_fp2 = int(np.ceil(rho_out*fp2res))
        
        # pupil plane
        dx = (D/2)/N_A
        dy = dx
        xs = np.matrix(np.linspace(-N_A+0.5,N_A-0.5,2*N_A)*dx)
        ys = xs
        
        # FPM
        dmx = 1./fpmres
        dmy = dmx
        mxs = np.matrix(np.linspace(-M_fp1+0.5,M_fp1-0.5,2*M_fp1)*dmx)
        mys = mxs

        # Lyot plane
        du = (D/2)/N_L
        dv = du
        us = np.matrix(np.linspace(-N_L+0.5,N_L-0.5,2*N_L)*du)
        vs = us
        
        # FP2
        dxi = 1./fp2res
        xis = np.matrix(np.linspace(-M_fp2+0.5,M_fp2-0.5,2*M_fp2)*dxi)
        etas = xis
        
        # wavelength ratios
        wrs = np.linspace(1.-bw/2, 1.+bw/2, Nlam)

        intens_polychrom = np.zeros((Nlam, 2*M_fp2, 2*M_fp2))
        for wi, wr in enumerate(wrs):
            Psi_B = dx*dy/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(mxs.T, xs)), TelAp*A ),
                                           np.exp(-1j*2*np.pi/wr*np.dot(xs.T, mxs)))
            Psi_B_stop = np.multiply(Psi_B, FPM)
            Psi_C = dmx*dmx/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(us.T, mxs)), Psi_B_stop),
                                             np.exp(-1j*2*np.pi/wr*np.dot(mxs.T, us)))
            Psi_C_stop = np.multiply(Psi_C, LS)
            Psi_D = du*dv/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(xis.T, us)), Psi_C_stop),
                                           np.exp(-1j*2*np.pi/wr*np.dot(us.T, xis)))

            Psi_C_0 = dmx*dmx/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(us.T, mxs)), Psi_B),
                                               np.exp(-1j*2*np.pi/wr*np.dot(mxs.T, us)))
            Psi_C_0_stop = np.multiply(Psi_C_0, LS)
            Psi_D_0_peak = np.sum(Psi_C_0_stop)*du*dv/wr
            intens_polychrom[wi,:,:] = np.power(np.absolute(Psi_D)/np.absolute(Psi_D_0_peak), 2)
             
        #seps = np.arange(self.design['FPM']['R0']+self.design['Image']['dR'], rho_out, rho_inc)
        seps = np.arange(0, rho_out, rho_inc)
        radial_intens_polychrom = np.zeros((len(wrs), len(seps)))
        XXs = np.asarray(np.dot(np.matrix(np.ones(xis.shape)).T, xis))
        YYs = np.asarray(np.dot(etas.T, np.matrix(np.ones(etas.shape))))
        RRs = np.sqrt(XXs**2 + YYs**2)
        theta_quad = np.rad2deg(np.arctan2(YYs[M_fp2:,M_fp2:], XXs[M_fp2:,M_fp2:]))
        #ang_mask = np.ones((2*M_fp2, 2*M_fp2))
        rad_mask_1d = np.logical_and(np.greater_equal(RRs, self.design['FPM']['R0']).ravel(),
                                     np.less_equal(RRs, self.design['FPM']['R1']).ravel())
        rad_mask = rad_mask_1d.reshape(RRs.shape)
        if self.design['FPM']['openang'] < 180:
            if self.design['FPM']['orient'] == 'V':
                theta_quad_mask = np.greater_equal(theta_quad, self.design['FPM']['openang']/2)
            else:
                theta_quad_mask = np.less_equal(theta_quad, self.design['FPM']['openang']/2)
            theta_rhs_mask = np.concatenate((theta_quad_mask[::-1,:], theta_quad_mask), axis=0)
            theta_mask = np.concatenate((theta_rhs_mask[:,::-1], theta_rhs_mask), axis=1)
            FoV_mask = theta_mask*rad_mask
        else:
            FoV_mask = rad_mask

        for si, sep in enumerate(seps):
            r_in = np.max([seps[0], sep-0.5])
            r_out = np.min([seps[-1], sep+0.5])
            meas_ann_ind = np.nonzero(np.logical_and(np.greater_equal(RRs, r_in).ravel(),
                                                     np.less_equal(RRs, r_out).ravel(),
                                                     np.greater(theta_mask, 0).ravel()))[0]
            for wi, wr in enumerate(wrs):
                radial_intens_polychrom[wi, si] = np.mean(np.ravel(intens_polychrom[wi,:,:])[meas_ann_ind])

        #pdb.set_trace()

        return intens_polychrom, seps, radial_intens_polychrom, FoV_mask

    def get_metrics(self, fp1res=8, fp2res=16, rho_out=None, Nlam=None, verbose=True): # for SPLC class
        TelAp_basename = os.path.basename(self.fileorg['TelAp fname'])
        gapstr_beg = TelAp_basename.find('gap')
        TelAp_nopad_basename = TelAp_basename.replace(TelAp_basename[gapstr_beg:gapstr_beg+4], 'gap0')
        TelAp_nopad_fname = os.path.join( os.path.dirname(self.fileorg['TelAp fname']), TelAp_nopad_basename )
        TelAp_p = np.loadtxt(TelAp_nopad_fname)
        A_col = np.loadtxt(self.fileorg['sol fname'])[:,-1]
        LS_p = np.loadtxt(self.fileorg['LS fname'])
        A_p = A_col.reshape(TelAp_p.shape)
        if isinstance(self, QuarterplaneSPLC):
            TelAp = np.concatenate((np.concatenate((TelAp_p[::-1,::-1], TelAp_p[:,::-1]),axis=0),
                                    np.concatenate((TelAp_p[::-1,:], TelAp_p),axis=0)), axis=1)
            A = np.concatenate((np.concatenate((A_p[::-1,::-1], A_p[:,::-1]),axis=0),
                                np.concatenate((A_p[::-1,:], A_p),axis=0)), axis=1)
            LS = np.concatenate((np.concatenate((LS_p[::-1,::-1], LS_p[:,::-1]),axis=0),
                                 np.concatenate((LS_p[::-1,:], LS_p),axis=0)), axis=1)
        elif isinstance(self, HalfplaneSPLC):
            TelAp = np.concatenate((TelAp_p[:,::-1], TelAp_p), axis=1)
            A = np.concatenate((A_p[:,::-1], A_p), axis=1)
            LS = np.concatenate((LS_p[:,::-1], LS_p), axis=1)
        else:
            TelAp = TelAp_p
            A = A_p
            LS = LS_p
            
        D = 1.
        N_A = self.design['Pupil']['N']
        N_L = self.design['LS']['N']
        self.eval_metrics['apod nb res ratio'] = np.sum(np.abs(A - np.round(A)))/np.sum(TelAp)
        if Nlam is None:
            Nlam = self.design['Image']['Nlam']
        if rho_out is None:
            rho_out = self.design['Image']['oda']

        # Pupil plane
        dx = (D/2)/N_A
        dy = dx
        xs = np.matrix(np.linspace(-N_A+0.5, N_A-0.5, 2*N_A)*dx)
        ys = xs
        # FP1
        M_fp1 = int(np.ceil(rho_out*fp1res))
        dmx = 1./fp1res
        dmy = dmx
        mxs = np.matrix(np.linspace(-M_fp1+0.5, M_fp1-0.5, 2*M_fp1)*dmx)
        mys = mxs
        # Lyot plane
        du = (D/2)/N_L
        dv = du
        us = np.matrix(np.linspace(-N_L+0.5,N_L-0.5,2*N_L)*du)
        vs = us
        # FP2
        M_fp2 = int(np.ceil(rho_out*fp2res))
        dxi = 1./fp2res
        xis = np.matrix(np.linspace(-M_fp2+0.5,M_fp2-0.5,2*M_fp2)*dxi)
        etas = xis.copy()
        # Compute unocculted PSF
        wrs = np.linspace(1.-self.design['Image']['bw']/2, 1.+self.design['Image']['bw']/2, Nlam)
        XXs = np.asarray(np.dot(np.matrix(np.ones(xis.shape)).T, xis))
        YYs = np.asarray(np.dot(etas.T, np.matrix(np.ones(etas.shape))))
        RRs = np.sqrt(XXs**2 + YYs**2)
        p7ap_ind = np.less_equal(RRs, 0.7)

        intens_D_0_polychrom = np.zeros((Nlam, 2*M_fp2, 2*M_fp2))
        intens_D_0_peak_polychrom = np.zeros((Nlam, 1))
        intens_TelAp_polychrom = np.zeros((Nlam, 2*M_fp2, 2*M_fp2))
        intens_TelAp_peak_polychrom = np.zeros((Nlam, 1))
        for wi, wr in enumerate(wrs):
            Psi_B_0 = dx*dy/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(mxs.T, xs)), TelAp*A),
                                             np.exp(-1j*2*np.pi/wr*np.dot(xs.T, mxs)))
            Psi_C_0 = dmx*dmy/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(us.T, mxs)), Psi_B_0),
                                               np.exp(-1j*2*np.pi/wr*np.dot(mxs.T, us)))
            Psi_C_0_stop = np.multiply(Psi_C_0, LS)
            Psi_D_0 = du*dv/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(xis.T, us)), Psi_C_0_stop),
                                             np.exp(-1j*2*np.pi/wr*np.dot(us.T, xis)))
            Psi_D_0_peak = du*dv/wr*np.sum(Psi_C_0_stop)

            intens_D_0_polychrom[wi] = np.power(np.absolute(Psi_D_0), 2)
            intens_D_0_peak_polychrom[wi] = np.power(np.absolute(Psi_D_0_peak), 2)
            Psi_TelAp = dx*dy/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(xis.T, xs)), TelAp),
                                               np.exp(-1j*2*np.pi/wr*np.dot(xs.T, xis)))
            intens_TelAp_polychrom[wi] = np.power(np.absolute(Psi_TelAp), 2)
            intens_TelAp_peak_polychrom[wi] = (np.sum(TelAp)*dx*dy/wr)**2

        intens_D_0 = np.mean(intens_D_0_polychrom, axis=0)
        intens_D_0_peak = np.mean(intens_D_0_peak_polychrom)
        intens_TelAp = np.mean(intens_TelAp_polychrom, axis=0)
        intens_TelAp_peak = np.mean(intens_TelAp_peak_polychrom)

        fwhm_ind_SPLC = np.greater_equal(intens_D_0, intens_D_0_peak/2)
        fwhm_ind_TelAp = np.greater_equal(intens_TelAp, intens_TelAp_peak/2)
        fwhm_area_polychrom.append(np.sum(fwhm_ind_SPLC)*dxi*dxi)

        fwhm_sum_TelAp = np.sum(intens_TelAp[fwhm_ind_TelAp])*dxi*dxi
        fwhm_sum_SPLC = np.sum(intens_D_0[fwhm_ind_SPLC])*dxi*dxi
        p7ap_sum_TelAp = np.sum(intens_TelAp[p7ap_ind])*dxi*dxi
        p7ap_sum_SPLC = np.sum(intens_D_0[p7ap_ind])*dxi*dxi

        self.eval_metrics['inc energy'] = np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['tot thrupt'] = np.sum(intens_D_0*dxi*dxi)/np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['fwhm thrupt'] = fwhm_sum_SPLC/np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['p7ap thrupt'] = p7ap_sum_SPLC/np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['rel fwhm thrupt'] = fwhm_sum_SPLC/fwhm_sum_TelAp
        self.eval_metrics['rel p7ap thrupt'] = p7ap_sum_SPLC/p7ap_sum_TelAp
        self.eval_metrics['fwhm area'] = np.sum(fwhm_ind_SPLC)*dxi*dxi
        if verbose:
            print("////////////////////////////////////////////////////////")
            print("{:s}".format(self.fileorg['job name']))
            print("Incident energy on aperture (dimensionless): {:.3f}".format(self.eval_metrics['inc energy']))
            print("Non-binary residuals, as a percentage of clear telescope aperture area: {:.2f}%".format(100*self.eval_metrics['apod nb res ratio']))
            print("Band-averaged total throughput: {:.2f}%".format(100*self.eval_metrics['tot thrupt']))
            print("Band-averaged half-max throughput: {:.2f}%".format(100*self.eval_metrics['fwhm thrupt']))
            print("Band-averaged r=.7 lam/D throughput: {:.2f}%".format(100*self.eval_metrics['p7ap thrupt']))
            print("Band-averaged relative half-max throughput: {:.2f}%".format(100*self.eval_metrics['rel fwhm thrupt']))
            print("Band-averaged relative r=0.7 lam/D throughput: {:.2f}%".format(100*self.eval_metrics['rel p7ap thrupt']))
            print("Band-averaged FWHM PSF area / (lambda0/D)^2: {:.2f}".format(self.eval_metrics['fwhm area']))

class QuarterplaneSPLC(SPLC): # Zimmerman SPLC subclass for the quarter-plane symmetry case
    def __init__(self, **kwargs):
        super(QuarterplaneSPLC, self).__init__(**kwargs)
        self.amplname_coron = "SPLC_quart"
        if 'ampl src fname' not in self.fileorg or self.fileorg['ampl src fname'] is None:
            ampl_src_fname_tail = self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                                  self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".mod"
            self.fileorg['ampl src fname'] = os.path.join(self.fileorg['ampl src dir'], ampl_src_fname_tail)

        if 'job name' not in self.fileorg or self.fileorg['job name'] is None:
            self.fileorg['job name'] = os.path.basename(self.fileorg['ampl src fname'])[:-4]

        if 'sol fname' not in self.fileorg or self.fileorg['sol fname'] is None:
            sol_fname_tail = "ApodSol_" + self.fileorg['job name'] + ".dat"
            self.fileorg['sol fname'] = os.path.join(self.fileorg['sol dir'], sol_fname_tail)

        if 'slurm fname' not in self.fileorg or self.fileorg['slurm fname'] is None:
            exec_script_fname_tail = self.fileorg['job name'] + ".sh"
            self.fileorg['slurm fname'] = os.path.join(self.fileorg['slurm dir'], exec_script_fname_tail)

        if 'log fname' not in self.fileorg or self.fileorg['log fname'] is None:
            log_fname_tail = self.fileorg['job name'] + ".log"
            self.fileorg['log fname'] = os.path.join(self.fileorg['log dir'], log_fname_tail)

        if 'TelAp fname' not in self.fileorg or self.fileorg['TelAp fname'] is None:
            self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_quart_" + self.telap_descrip + ".dat") )

        if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
            self.fileorg['FPM fname'] = \
              os.path.join( self.fileorg['FPM dir'], "FPM_quart_diaphragm_{0:03d}M{1:03d}_{2:s}{3:03d}deg.dat".format(
                                                      int(round(self.design['FPM']['fpmres']*self.design['FPM']['R0'])),
                                                      int(np.ceil(self.design['FPM']['fpmres']*self.design['FPM']['R1'])),
                                                      self.design['FPM']['orient'], self.design['FPM']['openang']) )

        if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_quart_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}_N{8:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['prim'], self.design['Pupil']['secobs'],
                                                         self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'],
                                                         self.design['LS']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_quart_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}Pad{5:02d}_N{6:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['secobs'], self.design['Pupil']['thick'],
                                                         self.design['LS']['pad'], self.design['LS']['N'])) )
            else:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_quart_" + \
                                                         "{0:s}{1:02d}D{2:02d}_clear_N{3:04d}.dat".format(self.design['LS']['shape'],
                                                         self.design['LS']['id'], self.design['LS']['od'], self.design['LS']['N'])) )

        if self.design['LS']['aligntol'] is not None and ('LDZ fname' not in self.fileorg or self.fileorg['LDZ fname'] is None):
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_quart_" + \
                                                          "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}_Tol{8:02d}_N{9:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['Pupil']['prim'], self.design['Pupil']['secobs'],
                                                          self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'],
                                                          self.design['LS']['aligntol'], self.design['LS']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_quart_" + \
                                                          "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}Pad{5:02d}_Tol{6:02d}_N{7:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['Pupil']['secobs'], self.design['Pupil']['thick'],
                                                          self.design['LS']['pad'], self.design['LS']['aligntol'], self.design['LS']['N'])) )
            else:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_quart_" + \
                                                          "{0:s}{1:02d}D{2:02d}_clear_Tol{3:02d}_N{4:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['LS']['aligntol'], self.design['LS']['N'])) )
        self.check_ampl_input_files()
    def write_ampl(self, overwrite=False, override_infile_status=False, ampl_src_fname=None, verbose=True):
        if self.ampl_infile_status is False and not override_infile_status:
            if verbose:
                logging.warning("Error: the most recent input file check for this design configuration failed.")
                logging.warning("The override_infile_status switch is off, so write_ampl() will now abort.")
                logging.warning("See previous warnings in the log to see what file was missing during the initialization")
            return 2
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
                if verbose:
                    logging.warning("Error: {0} already exists and overwrite switch is off, so write_ampl() will now abort".format(self.fileorg['ampl src fname']))
                return 1
        elif not os.path.exists(self.fileorg['ampl src dir']):
            os.mkdir(self.fileorg['ampl src dir'])
            if verbose:
                logging.info("Created new AMPL source code directory, {0:s}".format(self.fileorg['ampl src dir']))
        mod_fobj = open(self.fileorg['ampl src fname'], "w")
 
        header = """\
        # AMPL program to optimize a quarter-plane symmetric SPLC
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
            param rho0 := {2:0.2f};
            param rho1 := {3:0.2f};
            param rho2 := {4:0.2f};
            param ang := {5:d};
            
            #---------------------
            param N_A := {6:d};				# discretization parameter (apodizer)
            param M := {7:d};				# discretization parameter (FPM)
            param fpmres := {8:d};			# discretization parameter (FPM)
            param N_L := {9:d};				# discretization parameter (Lyot plane)
            param Nimg := {10:d};           # discretization parameter (image)
                                  
            #---------------------
            param bw := {11:0.2f};
            param Nlam := {12:d};
            
            #---------------------
            """.format(self.design['Image']['c'], self.design['LS']['s'], self.design['FPM']['R0']+self.design['Image']['dR'],
                       self.design['FPM']['R0'], self.design['FPM']['R1'], self.design['FPM']['openang'],
                       self.design['Pupil']['N'], self.design['FPM']['M'], self.design['FPM']['fpmres'], self.design['LS']['N'],
                       self.design['Image']['Nimg'], self.design['Image']['bw+'], self.design['Image']['Nlam'])
        else:
            params = """
            #---------------------
 
            param pi:= 4*atan(1);
 
            #---------------------
            param c := {0:.2f};
 
            #---------------------
            param rho0 := {1:0.2f};
            param rho1 := {2:0.2f};
            param rho2 := {3:0.2f};
            param ang := {4:d};
            
            #---------------------
            param N_A := {5:d};				# discretization parameter (apodizer)
            param M := {6:d};				# discretization parameter (FPM)
            param fpmres := {7:d};			# discretization parameter (FPM)
            param N_L := {8:d};				# discretization parameter (Lyot plane)
            param Nimg := {9:d};           # discretization parameter (image)
                                  
            #---------------------
            param bw := {10:0.2f};
            param Nlam := {11:d};
            
            #---------------------
            """.format(self.design['Image']['c'], self.design['FPM']['R0']+self.design['Image']['dR'], self.design['FPM']['R0'],
                       self.design['FPM']['R1'], self.design['FPM']['openang'],
                       self.design['Pupil']['N'], self.design['FPM']['M'], self.design['FPM']['fpmres'],
                       self.design['LS']['N'], self.design['Image']['Nimg'], self.design['Image']['bw+'], self.design['Image']['Nlam'])

        define_coords = """
        #---------------------
        # steps in each plane
        param dx := 1/(2*N_A);
        param dy := dx;
        
        param dmx := 1/fpmres;
        param dmy := dmx;

        param du := 1/(2*N_L);
        param dv := du;
        
        param dxi := rho2/Nimg;
        param deta := dxi;
 
        #---------------------
        # coordinate vectors in each plane
        set Xs := setof {i in 0.5..N_A-0.5 by 1} i*dx;
        set Ys := setof {j in 0.5..N_A-0.5 by 1} j*dy;
        
        set MXs := setof {i in 0.5..M-0.5 by 1} i*dmx;
        set MYs := setof {j in 0.5..M-0.5 by 1} j*dmy;

        set Us := setof {i in 0.5..N_L-0.5 by 1} i*du;
        set Vs := setof {j in 0.5..N_L-0.5 by 1} j*dv;

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
            param LS {{u in Us, v in Vs}};
            read {{v in Vs, u in Us}} LS[u,v] < "{4:s}";
            close "{5:s}";
            
            # Load Lyot dark zone
            param LDZ {{u in Us, v in Vs}};
            read {{v in Vs, u in Us}} LDZ[u,v] < "{6:s}";
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
            param LS {{u in Us, v in Vs}};
            
            read {{v in Vs, u in Us}} LS[u,v] < "{4:s}";
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

        if self.design['Pupil']['edge'] is 'floor': # floor to binary
            define_pupil_and_telap = """
            #---------------------

            set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] == 1} (x,y);
            param TelApProp {x in Xs, y in Ys};
            let {x in Xs, y in Ys} TelApProp[x,y] := 0;
            let {(x,y) in Pupil} TelApProp[x,y] := 1;
            """
        elif self.design['Pupil']['edge'] is 'round': # round to binary
            define_pupil_and_telap = """
            #---------------------

            set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] > 0.5} (x,y);
            param TelApProp {x in Xs, y in Ys};
            let {x in Xs, y in Ys} TelApProp[x,y] := 0;
            let {(x,y) in Pupil} TelApProp[x,y] := 1;
            """
        else: # gray, default
            define_pupil_and_telap = """
            #---------------------

            set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] > 0} (x,y);
            param TelApProp {x in Xs, y in Ys};
            let {x in Xs, y in Ys} TelApProp[x,y] := 0;
            let {(x,y) in Pupil} TelApProp[x,y] := TelAp[x,y];
            """

        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon'] is not None: 
            sets_and_arrays_part1 = """
            #---------------------

            set FPMtrans := setof {mx in MXs, my in MYs: FPM[mx,my] > 0} (mx,my);
            set FPMall := setof {mx in MXs, my in MYs: FPM[mx,my] >= 0} (mx,my);
            set Lyot := setof {u in Us, v in Vs: LS[u,v] > 0} (u,v);
            set LyotDarkZone := setof {u in Us, v in Vs: LDZ[u,v] > 0} (u,v);

            param TR := sum {(x,y) in Pupil} TelApProp[x,y]*dx*dy; # Transmission of the Pupil. Used for calibration.
            
            var A {x in Xs, y in Ys} >= 0, <= 1, := 0.5;
            """
        else:
            sets_and_arrays_part1 = """
            #---------------------

            set FPMtrans := setof {mx in MXs, my in MYs: FPM[mx,my] > 0} (mx,my);
            set FPMall := setof {mx in MXs, my in MYs: FPM[mx,my] >= 0} (mx,my);
            set Lyot := setof {u in Us, v in Vs: LS[u,v] > 0} (u,v);

            param TR := sum {(x,y) in Pupil} TelApProp[x,y]*dx*dy; # Transmission of the Pupil. Used for calibration.
            
            var A {x in Xs, y in Ys} >= 0, <= 1, := 0.5;
            """

        if self.design['FPM']['orient'] is 'H': # Horizontal bowtie FoV
            sets_and_arrays_part2 = """
            #---------------------
            set DarkHole := setof {xi in Xis, eta in Etas:
                sqrt(xi^2+eta^2) >= rho0 && 
                sqrt(xi^2+eta^2) <= rho2 &&
                eta <= xi*tan(ang/2*pi/180)} (xi,eta);
            """
        else: # Vertical bowtie FoV
            sets_and_arrays_part2 = """
            #---------------------
            set DarkHole := setof {xi in Xis, eta in Etas:
                sqrt(xi^2+eta^2) >= rho0 && 
                sqrt(xi^2+eta^2) <= rho2 &&
                eta >= xi*tan(ang/2*pi/180)} (xi,eta);
            """
 
        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon'] is not None: 
            field_propagation = """
            #---------------------
            var EB_real_X {mx in MXs, y in Ys, lam in Ls};
            var EB_real {mx in MXs, my in MYs, lam in Ls};
            
            subject to st_EB_real_X {mx in MXs, y in Ys, lam in Ls}:
                EB_real_X[mx,y,lam] = 2*sum {x in Xs: (x,y) in Pupil} TelApProp[x,y]*A[x,y]*cos(2*pi*x*mx/lam)*dx;
            subject to st_EB_real {(mx, my) in FPMtrans, lam in Ls}:
                EB_real[mx,my,lam] = 2/lam*sum {y in Ys} EB_real_X[mx,y,lam]*cos(2*pi*y*my/lam)*dy;
            
            #---------------------
            var EC_real_X {u in Us, my in MYs, lam in Ls};
            var EC_real {u in Us, v in Vs, lam in Ls};
            
            subject to st_EC_real_X {u in Us, my in MYs, lam in Ls}:
                EC_real_X[u,my,lam] = 2*sum {mx in MXs: (mx,my) in FPMtrans} FPM[mx,my]*EB_real[mx,my,lam]*cos(2*pi*u*mx/lam)*dmx;
            subject to st_EC_real {(u,v) in Lyot union LyotDarkZone, lam in Ls}:
                EC_real[u,v,lam] = 2/lam*sum {my in MYs} EC_real_X[u,my,lam]*cos(2*pi*v*my/lam)*dmy;
            
            #---------------------
            var ED_real_X {xi in Xis, v in Vs, lam in Ls};
            var ED_real {xi in Xis, eta in Etas, lam in Ls};
            
            subject to st_ED_real_X {xi in Xis, v in Vs, lam in Ls}: 
                ED_real_X[xi,v,lam] = 2*sum {u in Us: (u,v) in Lyot} LS[u,v]*EC_real[u,v,lam]*cos(2*pi*u*xi/lam)*du;
            subject to st_ED_real {(xi, eta) in DarkHole, lam in Ls}: 
                ED_real[xi,eta,lam] = 2/lam*sum {v in Vs} ED_real_X[xi,v,lam]*cos(2*pi*v*eta/lam)*dv;
            
            #---------------------
            var EB00_real_X {mx in MXs, y in Ys};
            var EB00_real {mx in MXs, my in MYs};
            var EC00_real_X {u in Us, my in MYs};
            var EC00_real {u in Us, v in Vs};
            var ED00_real := 0.0;

            subject to st_EB00_real_X {mx in MXs, y in Ys}:
                EB00_real_X[mx,y] = 2*sum {x in Xs: (x,y) in Pupil} TelApProp[x,y]*A[x,y]*cos(2*pi*x*mx)*dx;
            subject to st_EB00_real {(mx, my) in FPMall}: 
                EB00_real[mx,my] = 2*sum {y in Ys} EB00_real_X[mx,y]*cos(2*pi*y*my)*dy;
            subject to st_EC00_real_X {u in Us, my in MYs}:
                EC00_real_X[u,my] = 2*sum {mx in MXs: (mx,my) in FPMall} EB00_real[mx,my]*cos(2*pi*u*mx)*dmx;
            subject to st_EC00_real {(u,v) in Lyot}:
                EC00_real[u,v] = 2*sum {my in MYs} EC00_real_X[u,my]*cos(2*pi*v*my)*dmy;
            subject to st_ED00_real:
                ED00_real = 4.*sum {u in Us, v in Vs: (u,v) in Lyot} LS[u,v]*EC00_real[u,v]*du*dv;
            """
        else:
            field_propagation = """
            #---------------------
            var EB_real_X {mx in MXs, y in Ys, lam in Ls};
            var EB_real {mx in MXs, my in MYs, lam in Ls};
            
            subject to st_EB_real_X {mx in MXs, y in Ys, lam in Ls}:
                EB_real_X[mx,y,lam] = 2*sum {x in Xs: (x,y) in Pupil} TelApProp[x,y]*A[x,y]*cos(2*pi*x*mx/lam)*dx;
            subject to st_EB_real {(mx, my) in FPMtrans, lam in Ls}:
                EB_real[mx,my,lam] = 2/lam*sum {y in Ys} EB_real_X[mx,y,lam]*cos(2*pi*y*my/lam)*dy;
            
            #---------------------
            var EC_real_X {u in Us, my in MYs, lam in Ls};
            var EC_real {u in Us, v in Vs, lam in Ls};
            
            subject to st_EC_real_X {u in Us, my in MYs, lam in Ls}:
                EC_real_X[u,my,lam] = 2*sum {mx in MXs: (mx,my) in FPMtrans} FPM[mx,my]*EB_real[mx,my,lam]*cos(2*pi*u*mx/lam)*dmx;
            subject to st_EC_real {(u,v) in Lyot, lam in Ls}:
                EC_real[u,v,lam] = 2/lam*sum {my in MYs} EC_real_X[u,my,lam]*cos(2*pi*v*my/lam)*dmy;
            
            #---------------------
            var ED_real_X {xi in Xis, v in Vs, lam in Ls};
            var ED_real {xi in Xis, eta in Etas, lam in Ls};
            
            subject to st_ED_real_X {xi in Xis, v in Vs, lam in Ls}: 
                ED_real_X[xi,v,lam] = 2*sum {u in Us: (u,v) in Lyot} LS[u,v]*EC_real[u,v,lam]*cos(2*pi*u*xi/lam)*du;
            subject to st_ED_real {(xi, eta) in DarkHole, lam in Ls}: 
                ED_real[xi,eta,lam] = 2/lam*sum {v in Vs} ED_real_X[xi,v,lam]*cos(2*pi*v*eta/lam)*dv;
            
            #---------------------
            var EB00_real_X {mx in MXs, y in Ys};
            var EB00_real {mx in MXs, my in MYs};
            var EC00_real_X {u in Us, my in MYs};
            var EC00_real {u in Us, v in Vs};
            var ED00_real := 0.0;

            subject to st_EB00_real_X {mx in MXs, y in Ys}:
                EB00_real_X[mx,y] = 2*sum {x in Xs: (x,y) in Pupil} TelApProp[x,y]*A[x,y]*cos(2*pi*x*mx)*dx;
            subject to st_EB00_real {(mx, my) in FPMall}: 
                EB00_real[mx,my] = 2*sum {y in Ys} EB00_real_X[mx,y]*cos(2*pi*y*my)*dy;
            subject to st_EC00_real_X {u in Us, my in MYs}:
                EC00_real_X[u,my] = 2*sum {mx in MXs: (mx,my) in FPMall} EB00_real[mx,my]*cos(2*pi*u*mx)*dmx;
            subject to st_EC00_real {(u,v) in Lyot}:
                EC00_real[u,v] = 2*sum {my in MYs} EC00_real_X[u,my]*cos(2*pi*v*my)*dmy;
            subject to st_ED00_real:
                ED00_real = 4.*sum {u in Us, v in Vs: (u,v) in Lyot} LS[u,v]*EC00_real[u,v]*du*dv;
            """

        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon'] is not None:
            constraints = """
            #---------------------
            maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;
           
            subject to Lyot_aligntol_constr_pos {(x,y) in LyotDarkZone, lam in Ls}:
                EC_real[x,y,lam] <= 10^-s;
            subject to Lyot_aligntol_constr_neg {(x,y) in LyotDarkZone, lam in Ls}:
                EC_real[x,y,lam] >= -10^-s;
            subject to sidelobe_zero_real_pos {(xi,eta) in DarkHole, lam in Ls}:
                ED_real[xi,eta,lam] <= 10^(-c/2)*ED00_real/lam/sqrt(2.);
            subject to sidelobe_zero_real_neg {(xi,eta) in DarkHole, lam in Ls}:
                ED_real[xi,eta,lam] >= -10^(-c/2)*ED00_real/lam/sqrt(2.);
            """
        else:
            constraints = """
            #---------------------
            maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;
            
            subject to sidelobe_zero_real_pos {(xi,eta) in DarkHole, lam in Ls}:
                ED_real[xi,eta,lam] <= 10^(-c/2)*ED00_real/lam/sqrt(2.);
            subject to sidelobe_zero_real_neg {(xi,eta) in DarkHole, lam in Ls}:
                ED_real[xi,eta,lam] >= -10^(-c/2)*ED00_real/lam/sqrt(2.);
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
            gurobi_opt_str += " lpmethod=2"
            if self.solver['convtol'] is not None:
                gurobi_opt_str += " barconvtol={0:.1e}".format(np.power(10,-self.solver['convtol']))
            if self.solver['method'] is 'barhom':
                gurobi_opt_str += " barhomogeneous=1"
            if self.solver['crossover'] is True:
                gurobi_opt_str += " crossoverbasis=1"
            else:
                gurobi_opt_str += " crossover=0"
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
        mod_fobj.write( textwrap.dedent(define_pupil_and_telap) )
        mod_fobj.write( textwrap.dedent(sets_and_arrays_part1) )
        mod_fobj.write( textwrap.dedent(sets_and_arrays_part2) )
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
        return 0

    def write_slurm_script(self, queue_spec='auto', account='s1649', email=None, arch=None, overwrite=False, verbose=True):
        if os.path.exists(self.fileorg['slurm fname']):
            if overwrite == True:
                if verbose:
                    logging.warning("Warning: Overwriting the existing copy of {0}".format(self.fileorg['slurm fname']))
            else:
                if verbose:
                    logging.warning("Error: {0} already exists and overwrite switch is off, so write_slurm_script() will now abort".format(self.fileorg['slurm fname']))
                return 1
        elif not os.path.exists(self.fileorg['slurm dir']):
            os.mkdir(self.fileorg['slurm dir'])
            if verbose:
                logging.info("Created new slurm script directory, {0:s}".format(self.fileorg['slurm dir']))

        bash_fobj = open(self.fileorg['slurm fname'], "w") 

        if email is not None: 
            header = """\
            #!/bin/bash

            #PBS -V
            #PBS -m e -M {0:s}
            """.format(email)
        else:
            header = """\
            #! /bin/bash

            """

        set_job = """\
        #SBATCH --job-name={0:s}
        #SBATCH -o {1:s}
        #SBATCH --account={2:s}
        """.format(self.fileorg['job name'], self.fileorg['log fname'], account)
       
        if arch is not None: # can be 'hasw' for Haswell only 
            set_node = """\
            #SBATCH --constraint={0:s}
            #SBATCH --ntasks=1 --nodes=1
            """.format(arch)
        else:
            set_node = """\
            #SBATCH --ntasks=1 --nodes=1
            """.format(arch)

        if queue_spec is 'auto':
            if self.design['LS']['aligntol'] is None:
                time_est_hrs = int(np.ceil(1.5*(self.design['Pupil']['N']/125.)**2*(self.design['Image']['Nlam']/3.)**3))
            else:
                time_est_hrs = int(np.ceil(3*(self.design['Pupil']['N']/125.)**2*(self.design['Image']['Nlam']/3.)**3))
            if time_est_hrs > 12:
                set_queue = """
                #SBATCH --qos=long
                #SBATCH --time={0:02d}:00:00
                """.format(np.min([24, time_est_hrs]))
            else:
                set_queue = """
                #SBATCH --qos=allnccs
                #SBATCH --time={0:02d}:00:00
                """.format(time_est_hrs)
        elif queue_spec is '1h':
            set_queue = """
            #SBATCH --qos=debug
            #SBATCH --time=1:00:00
            """
        elif queue_spec is '12h':
            set_queue = """
            #SBATCH --qos=allnccs
            #SBATCH --time=12:00:00
            """
        else:
            set_queue = """
            #SBATCH --qos=long
            #SBATCH --time=24:00:00
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
        
        exit 0
        """.format(self.fileorg['ampl src fname'])

        bash_fobj.write( textwrap.dedent(header) )
        bash_fobj.write( textwrap.dedent(set_job) )
        bash_fobj.write( textwrap.dedent(set_node) )
        bash_fobj.write( textwrap.dedent(set_queue) )
        bash_fobj.write( textwrap.dedent(intel_module) )
        bash_fobj.write( textwrap.dedent(monitor_mem) )
        bash_fobj.write( textwrap.dedent(call_ampl) )

        bash_fobj.close()
        if verbose:
            logging.info("Wrote %s"%self.fileorg['slurm fname'])
        return 0

class NdiayeAPLC(LyotCoronagraph): # Image-constrained APLC following N'Diaye et al. (2015, 2016)
    _design_fields = OrderedDict([ ( 'Pupil', OrderedDict([('N',(int, 250)), ('prim',(str, 'hex3')), ('secobs',(str, None)), 
                                                           ('thick',(str, '025')), ('centobs',(bool, True)),
                                                           ('gap',(int, 1)), ('edge',(str, 'gray'))]) ),
                                   ( 'FPM', OrderedDict([('rad',(float, 4.)), ('M',(int, 60))]) ),
                                   ( 'LS', OrderedDict([('shape',(str, 'ann')), ('id',(int, 20)), ('od',(int, None)), ('obscure',(int, 0)),
                                                        ('pad',(int, 0)), ('aligntol',(int, None)), ('aligntolcon',(float, 3.))]) ),
                                   ( 'Image', OrderedDict([('c',(float, 10.)), ('ida',(float, -0.5)), ('oda',(float, 10.)),
                                                           ('bw',(float, 0.10)), ('Nlam',(int, 1)), ('fpres',(int,2)),
                                                           ('wingang',(float, None)), ('incon',(float, None)), ('wingcon',(float, None))]) ) ])
    _eval_fields =   { 'Pupil': _design_fields['Pupil'], 'FPM': _design_fields['FPM'], \
                       'LS': _design_fields['LS'], 'Image': _design_fields['Image'], \
                       'Tel': {'TelAp diam':(float, 12.)}, 'Target': {}, 'Aber': {}, 'WFSC': {} }
    _LS_OD_map = {'hex1':76, 'hex2':82, 'hex3':81, 'hex4':82, 'pie08':90, 'pie12':90, 'key24':90}
    _prim_secobs_map = {'hex1':'X', 'hex2':'X', 'hex3':'X', 'hex4':'X', 'pie08':'Cross', 'pie12':'Cross', 'key24':'Cross'}

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
        # If secondary obscuration is not specified, look up default based on primary key
        if self.design['Pupil']['secobs'] is None:
            self.design['Pupil']['secobs'] = self._prim_secobs_map[self.design['Pupil']['prim']]
        # Unless Nlam is explicitly specified, set the number of wavelength samples according to the bandwidth
        if self.design['Image']['Nlam'] == 1 and self.design['Image']['bw'] > 0:
            self.design['Image']['Nlam'] = int(np.ceil(self.design['Image']['bw']/(0.10/3)))
        # Finally, set a private attribute for the number of image plane samples between the center and the outer constraint angle
        if self.design['Image']['wingang'] is not None:
            self.design['Image']['Nimg'] = int( np.ceil( self.design['Image']['fpres']*self.design['Image']['wingang']/(1. - self.design['Image']['bw']/2) ) )
        else:
            self.design['Image']['Nimg'] = int( np.ceil( self.design['Image']['fpres']*self.design['Image']['oda']/(1. - self.design['Image']['bw']/2) ) )
        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon']:
            # The Lyot dark zone field suppression factor decreases with the square of pupil array size. The units of input parameter are arbitrarily normalized to N=125.
            self.design['LS']['s'] = self.design['LS']['aligntolcon'] - 2*np.log10(self.design['Pupil']['N']/125.)
        if self.design['LS']['od'] is None:
            if self.design['Pupil']['prim'] in self._LS_OD_map:
                self.design['LS']['od'] = self._LS_OD_map[self.design['Pupil']['prim']]
            else:
                self.design['LS']['od'] = 0.80
        if verbose: # Print summary of the set parameters
            logging.info("Design parameters: {}".format(self.design))
            logging.info("Optimization and solver parameters: {}".format(self.solver))
            logging.info("File organization parameters: {}".format(self.fileorg))
     
        self.amplname_coron = "APLC_full"
        self.telap_descrip = "{0:s}{1:s}{2:s}cobs{3:d}gap{4:d}_N{5:04d}".format(self.design['Pupil']['prim'], self.design['Pupil']['secobs'], self.design['Pupil']['thick'], \
                                                                                int(self.design['Pupil']['centobs']), self.design['Pupil']['gap'], self.design['Pupil']['N'])
        self.amplname_pupil = "{0:s}{1:s}".format(self.telap_descrip, self.design['Pupil']['edge'][0])

        self.amplname_fpm = "FPM{:02}M{:03}".format(int(round(100*self.design['FPM']['rad'])), self.design['FPM']['M'])
        if self.design['LS']['obscure'] == 2: # LS includes primary and secondary aperture features
            self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}".format(self.design['LS']['shape'], self.design['LS']['id'], \
                               self.design['LS']['od'], self.design['Pupil']['prim'], self.design['Pupil']['secobs'], self.design['Pupil']['thick'], \
                               int(self.design['Pupil']['centobs']), self.design['LS']['pad'])
        elif self.design['LS']['obscure'] == 1: # LS includes secondary aperture features
            self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}{3:s}{4:s}Pad{5:02d}".format(self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'], \
                               self.design['Pupil']['secobs'], self.design['Pupil']['thick'], self.design['LS']['pad'])
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
        if self.solver['crossover'] is not None:
            self.amplname_solver += "cross"
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
                self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_full_" + self.telap_descrip + ".dat") )

            if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
                self.fileorg['FPM fname'] = os.path.join( self.fileorg['FPM dir'], "FPM_full_occspot_M{:03}.dat".format(self.design['FPM']['M']) )

            if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
                if self.design['LS']['obscure'] == 2:
                    self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}_N{8:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['Pupil']['prim'], self.design['Pupil']['secobs'],
                                                             self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'],
                                                             self.design['Pupil']['N'])) )
                elif self.design['LS']['obscure'] == 1:
                    self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}Pad{5:02d}_N{6:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['Pupil']['secobs'], self.design['Pupil']['thick'],
                                                             self.design['LS']['pad'], self.design['Pupil']['N'])) )
                else:
                    self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_full_" + \
                                                             "{0:s}{1:02d}D{2:02d}_clear_N{3:04d}.dat".format(self.design['LS']['shape'],
                                                             self.design['LS']['id'], self.design['LS']['od'], self.design['Pupil']['N'])) )

            if self.design['LS']['aligntol'] is not None and ('LDZ fname' not in self.fileorg or self.fileorg['LDZ fname'] is None):
                if self.design['LS']['obscure'] == 2:
                    self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_full_" + \
                                                              "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}_Tol{8:02d}_N{9:04d}.dat".format(
                                                              self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                              self.design['Pupil']['prim'], self.design['Pupil']['secobs'],
                                                              self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'],
                                                              self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
                elif self.design['LS']['obscure'] == 1:
                    self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_full_" + \
                                                              "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}Pad{5:02d}_Tol{6:02d}_N{7:04d}.dat".format(
                                                              self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                              self.design['Pupil']['secobs'], self.design['Pupil']['thick'],
                                                              self.design['LS']['pad'], self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
                else:
                    self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_full_" + \
                                                              "{0:s}{1:02d}D{2:02d}_clear_Tol{3:02d}_N{4:04d}.dat".format(
                                                              self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                              self.design['LS']['aligntol'], self.design['Pupil']['N'])) )

            self.check_ampl_input_files()
                                                              
    def write_ampl(self, overwrite=False):
        logging.info("Writing the AMPL program")

    def get_onax_psf(self, fp2res=8, rho_inc=0.25, rho_out=None, Nlam=None): # for APLC class
        if self.design['Pupil']['edge'] is 'floor': # floor to binary
            TelAp_p = np.floor(np.loadtxt(self.fileorg['TelAp fname'])).astype(int)
        elif self.design['Pupil']['edge'] is 'round': # round to binary
            TelAp_p = np.round(np.loadtxt(self.fileorg['TelAp fname'])).astype(int)
        else: # keey it gray
            TelAp_p = np.loadtxt(self.fileorg['TelAp fname'])
        A_col = np.loadtxt(self.fileorg['sol fname'])[:,-1]
        FPM_p = np.loadtxt(self.fileorg['FPM fname'])
        LS_p = np.loadtxt(self.fileorg['LS fname'])
        A_p = A_col.reshape(TelAp_p.shape)
        if isinstance(self, QuarterplaneAPLC):
            TelAp = np.concatenate((np.concatenate((TelAp_p[::-1,::-1], TelAp_p[:,::-1]),axis=0),
                                    np.concatenate((TelAp_p[::-1,:], TelAp_p),axis=0)), axis=1)
            A = np.concatenate((np.concatenate((A_p[::-1,::-1], A_p[:,::-1]),axis=0),
                                np.concatenate((A_p[::-1,:], A_p),axis=0)), axis=1)
            FPM = np.concatenate((np.concatenate((FPM_p[::-1,::-1], FPM_p[:,::-1]),axis=0),
                                  np.concatenate((FPM_p[::-1,:], FPM_p),axis=0)), axis=1)
            LS = np.concatenate((np.concatenate((LS_p[::-1,::-1], LS_p[:,::-1]),axis=0),
                                 np.concatenate((LS_p[::-1,:], LS_p),axis=0)), axis=1)
        elif isinstance(self, HalfplaneAPLC):
            TelAp = np.concatenate((TelAp_p[:,::-1], TelAp_p), axis=1)
            A = np.concatenate((A_p[:,::-1], A_p), axis=1)
            FPM = np.concatenate((np.concatenate((FPM_p[::-1,::-1], FPM_p[:,::-1]),axis=0),
                                  np.concatenate((FPM_p[::-1,:], FPM_p),axis=0)), axis=1)
            LS = np.concatenate((LS_p[:,::-1], LS_p), axis=1)
        else:
            TelAp = TelAp_p
            A = A_p
            FPM = FPM_p
            LS = LS_p

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

    def get_metrics(self, fp2res=16, rho_out=None, Nlam=None, verbose=True): # for APLC class
        TelAp_basename = os.path.basename(self.fileorg['TelAp fname'])
        gapstr_beg = TelAp_basename.find('gap')
        TelAp_nopad_basename = TelAp_basename.replace(TelAp_basename[gapstr_beg:gapstr_beg+4], 'gap0')
        TelAp_nopad_fname = os.path.join( os.path.dirname(self.fileorg['TelAp fname']), TelAp_nopad_basename )
        #if self.design['Pupil']['edge'] is 'floor': # floor to binary
        #    TelAp_p = np.floor(np.loadtxt(self.fileorg['TelAp fname'])).astype(int)
        #elif self.design['Pupil']['edge'] is 'round': # round to binary
        #    TelAp_p = np.round(np.loadtxt(self.fileorg['TelAp fname'])).astype(int)
        #else:
        #    TelAp_p = np.loadtxt(self.fileorg['TelAp fname'])
        TelAp_p = np.loadtxt(TelAp_nopad_fname)
        A_col = np.loadtxt(self.fileorg['sol fname'])[:,-1]
        LS_p = np.loadtxt(self.fileorg['LS fname'])
        A_p = A_col.reshape(TelAp_p.shape)
        if isinstance(self, QuarterplaneAPLC):
            TelAp = np.concatenate((np.concatenate((TelAp_p[::-1,::-1], TelAp_p[:,::-1]),axis=0),
                                    np.concatenate((TelAp_p[::-1,:], TelAp_p),axis=0)), axis=1)
            A = np.concatenate((np.concatenate((A_p[::-1,::-1], A_p[:,::-1]),axis=0),
                                np.concatenate((A_p[::-1,:], A_p),axis=0)), axis=1)
            LS = np.concatenate((np.concatenate((LS_p[::-1,::-1], LS_p[:,::-1]),axis=0),
                                 np.concatenate((LS_p[::-1,:], LS_p),axis=0)), axis=1)
        elif isinstance(self, HalfplaneAPLC):
            TelAp = np.concatenate((TelAp_p[:,::-1], TelAp_p), axis=1)
            A = np.concatenate((A_p[:,::-1], A_p), axis=1)
            LS = np.concatenate((LS_p[:,::-1], LS_p), axis=1)
        else:
            TelAp = TelAp_p
            A = A_p
            LS = LS_p

        self.eval_metrics['apod nb res ratio'] = np.sum(np.abs(A - np.round(A)))/np.sum(TelAp)
        D = 1.
        N = self.design['Pupil']['N']
        if Nlam is None:
            Nlam = self.design['Image']['Nlam']
        if rho_out is None:
            rho_out = self.design['Image']['oda']
        dx = (D/2)/N
        dy = dx
        xs = np.matrix(np.linspace(-N+0.5, N-0.5, 2*N)*dx)
        ys = xs.copy()
        M_fp2 = int(np.ceil(rho_out*fp2res))
        dxi = 1./fp2res
        xis = np.matrix(np.linspace(-M_fp2+0.5, M_fp2-0.5, 2*M_fp2)*dxi)
        etas = xis.copy()
        wrs = np.linspace(1.-self.design['Image']['bw']/2, 1.+self.design['Image']['bw']/2, Nlam)
        XXs = np.asarray(np.dot(np.matrix(np.ones(xis.shape)).T, xis))
        YYs = np.asarray(np.dot(etas.T, np.matrix(np.ones(etas.shape))))
        RRs = np.sqrt(XXs**2 + YYs**2)
        p7ap_ind = np.less_equal(RRs, 0.7)

        intens_D_0_polychrom = np.zeros((Nlam, 2*M_fp2, 2*M_fp2))
        intens_D_0_peak_polychrom = np.zeros((Nlam, 1))
        intens_TelAp_polychrom = np.zeros((Nlam, 2*M_fp2, 2*M_fp2))
        intens_TelAp_peak_polychrom = np.zeros((Nlam, 1))
        for wi, wr in enumerate(wrs):
            Psi_D_0 = dx*dy/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(xis.T, xs)), TelAp*A*LS),
                                             np.exp(-1j*2*np.pi/wr*np.dot(xs.T, xis)))
            intens_D_0_polychrom[wi] = np.power(np.absolute(Psi_D_0), 2)
            intens_D_0_peak_polychrom[wi] = (np.sum(TelAp*A*LS)*dx*dy/wr)**2
            Psi_TelAp = dx*dy/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(xis.T, xs)), TelAp),
                                               np.exp(-1j*2*np.pi/wr*np.dot(xs.T, xis)))
            intens_TelAp_polychrom[wi] = np.power(np.absolute(Psi_TelAp), 2)
            intens_TelAp_peak_polychrom[wi] = (np.sum(TelAp)*dx*dy/wr)**2

        intens_D_0 = np.mean(intens_D_0_polychrom, axis=0)
        intens_D_0_peak = np.mean(intens_D_0_peak_polychrom)
        intens_TelAp = np.mean(intens_TelAp_polychrom, axis=0)
        intens_TelAp_peak = np.mean(intens_TelAp_peak_polychrom)

        fwhm_ind_APLC = np.greater_equal(intens_D_0, intens_D_0_peak/2)
        fwhm_ind_TelAp = np.greater_equal(intens_TelAp, intens_TelAp_peak/2)

        fwhm_sum_TelAp = np.sum(intens_TelAp[fwhm_ind_TelAp])*dxi*dxi
        fwhm_sum_APLC = np.sum(intens_D_0[fwhm_ind_APLC])*dxi*dxi
        p7ap_sum_TelAp = np.sum(intens_TelAp[p7ap_ind])*dxi*dxi
        p7ap_sum_APLC = np.sum(intens_D_0[p7ap_ind])*dxi*dxi

        self.eval_metrics['inc energy'] = np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['tot thrupt'] = np.sum(intens_D_0*dxi*dxi)/np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['fwhm thrupt'] = fwhm_sum_APLC/np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['p7ap thrupt'] = p7ap_sum_APLC/np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['rel fwhm thrupt'] = fwhm_sum_APLC/fwhm_sum_TelAp
        self.eval_metrics['rel p7ap thrupt'] = p7ap_sum_APLC/p7ap_sum_TelAp
        self.eval_metrics['fwhm area'] = np.sum(fwhm_ind_APLC)*dxi*dxi
        if verbose:
            print("////////////////////////////////////////////////////////")
            print("{:s}".format(self.fileorg['job name']))
            print("Incident energy on aperture (dimensionless): {:.3f}".format(self.eval_metrics['inc energy']))
            print("Non-binary residuals, as a percentage of clear telescope aperture area: {:.2f}%".format(100*self.eval_metrics['apod nb res ratio']))
            print("Band-averaged total throughput: {:.2f}%".format(100*self.eval_metrics['tot thrupt']))
            print("Band-averaged half-max throughput: {:.2f}%".format(100*self.eval_metrics['fwhm thrupt']))
            print("Band-averaged r=.7 lam/D throughput: {:.2f}%".format(100*self.eval_metrics['p7ap thrupt']))
            print("Band-averaged relative half-max throughput: {:.2f}%".format(100*self.eval_metrics['rel fwhm thrupt']))
            print("Band-averaged relative r=0.7 lam/D throughput: {:.2f}%".format(100*self.eval_metrics['rel p7ap thrupt']))
            print("Band-averaged FWHM PSF area / (lambda0/D)^2: {:.2f}".format(self.eval_metrics['fwhm area']))
    
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

        if 'job name' not in self.fileorg or self.fileorg['job name'] is None:
            self.fileorg['job name'] = os.path.basename(self.fileorg['ampl src fname'])[:-4]
 
        if 'sol fname' not in self.fileorg or self.fileorg['sol fname'] is None:
            sol_fname_tail = "ApodSol_" + self.fileorg['job name'] + ".dat"
            self.fileorg['sol fname'] = os.path.join(self.fileorg['sol dir'], sol_fname_tail)

        if 'slurm fname' not in self.fileorg or self.fileorg['slurm fname'] is None:
            exec_script_fname_tail = self.fileorg['job name'] + ".sh"
            self.fileorg['slurm fname'] = os.path.join(self.fileorg['slurm dir'], exec_script_fname_tail)

        if 'log fname' not in self.fileorg or self.fileorg['log fname'] is None:
            log_fname_tail = self.fileorg['job name'] + ".log"
            self.fileorg['log fname'] = os.path.join(self.fileorg['log dir'], log_fname_tail)
 
        if 'TelAp fname' not in self.fileorg or self.fileorg['TelAp fname'] is None:
            self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_half_" + self.telap_descrip + ".dat") )

        if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
            self.fileorg['FPM fname'] = os.path.join( self.fileorg['FPM dir'], "FPM_quart_occspot_M{:03}.dat".format(self.design['FPM']['M']) )

        if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}_N{8:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['prim'], self.design['Pupil']['secobs'],
                                                         self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'],
                                                         self.design['Pupil']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}Pad{5:02d}_N{6:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['secobs'], self.design['Pupil']['thick'],
                                                         self.design['LS']['pad'], self.design['Pupil']['N'])) )
            else:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_clear_N{3:04d}.dat".format(self.design['LS']['shape'],
                                                         self.design['LS']['id'], self.design['LS']['od'], self.design['Pupil']['N'])) )

        if self.design['LS']['aligntol'] is not None and ('LDZ fname' not in self.fileorg or self.fileorg['LDZ fname'] is None):
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_half_" + \
                                                          "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}_Tol{8:02d}_N{9:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['Pupil']['prim'], self.design['Pupil']['secobs'],
                                                          self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'],
                                                          self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_half_" + \
                                                          "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}Pad{5:02d}_Tol{6:02d}_N{7:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['Pupil']['secobs'], self.design['Pupil']['thick'],
                                                          self.design['LS']['pad'], self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
            else:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_half_" + \
                                                          "{0:s}{1:02d}D{2:02d}_clear_Tol{3:02d}_N{4:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['LS']['aligntol'], self.design['Pupil']['N'])) )

        self.check_ampl_input_files()

    def write_ampl(self, overwrite=False, override_infile_status=False, ampl_src_fname=None, verbose=True):
        if self.ampl_infile_status is False and not override_infile_status:
            if verbose:
                logging.warning("Error: the most recent input file check for this design configuration failed.")
                logging.warning("The override_infile_status switch is off, so write_ampl() will now abort.")
                logging.warning("See previous warnings in the log to see what file was missing during the initialization")
            return 2
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
                if verbose:
                    logging.warning("Error: {0} already exists and overwrite switch is off, so write_ampl() will now abort".format(self.fileorg['ampl src fname']))
                return 1
        elif not os.path.exists(self.fileorg['ampl src dir']):
            os.mkdir(self.fileorg['ampl src dir'])
            if verbose:
                logging.info("Created new AMPL source code directory, {0:s}".format(self.fileorg['ampl src dir']))
        mod_fobj = open(self.fileorg['ampl src fname'], "w")

        header = """\
        # AMPL program to optimize a half-plane symmetric APLC
        # Created by {0:s} with {1:s} on {2:s} at {3:s}
        load amplgsl.dll;
        """.format(getpass.getuser(), os.path.basename(__file__), socket.gethostname(), datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))

        if True: # no Lyot plane tolerancing constraints implemented yet
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
        set Ys := setof {j in -N+0.5..N-0.5 by 1} j*dy;
        
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
            set Lyot := setof {x in Xs, y in Ys: LS[x,y] >= 0.5} (x,y);
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
            set Lyot := setof {x in Xs, y in Ys: LS[x,y] >= 0.5} (x,y);

            param TR := sum {(x,y) in Pupil} TelAp[x,y]*dx*dy; # Transmission of the Pupil. Used for calibration.
            
            var A {x in Xs, y in Ys} >= 0, <= 1, := 0.5;
            
            #---------------------
            set DarkHole := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= rho0 && sqrt(xi^2+eta^2) <= rho1} (xi,eta);
            """

        field_propagation = """
        #---------------------
        var EBm_part {mx in MXs, y in Ys, lam in Ls};
        var EBm_real {mx in MXs, my in MYs, lam in Ls};
        var EBm_imag {mx in MXs, my in MYs, lam in Ls};

        subject to st_EBm_part {mx in MXs, y in Ys, lam in Ls}:
            EBm_part[mx,y,lam] = 2*sum {x in Xs: (x,y) in Pupil} TelAp[x,y]*A[x,y]*cos(2*pi*x*mx/lam)*dx;
        subject to st_EBm_real {(mx,my) in Mask, lam in Ls}:
            EBm_real[mx,my,lam] = 1/lam*sum {y in Ys} EBm_part[mx,y,lam]*cos(2*pi*y*my/lam)*dy;
        subject to st_EBm_imag {(mx,my) in Mask, lam in Ls}:
            EBm_imag[mx,my,lam] = 1/lam*sum {y in Ys} EBm_part[mx,y,lam]*sin(2*pi*y*my/lam)*dy;
        
        #---------------------
        var EC_part_real {x in Xs, my in MYs, lam in Ls};
        var EC_part_imag {x in Xs, my in MYs, lam in Ls};
        var EC {x in Xs, y in Ys, lam in Ls};
        
        subject to st_EC_part_real {x in Xs, my in MYs, lam in Ls}:
            EC_part_real[x,my,lam] = 2*sum {mx in MXs: (mx,my) in Mask} FPM[mx,my]*EBm_real[mx,my,lam]*cos(2*pi*x*mx/lam)*dmx;
        subject to st_EC_part_imag {x in Xs, my in MYs, lam in Ls}:
            EC_part_imag[x,my,lam] = 2*sum {mx in MXs: (mx,my) in Mask} FPM[mx,my]*EBm_imag[mx,my,lam]*cos(2*pi*x*mx/lam)*dmx;
        subject to st_EC {(x,y) in Lyot, lam in Ls}:
            EC[x,y,lam] = TelAp[x,-y]*A[x,-y] - 2/lam*sum{my in MYs} ( EC_part_real[x,my,lam]*cos(2*pi*my*y/lam) + EC_part_imag[x,my,lam]*sin(2*pi*my*y/lam) )*dmy;
        
        #---------------------
        var ED_part {xi in Xis, y in Ys, lam in Ls};
        var ED_real {xi in Xis, eta in Etas, lam in Ls};
        var ED_imag {xi in Xis, eta in Etas, lam in Ls};
        
        subject to st_ED_part {xi in Xis, y in Ys, lam in Ls}:
            ED_part[xi,y,lam] = 2*sum {x in Xs: (x,y) in Lyot} LS[x,y]*EC[x,y,lam]*cos(2*pi*x*xi/lam)*dx;
        subject to st_ED_real {xi in Xis, eta in Etas, lam in Ls}:
            ED_real[xi,eta,lam] = 1/lam*sum {y in Ys} ED_part[xi,y,lam]*cos(2*pi*y*eta/lam)*dy;
        subject to st_ED_imag {xi in Xis, eta in Etas, lam in Ls}:
            ED_imag[xi,eta,lam] = 1/lam*sum {y in Ys} ED_part[xi,y,lam]*sin(2*pi*y*eta/lam)*dy;
        
        #---------------------
        var ED00_real := 0.0;
        subject to st_ED00_real: ED00_real = 2*sum {x in Xs, y in Ys: (x,y) in Lyot} (TelAp[x,y]*A[x,y])*dx*dy;
        """

        constraints = """
        #---------------------
        maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;

        subject to sidelobe_real_pos {(xi,eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] <= 10^(-c/2)*ED00_real/lam/sqrt(2.);
        subject to sidelobe_real_neg {(xi,eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] >= -10^(-c/2)*ED00_real/lam/sqrt(2.);
        subject to sidelobe_imag_pos {(xi,eta) in DarkHole, lam in Ls}: ED_imag[xi,eta,lam] <= 10^(-c/2)*ED00_real/lam/sqrt(2.);
        subject to sidelobe_imag_neg {(xi,eta) in DarkHole, lam in Ls}: ED_imag[xi,eta,lam] >= -10^(-c/2)*ED00_real/lam/sqrt(2.);
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
            gurobi_opt_str += " lpmethod=2"
            if self.solver['convtol'] is not None:
                gurobi_opt_str += " barconvtol={0:.1e}".format(np.power(10,-self.solver['convtol']))
            if self.solver['method'] is 'barhom':
                gurobi_opt_str += " barhomogeneous=1"
            if self.solver['crossover'] is True:
                gurobi_opt_str += " crossoverbasis=1"
            else:
                gurobi_opt_str += " crossover=0"
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

        param A_fin {{x in Xs, y in Ys}};
        let {{x in Xs, y in Ys}} A_fin[x,y] := 0;
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
        return 0

    def write_slurm_script(self, queue_spec='auto', account='s1649', email=None, arch=None, overwrite=False, verbose=True):
        if os.path.exists(self.fileorg['slurm fname']):
            if overwrite == True:
                if verbose:
                    logging.warning("Warning: Overwriting the existing copy of {0}".format(self.fileorg['slurm fname']))
            else:
                if verbose:
                    logging.warning("Error: {0} already exists and overwrite switch is off, so write_slurm_script() will now abort".format(self.fileorg['slurm fname']))
                return 1
        elif not os.path.exists(self.fileorg['slurm dir']):
            os.mkdir(self.fileorg['slurm dir'])
            if verbose:
                logging.info("Created new slurm script directory, {0:s}".format(self.fileorg['slurm dir']))

        bash_fobj = open(self.fileorg['slurm fname'], "w") 

        if email is not None: 
            header = """\
            #!/bin/bash

            #PBS -V
            #PBS -m e -M {0:s}
            """.format(email)
        else:
            header = """\
            #! /bin/bash

            """

        set_job = """\
        #SBATCH --job-name={0:s}
        #SBATCH -o {1:s}
        #SBATCH --account={2:s}
        """.format(self.fileorg['job name'], self.fileorg['log fname'], account)
       
        if arch is not None: # can be 'hasw' for Haswell only 
            set_node = """\
            #SBATCH --constraint={0:s}
            #SBATCH --ntasks=1 --nodes=1
            """.format(arch)
        else:
            set_node = """\
            #SBATCH --ntasks=1 --nodes=1
            """.format(arch)

        if queue_spec is 'auto':
            if self.design['LS']['aligntol'] is None:
                time_est_hrs = int(np.ceil(1.5*(self.design['Pupil']['N']/125.)**2*(self.design['Image']['Nlam']/3.)**3))
            else:
                time_est_hrs = int(np.ceil(3*(self.design['Pupil']['N']/125.)**2*(self.design['Image']['Nlam']/3.)**3))
            if time_est_hrs > 12:
                set_queue = """
                #SBATCH --qos=long
                #SBATCH --time={0:02d}:00:00
                """.format(np.min([24, time_est_hrs]))
            else:
                set_queue = """
                #SBATCH --qos=allnccs
                #SBATCH --time={0:02d}:00:00
                """.format(time_est_hrs)
        elif queue_spec is '1h':
            set_queue = """
            #SBATCH --qos=debug
            #SBATCH --time=1:00:00
            """
        elif queue_spec is '12h':
            set_queue = """
            #SBATCH --qos=allnccs
            #SBATCH --time=12:00:00
            """
        else:
            set_queue = """
            #SBATCH --qos=long
            #SBATCH --time=24:00:00
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
        
        exit 0
        """.format(self.fileorg['ampl src fname'])

        bash_fobj.write( textwrap.dedent(header) )
        bash_fobj.write( textwrap.dedent(set_job) )
        bash_fobj.write( textwrap.dedent(set_node) )
        bash_fobj.write( textwrap.dedent(set_queue) )
        bash_fobj.write( textwrap.dedent(intel_module) )
        bash_fobj.write( textwrap.dedent(monitor_mem) )
        bash_fobj.write( textwrap.dedent(call_ampl) )

        bash_fobj.close()
        if verbose:
            logging.info("Wrote %s"%self.fileorg['slurm fname'])
        return 0

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

        if 'slurm fname' not in self.fileorg or self.fileorg['slurm fname'] is None:
            exec_script_fname_tail = self.fileorg['job name'] + ".sh"
            self.fileorg['slurm fname'] = os.path.join(self.fileorg['slurm dir'], exec_script_fname_tail)

        if 'log fname' not in self.fileorg or self.fileorg['log fname'] is None:
            log_fname_tail = self.fileorg['job name'] + ".log"
            self.fileorg['log fname'] = os.path.join(self.fileorg['log dir'], log_fname_tail)
 
        if 'TelAp fname' not in self.fileorg or self.fileorg['TelAp fname'] is None:
            self.fileorg['TelAp fname'] = os.path.join( self.fileorg['TelAp dir'], ("TelAp_quart_" + self.telap_descrip + ".dat") )
 
        if 'FPM fname' not in self.fileorg or self.fileorg['FPM fname'] is None:
            self.fileorg['FPM fname'] = os.path.join( self.fileorg['FPM dir'], "FPM_quart_occspot_M{:03}.dat".format(self.design['FPM']['M']) )

        if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_quart_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}_N{8:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['prim'], self.design['Pupil']['secobs'],
                                                         self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'],
                                                         self.design['Pupil']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_quart_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}Pad{5:02d}_N{6:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['secobs'], self.design['Pupil']['thick'],
                                                         self.design['LS']['pad'], self.design['Pupil']['N'])) )
            else:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_quart_" + \
                                                         "{0:s}{1:02d}D{2:02d}_clear_N{3:04d}.dat".format(self.design['LS']['shape'],
                                                         self.design['LS']['id'], self.design['LS']['od'], self.design['Pupil']['N'])) )

        if self.design['LS']['aligntol'] is not None and ('LDZ fname' not in self.fileorg or self.fileorg['LDZ fname'] is None):
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_quart_" + \
                                                          "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}_Tol{8:02d}_N{9:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['Pupil']['prim'], self.design['Pupil']['secobs'],
                                                          self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'],
                                                          self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_quart_" + \
                                                          "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}Pad{5:02d}_Tol{6:02d}_N{7:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['Pupil']['secobs'], self.design['Pupil']['thick'],
                                                          self.design['LS']['pad'], self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
            else:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_quart_" + \
                                                          "{0:s}{1:02d}D{2:02d}_clear_Tol{3:02d}_N{4:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['LS']['aligntol'], self.design['Pupil']['N'])) )
        self.check_ampl_input_files()
    def write_ampl(self, overwrite=False, override_infile_status=False, ampl_src_fname=None, verbose=True):
        if self.ampl_infile_status is False and not override_infile_status:
            if verbose:
                logging.warning("Error: the most recent input file check for this design configuration failed.")
                logging.warning("The override_infile_status switch is off, so write_ampl() will now abort.")
                logging.warning("See previous warnings in the log to see what file was missing during the initialization")
            return 2
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
                if verbose:
                    logging.warning("Error: {0} already exists and overwrite switch is off, so write_ampl() will now abort".format(self.fileorg['ampl src fname']))
                return 1
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

        if self.design['Pupil']['edge'] is 'floor': # floor to binary
            define_pupil_and_telap = """
            #---------------------

            set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] == 1} (x,y);
            param TelApProp {x in Xs, y in Ys};
            let {x in Xs, y in Ys} TelApProp[x,y] := 0;
            let {(x,y) in Pupil} TelApProp[x,y] := 1;
            """
        elif self.design['Pupil']['edge'] is 'round': # round to binary
            define_pupil_and_telap = """
            #---------------------

            set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] > 0.5} (x,y);
            param TelApProp {x in Xs, y in Ys};
            let {x in Xs, y in Ys} TelApProp[x,y] := 0;
            let {(x,y) in Pupil} TelApProp[x,y] := 1;
            """
        else: # gray, default
            define_pupil_and_telap = """
            #---------------------

            set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] > 0} (x,y);
            param TelApProp {x in Xs, y in Ys};
            let {x in Xs, y in Ys} TelApProp[x,y] := 0;
            let {(x,y) in Pupil} TelApProp[x,y] := TelAp[x,y];
            """

        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon'] is not None: 
            sets_and_arrays = """
            set Mask := setof {mx in MXs, my in MYs: FPM[mx,my] > 0} (mx,my);
            set Lyot := setof {x in Xs, y in Ys: LS[x,y] > 0} (x,y);
            set LyotDarkZone := setof {x in Xs, y in Ys: LDZ[x,y] == 1 && TelApProp[x,y] > 0} (x,y);

            param TR := sum {(x,y) in Pupil} TelApProp[x,y]*dx*dy; # Transmission of the Pupil. Used for calibration.
            
            var A {x in Xs, y in Ys} >= 0, <= 1, := 0.5;
            
            #---------------------
            set DarkHole := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= rho0 && sqrt(xi^2+eta^2) <= rho1} (xi,eta);
            """
        else:
            sets_and_arrays = """
            set Mask := setof {mx in MXs, my in MYs: FPM[mx,my] > 0} (mx,my);
            set Lyot := setof {x in Xs, y in Ys: LS[x,y] > 0} (x,y);

            param TR := sum {(x,y) in Pupil} TelApProp[x,y]*dx*dy; # Transmission of the Pupil. Used for calibration.
            
            var A {x in Xs, y in Ys} >= 0, <= 1, := 0.5;
            
            #---------------------
            set DarkHole := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= rho0 && sqrt(xi^2+eta^2) <= rho1} (xi,eta);
            """
 
        if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon'] is not None: 
            field_propagation = """
            #---------------------
            var EBm_real_X {mx in MXs, y in Ys, lam in Ls};
            var EBm_real {mx in MXs, my in MYs, lam in Ls};
            
            subject to st_EBm_real_X {mx in MXs, y in Ys, lam in Ls}: EBm_real_X[mx,y,lam] = 2*sum {x in Xs: (x,y) in Pupil} TelApProp[x,y]*A[x,y]*cos(2*pi*x*mx/lam)*dx;
            subject to st_EBm_real {(mx, my) in Mask, lam in Ls}: EBm_real[mx,my,lam] = 2/lam*sum {y in Ys} EBm_real_X[mx,y,lam]*cos(2*pi*y*my/lam)*dy;
            
            #---------------------
            var ECm_real_X {x in Xs, my in MYs, lam in Ls};
            var EC_real {x in Xs, y in Ys, lam in Ls};
            
            subject to st_ECm_real_X {x in Xs, my in MYs, lam in Ls}: ECm_real_X[x,my,lam] = 2*sum {mx in MXs: (mx,my) in Mask} FPM[mx,my]*EBm_real[mx,my,lam]*cos(2*pi*x*mx/lam)*dmx;
            subject to st_EC_real {(x,y) in Lyot union LyotDarkZone, lam in Ls}: EC_real[x,y,lam] = TelApProp[x,y]*A[x,y] - 2/lam*sum {my in MYs} ECm_real_X[x,my,lam]*cos(2*pi*y*my/lam)*dmy;
            
            #---------------------
            var ED_real_X {xi in Xis, y in Ys, lam in Ls};
            var ED_real {xi in Xis, eta in Etas, lam in Ls};
            
            subject to st_ED_real_X {xi in Xis, y in Ys, lam in Ls}: ED_real_X[xi,y,lam] = 2*sum {x in Xs: (x,y) in Lyot} EC_real[x,y,lam]*cos(2*pi*x*xi/lam)*dx;
            subject to st_ED_real {(xi, eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] = 2/lam*sum {y in Ys} ED_real_X[xi,y,lam]*cos(2*pi*y*eta/lam)*dy;
            
            #---------------------
            var ED00_real := 0.0;
            subject to st_ED00_real: ED00_real = 4*sum {x in Xs, y in Ys: (x,y) in Lyot} (A[x,y]*TelApProp[x,y])*dx*dy;
            """
        else:
            field_propagation = """
            #---------------------
            var EBm_real_X {mx in MXs, y in Ys, lam in Ls};
            var EBm_real {mx in MXs, my in MYs, lam in Ls};
            
            subject to st_EBm_real_X {mx in MXs, y in Ys, lam in Ls}: EBm_real_X[mx,y,lam] = 2*sum {x in Xs: (x,y) in Pupil} TelApProp[x,y]*A[x,y]*cos(2*pi*x*mx/lam)*dx;
            subject to st_EBm_real {(mx, my) in Mask, lam in Ls}: EBm_real[mx,my,lam] = 2/lam*sum {y in Ys} EBm_real_X[mx,y,lam]*cos(2*pi*y*my/lam)*dy;
            
            #---------------------
            var ECm_real_X {x in Xs, my in MYs, lam in Ls};
            var ECm_real {x in Xs, y in Ys, lam in Ls};
            
            subject to st_ECm_real_X {x in Xs, my in MYs, lam in Ls}: ECm_real_X[x,my,lam] = 2*sum {mx in MXs: (mx,my) in Mask} FPM[mx,my]*EBm_real[mx,my,lam]*cos(2*pi*x*mx/lam)*dmx;
            subject to st_ECm_real {(x,y) in Lyot, lam in Ls}: ECm_real[x,y,lam] = 2/lam*sum {my in MYs} ECm_real_X[x,my,lam]*cos(2*pi*y*my/lam)*dmy;
            
            #---------------------
            var ED_real_X {xi in Xis, y in Ys, lam in Ls};
            var ED_real {xi in Xis, eta in Etas, lam in Ls};
            
            subject to st_ED_real_X {xi in Xis, y in Ys, lam in Ls}: ED_real_X[xi,y,lam] = 2*sum {x in Xs: (x,y) in Lyot} (TelApProp[x,y]*A[x,y]-ECm_real[x,y,lam])*cos(2*pi*x*xi/lam)*dx;
            subject to st_ED_real {(xi, eta) in DarkHole, lam in Ls}: ED_real[xi,eta,lam] = 2/lam*sum {y in Ys} ED_real_X[xi,y,lam]*cos(2*pi*y*eta/lam)*dy;
            
            #---------------------
            var ED00_real := 0.0;
            subject to st_ED00_real: ED00_real = 4*sum {x in Xs, y in Ys: (x,y) in Lyot} (A[x,y]*TelApProp[x,y])*dx*dy;
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
            gurobi_opt_str += " lpmethod=2"
            if self.solver['convtol'] is not None:
                gurobi_opt_str += " barconvtol={0:.1e}".format(np.power(10,-self.solver['convtol']))
            if self.solver['method'] is 'barhom':
                gurobi_opt_str += " barhomogeneous=1"
            if self.solver['crossover'] is True:
                gurobi_opt_str += " crossoverbasis=1"
            else:
                gurobi_opt_str += " crossover=0"
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

        param A_fin {{x in Xs, y in Ys}};
        let {{x in Xs, y in Ys}} A_fin[x,y] := 0;
        let {{(x,y) in Pupil}} A_fin[x,y] := A[x,y];
 
        printf {{y in Ys, x in Xs}}: "%15g %15g %15g \\n", x, y, A_fin[x,y] > "{0:s}";
        """.format(self.fileorg['sol fname'])
 
        mod_fobj.write( textwrap.dedent(header) )
        mod_fobj.write( textwrap.dedent(params) )
        mod_fobj.write( textwrap.dedent(define_coords) )
        mod_fobj.write( textwrap.dedent(load_masks) )
        mod_fobj.write( textwrap.dedent(define_wavelengths) )
        mod_fobj.write( textwrap.dedent(define_pupil_and_telap) )
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
        return 0

    def write_slurm_script(self, queue_spec='auto', account='s1649', email=None, arch=None, overwrite=False, verbose=True):
        if os.path.exists(self.fileorg['slurm fname']):
            if overwrite == True:
                if verbose:
                    logging.warning("Warning: Overwriting the existing copy of {0}".format(self.fileorg['slurm fname']))
            else:
                if verbose:
                    logging.warning("Error: {0} already exists and overwrite switch is off, so write_slurm_script() will now abort".format(self.fileorg['slurm fname']))
                return 1
        elif not os.path.exists(self.fileorg['slurm dir']):
            os.mkdir(self.fileorg['slurm dir'])
            if verbose:
                logging.info("Created new slurm script directory, {0:s}".format(self.fileorg['slurm dir']))

        bash_fobj = open(self.fileorg['slurm fname'], "w") 

        if email is not None: 
            header = """\
            #!/bin/bash

            #PBS -V
            #PBS -m e -M {0:s}
            """.format(email)
        else:
            header = """\
            #! /bin/bash

            """

        set_job = """\
        #SBATCH --job-name={0:s}
        #SBATCH -o {1:s}
        #SBATCH --account={2:s}
        """.format(self.fileorg['job name'], self.fileorg['log fname'], account)
       
        if arch is not None: # can be 'hasw' for Haswell only 
            set_node = """\
            #SBATCH --constraint={0:s}
            #SBATCH --ntasks=1 --nodes=1
            """.format(arch)
        else:
            set_node = """\
            #SBATCH --ntasks=1 --nodes=1
            """.format(arch)

        if queue_spec is 'auto':
            if self.design['LS']['aligntol'] is None:
                time_est_hrs = int(np.ceil(1.5*(self.design['Pupil']['N']/125.)**2*(self.design['Image']['Nlam']/3.)**3))
            else:
                time_est_hrs = int(np.ceil(3*(self.design['Pupil']['N']/125.)**2*(self.design['Image']['Nlam']/3.)**3))
            if time_est_hrs > 12:
                set_queue = """
                #SBATCH --qos=long
                #SBATCH --time={0:02d}:00:00
                """.format(np.min([24, time_est_hrs]))
            else:
                set_queue = """
                #SBATCH --qos=allnccs
                #SBATCH --time={0:02d}:00:00
                """.format(time_est_hrs)
        elif queue_spec is '1h':
            set_queue = """
            #SBATCH --qos=debug
            #SBATCH --time=1:00:00
            """
        elif queue_spec is '12h':
            set_queue = """
            #SBATCH --qos=allnccs
            #SBATCH --time=12:00:00
            """
        else:
            set_queue = """
            #SBATCH --qos=long
            #SBATCH --time=24:00:00
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
        
        exit 0
        """.format(self.fileorg['ampl src fname'])

        bash_fobj.write( textwrap.dedent(header) )
        bash_fobj.write( textwrap.dedent(set_job) )
        bash_fobj.write( textwrap.dedent(set_node) )
        bash_fobj.write( textwrap.dedent(set_queue) )
        bash_fobj.write( textwrap.dedent(intel_module) )
        bash_fobj.write( textwrap.dedent(monitor_mem) )
        bash_fobj.write( textwrap.dedent(call_ampl) )

        bash_fobj.close()
        if verbose:
            logging.info("Wrote %s"%self.fileorg['slurm fname'])
        return 0
