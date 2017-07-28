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
import pyfits
import matplotlib
matplotlib.use('Agg') # non-interactive
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches
matplotlib.rcParams['image.origin'] = 'lower'
matplotlib.rcParams['image.interpolation'] = 'nearest'
matplotlib.rcParams['image.cmap'] = 'gray'
matplotlib.rcParams['axes.linewidth'] = 1.
matplotlib.rcParams['lines.linewidth'] = 2.5
matplotlib.rcParams['font.size'] = 12

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
        if 'M' in design_params['FPM'] and isinstance(coron, SPLC): # M is only defined implicitly in the SPLC class
            design_params['FPM'].pop('M',None)
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

def merge_design_param_surveys(survey_list, merged_survey_fname=None):
    assert len(survey_list) >= 1, "The input list must contain at least one survey."
    for survey in survey_list:
        assert isinstance(survey, DesignParamSurvey), "Each element in the input list must be a DesignParamSurvey."
        assert survey.coron_class == survey_list[0].coron_class, "The coronagraph classes of the input surveys must match."
        # Until we take the time to write out a more intelligent merge routine,
        # require the fixed and varied parameter categories to match
        # This means, if there are parameters that are constant within an individual survey
        # but change across the surveys, they must be defined as single-element lists so that they
        # are classified as 'varied' in the input survey objects.
        assert survey.fixed_param_index == survey_list[0].fixed_param_index, "Fixed parameter categories of each survey must be the same."
        assert survey.fixed_param_vals == survey_list[0].fixed_param_vals, "Fixed parameter values of each survey must be the same."
        assert survey.varied_param_index == survey_list[0].varied_param_index, "Varied parameter categories of each survey must be the same."
 
    merged_survey = DesignParamSurvey(coron_class=survey_list[0].coron_class, survey_config=survey_list[0].survey_config,
                                      fileorg=survey_list[0].fileorg, solver=survey_list[0].solver)
    varied_param_combos_list_merged = list(merged_survey.varied_param_combos)

    for survey in survey_list[1:]:
        varied_param_combos_list_B = list(survey.varied_param_combos)
        for ii, combo in enumerate(varied_param_combos_list_B):
            if combo not in varied_param_combos_list_merged:
                varied_param_combos_list_merged.append(combo)
                merged_survey.coron_list.append(survey.coron_list[ii]) 
        
    merged_survey.varied_param_combos = tuple(varied_param_combos_list_merged)
    merged_survey.N_combos = len(merged_survey.varied_param_combos)
    if merged_survey_fname is not None: # unless specified, use the same name as the first input survey
        merged_survey.fileorg['survey fname'] = merged_survey_fname

    merged_survey_name = os.path.basename(merged_survey.fileorg['survey fname'][:-4])
    num_ID_digits = int(np.floor(np.log10(merged_survey.N_combos))) + 1
    ID_fmt_str = "{{:s}}-{{:0{:d}d}}".format(num_ID_digits)
    for idx, coron in enumerate(merged_survey.coron_list): # overwrite design IDs in coron_list
        coron.fileorg['design ID'] = ID_fmt_str.format(merged_survey_name, idx)
    
    return merged_survey

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
        if 'planeofconstr' not in self.solver or self.solver['planeofconstr'] is None: self.solver['planeofconstr'] = 'FP2'
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
        for idx, param_combo in enumerate(self.varied_param_combos): # TODO: Switch the coronagraph type depending on the symmetry of the telescope aperture and support struts 
            for (varied_keycat, varied_parname), current_val in zip(self.varied_param_index, param_combo):
                design[varied_keycat][varied_parname] = current_val
            coron_fileorg = self.fileorg.copy()
            if 'survey fname' in coron_fileorg:
                survey_name = os.path.basename(coron_fileorg['survey fname'][:-4])
                coron_fileorg.pop('survey fname')
            else:
                survey_name = os.path.basename(os.path.abspath(self.fileorg['work dir']))
            num_ID_digits = int(np.floor(np.log10(self.N_combos))) + 1
            ID_fmt_str = "{{:s}}-{{:0{:d}d}}".format(num_ID_digits)
            coron_fileorg['design ID'] = ID_fmt_str.format(survey_name, idx)
            self.coron_list.append( coron_class(design=design, fileorg=coron_fileorg, solver=self.solver) )
 
        setattr(self, 'ampl_infile_status', False)
        self.check_ampl_input_files()
        setattr(self, 'ampl_src_status', False)
        setattr(self, 'ampl_submission_status', False)
        setattr(self, 'solution_status', False)
        setattr(self, 'eval_status', False)

    def write_serial_bash(self, serial_bash_fname=None, overwrite=False, override_infile_status=False):
        # Write a bash script to sequentially run each program in a design survey
        if serial_bash_fname is None:
            coron_fileorg = self.fileorg.copy()
            if 'survey fname' in coron_fileorg:
                survey_name = os.path.basename(coron_fileorg['survey fname'][:-4])
            else:
                survey_name = os.path.basename(os.path.abspath(coron_fileorg['work dir']))
            serial_bash_fname = os.path.join(coron_fileorg['work dir'], "run_" + survey_name + "_serial.sh")
        else:
            basename = os.path.basename(serial_bash_fname)
            serial_bash_fname = os.path.join(coron_fileorg['work dir'], basename)

        if self.ampl_infile_status is False and not override_infile_status:
            logging.warning("Error: the most recent input file check for this survey configuration failed.")
            logging.warning("The override_infile_status switch is off, so write_serial_bash() will now abort.")
            return 2
        if not os.path.exists(serial_bash_fname) or overwrite:
            serial_bash_fobj = open(serial_bash_fname, "w")
            serial_bash_fobj.write("#! /bin/bash -x\n")
            for coron in self.coron_list:
                serial_bash_fobj.write("ampl {0:s} > {1:s}\n".format(coron.fileorg['ampl src fname'],
                                                                     coron.fileorg['log fname']))
            serial_bash_fobj.close()
            os.chmod(serial_bash_fname, 0775)
            logging.info("Wrote serial bash survey script to {:s}".format(serial_bash_fname))
            return serial_bash_fname
        else:
            logging.warning("Denied overwrite of serial bash survey script {:s}".format(serial_bash_fname))
            return 1

    def write_ampl_batch(self, overwrite=False, override_infile_status=False):
        write_count = 0
        overwrite_deny_count = 0
        infile_deny_count = 0
        for coron in self.coron_list:
            status = coron.write_ampl(overwrite, override_infile_status, verbose=False)
            if status == 2:
                infile_deny_count += 1
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

    def describe(self):
        print("This survey has {0:d} design parameter combinations.".format(self.N_combos))
        print("{0:d} parameters are varied: {1}".format(len(self.varied_param_index), self.varied_param_index))
        print("")
        print("File organization:")
        pprint.pprint(self.fileorg)
        print("")
        print("All input files exist? {}".format(self.check_ampl_input_files()))
        print("")
        print("Last coronagraph in survey list:")
        print("Telescope aperture file {:s}".format(self.coron_list[-1].fileorg['TelAp fname']))
        print("Focal plane mask file {:s}".format(self.coron_list[-1].fileorg['FPM fname']))
        print("Lyot stop file {:s}".format(self.coron_list[-1].fileorg['LS fname']))
        print("Job label {:s}".format(self.coron_list[-1].fileorg['job name']))
        print("Varied parameter combo tuple:")
        pprint.pprint(self.varied_param_combos[-1])

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
        telap_warning = False
        for coron in self.coron_list:
            if os.path.exists(coron.fileorg['sol fname']) and \
                (coron.eval_metrics['fwhm area'] is None \
                 or coron.eval_metrics['apod nb res ratio'] is None):
                telap_flag = coron.get_metrics(verbose=verbose)
                coron.eval_status = True
                if telap_flag > 0:
                    telap_warning = True
            if os.path.exists(coron.fileorg['sol fname']) and os.path.exists(coron.fileorg['log fname']) and \
               (not hasattr(coron, 'ampl_completion_time') or coron.ampl_completion_time is None):
                setattr(coron, 'ampl_completion_time', None)
                log = open(coron.fileorg['log fname'])
                lines = log.readlines()
                for line in lines:
                    if 'iterations' in line and 'seconds' in line:
                        split_line = line.split()
                        coron.ampl_completion_time = float(split_line[split_line.index('seconds')-1])/3600
                        break
        if telap_warning:
            logging.warning("No unpadded version of telescope aperture was found, so the optimization version was used to derive throughput metrics.")

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
                catrow.extend(['Design ID', 'AMPL program', '', '', '', 'Solution', '', 'Evaluation metrics', '', ''])
                paramrow.extend(['survey-index', 'src filename', 'src exists?', 'input files?', 'submitted?', 'sol filename', 'sol exists?', 'comp time (h)',
                                 'inc. energy', 'apodizer non-binarity', 'Tot thrupt', 'half-max thrupt', 'half-max circ thrupt', 'rel. half-max thrupt', 'r=0.7 thrupt',  'r=0.7 circ thrupt', 'rel. r=0.7 thrupt', 'PSF area'])
                surveywriter.writerow(catrow)
                surveywriter.writerow(paramrow)
                for ii, param_combo in enumerate(self.varied_param_combos):
                    param_combo_row = list(param_combo)
                    param_combo_row.append(self.coron_list[ii].fileorg['design ID'])
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
                    if hasattr(self.coron_list[ii], 'ampl_completion_time') and self.coron_list[ii].ampl_completion_time is not None:
                        param_combo_row.append(self.coron_list[ii].ampl_completion_time)
                    else:
                        param_combo_row.append('')
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
                    if self.coron_list[ii].eval_metrics['fwhm circ thrupt'] is not None:
                        param_combo_row.append(self.coron_list[ii].eval_metrics['fwhm circ thrupt'])
                    else:
                        param_combo_row.append('')
                    if self.coron_list[ii].eval_metrics['rel fwhm thrupt'] is not None:
                        param_combo_row.append(self.coron_list[ii].eval_metrics['rel fwhm thrupt'])
                    else:
                        param_combo_row.append('')
                    if self.coron_list[ii].eval_metrics['p7ap thrupt'] is not None:
                        param_combo_row.append(self.coron_list[ii].eval_metrics['p7ap thrupt'])
                    else:
                        param_combo_row.append('')
                    if self.coron_list[ii].eval_metrics['p7ap circ thrupt'] is not None:
                        param_combo_row.append(self.coron_list[ii].eval_metrics['p7ap circ thrupt'])
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
                                 'sol dir', 'log dir', 'eval dir', 'eval subdir', 'slurm dir',
                                 'ampl src fname', 'slurm fname', 'log fname', 'job name', 'design ID', 
                                 'TelAp fname', 'FPM fname', 'LS fname', 'LDZ fname', 'sol fname'],
                     'solver': ['planeofconstr', 'constr', 'method', 'presolve', 'threads', 'solver', 'crossover'] }

    _solver_menu = { 'planeofconstr': ['FP1', 'Lyot', 'FP2'], 
                     'constr': ['lin', 'quad'], 'solver': ['LOQO', 'gurobi', 'gurobix'], 
                     'method': ['bar', 'barhom', 'dualsimp'],
                     'presolve': [True, False], 'threads': [None]+range(1,33), 'crossover': [None]+[True, False] }

    _aperture_menu = { 'prim': ['hex1', 'hex2', 'hex3', 'hex4', 'key24', 'pie12', 'pie08', 'circ', 
                                'ochex1', 'ochex2', 'ochex3', 'ochex4', 'irisao', 'atlast',
                                'luvoir15m', 'luvoir15mwStruts'],
                       'secobs': ['Y60d','Yoff60d','X','Cross','T','Y90d', 'Y00d'],
                       'thick': ['025','100','125'],
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
            if namekey.endswith(' dir') and ( namekey not in self.fileorg or self.fileorg[namekey] is None ):
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
        if 'planeofconstr' not in self.solver or self.solver['planeofconstr'] is None: self.solver['planeofconstr'] = 'FP2'
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
        setattr(self, 'solution_status', False) # Only changed by the queue filler program
        setattr(self, 'ampl_completion_time', None)

        setattr(self, 'eval_metrics', {})
        self.eval_metrics['inc energy'] = None
        self.eval_metrics['tot thrupt'] = None
        self.eval_metrics['fwhm thrupt'] = None
        self.eval_metrics['fwhm circ thrupt'] = None
        self.eval_metrics['p7ap thrupt'] = None
        self.eval_metrics['p7ap circ thrupt'] = None
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

    def get_design_portrait(self, intens_maps, intens_curves, xis, seps, star_diams,
                            second_curve_diam=None, use_gray_gap_zero=False):
        TelAp, Apod, FPM, LS = self.get_coron_masks(use_gray_gap_zero=use_gray_gap_zero)
        portrait_fig = plt.figure(figsize=(10,7))
        gs1 = gridspec.GridSpec(2, 3)
        gs1.update(left=0.01, right=0.99, bottom=0.02, top=0.99, wspace=0.01)
        ax1 = plt.subplot(gs1[0,0])
        plt.imshow(TelAp)
        _ = plt.axis('off')
        ax2 = plt.subplot(gs1[0,1])
        plt.imshow(Apod)
        _ = plt.axis('off')
        ax3 = plt.subplot(gs1[0,2])
        #plt.imshow(LS.T) # why was this previously transposed?
        plt.imshow(LS)
        _ =plt.axis('off')
        gs2 = gridspec.GridSpec(2, 3)
        gs2.update(left=0.03, right=1.36, bottom=0.04, top=1.10, wspace=0.01)
        ax4 = plt.subplot(gs2[1,0])
        pixscale_lamoD = xis[1] - xis[0]
        tick_labels = np.arange(np.round(xis[0]), np.round(xis[-1]+pixscale_lamoD), 2.)
        xc_pix = intens_maps.shape[-1]/2 - 0.5
        if isinstance(self, NdiayeAPLC):
            fpm_rad = self.design['FPM']['rad']
        else:
            fpm_rad = self.design['FPM']['R0']
        fpm_rad_pix = fpm_rad/pixscale_lamoD
        tick_locs = tick_labels/pixscale_lamoD + xc_pix
        plt.imshow(np.log10(intens_maps[0,:,:]),
                   vmin=-(self.design['Image']['c']+2),
                   vmax=-(self.design['Image']['c']-1), cmap='CMRmap')
        fpm_circle = matplotlib.patches.Circle((xc_pix, xc_pix), fpm_rad_pix, facecolor='none',
                                               edgecolor='w', linewidth=2., alpha=1.,
                                               clip_on=False, linestyle='--')
        ax4.add_patch(fpm_circle)
        plt.xticks(tick_locs, tick_labels)
        plt.yticks(tick_locs, tick_labels)
        plt.tick_params(labelsize=9)
        plt.colorbar(orientation='vertical', shrink=0.75, pad=0.03)
        gs3 = gridspec.GridSpec(2, 3)
        gs3.update(left=0.36, right=0.98, bottom=0.08, top=1.05, wspace=0.01)
        ax5 = plt.subplot(gs3[1,1:])
        diam_1 = star_diams[0]
        diam_1_curve, = plt.semilogy(seps, intens_curves[0,:], 'b', zorder=3)
        if second_curve_diam is not None:
            ind_diam_2 = np.argmin(np.abs(np.array(star_diams) - second_curve_diam))
            diam_2 = star_diams[ind_diam_2]
            diam_2_curve, = plt.semilogy(seps, intens_curves[ind_diam_2,:], 'r', zorder=2)
        fpm_line = plt.vlines(fpm_rad, 10**-16, 1, linestyle='--', color='gray', zorder=1)
        plt.xlim([seps[0], seps[-1]])
        plt.ylim([10**-(self.design['Image']['c']+1), 5*10**-(self.design['Image']['c'])])
        plt.ylabel(r'$I/I_\star$',fontsize=16)
        plt.xlabel(r'Separation ($\lambda/D$)',fontsize=12)
        if second_curve_diam is not None:
            plt.legend([diam_1_curve, diam_2_curve, fpm_line],
                       [r'$\theta_\star=${:.2f} $\lambda/D$'.format(diam_1),
                        r'$\theta_\star=${:.2f} $\lambda/D$'.format(diam_2),
                        'FPM radius'], fontsize=12, loc='upper center')
        else:
            plt.legend([diam_1_curve, fpm_line],
                       [r'$\theta_\star=${:.2f} $\lambda/D$'.format(diam_1),
                        'FPM radius'], fontsize=12, loc='upper center')
        return portrait_fig

    def write_design_package(self, eval_path=None, pixscale_lamoD=0.25, Nlam=None, dpi=300):
        """
		Write a coronagraph design package including mask files (Telescope pupil,
        Apodizer, FPM, Lyot stop) in FITS format and a simple portrait viewgraph.
        TBA: example PSF evaluation scripts.
        """
        if eval_path is None:
            if 'design ID' in self.fileorg:
                design_label = "{:s}_{:s}".format(self.fileorg['design ID'],
                                                  self.fileorg['job name'])
            else:
                design_label = self.fileorg['job name']
            self.fileorg['eval subdir'] = os.path.join(self.fileorg['eval dir'], design_label)
        else:
            self.fileorg['eval subdir'] = os.path.normpath(eval_path)
            design_label = os.path.basename(eval_path)
        eval_path = self.fileorg['eval subdir']
        if not os.path.exists(eval_path):
            os.mkdir(eval_path)
        if Nlam is None:
            Nlam = 2*self.design['Image']['Nlam'] + 1

        TelAp, Apod, FPM, LS = self.get_coron_masks(use_gray_gap_zero=False)
        telap_hdu = pyfits.PrimaryHDU(TelAp)
        telap_fits_fname = os.path.join(eval_path, 'TelAp.fits')
        telap_hdu.writeto(telap_fits_fname, clobber=True)
        apod_hdu = pyfits.PrimaryHDU(Apod)
        apod_fits_fname = os.path.join(eval_path, 'Apod.fits')
        apod_hdu.writeto(apod_fits_fname, clobber=True)
        fpm_hdu = pyfits.PrimaryHDU(FPM)
        fpm_fits_fname = os.path.join(eval_path, 'FPM.fits')
        fpm_hdu.writeto(fpm_fits_fname, clobber=True)
        LS_hdu = pyfits.PrimaryHDU(LS)
        LS_fits_fname = os.path.join(eval_path, 'LS.fits')
        LS_hdu.writeto(LS_fits_fname, clobber=True)

        if isinstance(self, NdiayeAPLC):
            xis, intens_polychrom, seps, radial_intens_polychrom = self.get_onax_psf(Nlam=Nlam)
        elif isinstance(self, SPLC):
            xis, intens_polychrom, seps, radial_intens_polychrom, fov_mask = self.get_onax_psf(Nlam=Nlam)
        else:
            logging.error("Unsupported coronagraph class {}".format(self.__class__))
            return 1

        portrait_fig = \
          self.get_design_portrait(np.mean(intens_polychrom, axis=0).reshape(
                                    (1,intens_polychrom.shape[1],intens_polychrom.shape[2])),
                                   np.mean(radial_intens_polychrom, axis=0).reshape((1,len(seps))),
                                   np.array(xis.T), seps, [0.])
        portrait_fname = os.path.join(eval_path, 'DesignPortrait_simple_{:s}.png'.format(design_label))
        portrait_fig.savefig(portrait_fname, dpi=dpi)

        return eval_path

    def write_eval_products(self, pixscale_lamoD=0.25, star_diam_vec=None, Npts_star_diam=7, Nlam=None, 
                            norm='aperture', second_curve_diam=0.2, dpi=300):
        if 'eval subdir' not in self.fileorg or self.fileorg['eval subdir'] is None:
            if 'design ID' in self.fileorg:
                design_label = "{:s}_{:s}".format(self.fileorg['design ID'],
                                                  self.fileorg['job name'])
            else:
                design_label = self.fileorg['job name']
            self.fileorg['eval subdir'] = os.path.join(self.fileorg['eval dir'], design_label)
        if not os.path.exists(self.fileorg['eval dir']):
            os.mkdir(self.fileorg['eval dir'])
        if not os.path.exists(self.fileorg['eval subdir']):
            os.mkdir(self.fileorg['eval subdir'])

        if star_diam_vec is None:
            star_diam_vec = np.concatenate([np.linspace(0,0.09,10), np.linspace(0.1, 1, 10), np.array([2., 3., 4.])])
        if Nlam is None:
            Nlam = 2*self.design['Image']['Nlam'] + 1

        stellar_intens_map, stellar_intens_curves, xis, seps, stellar_intens_diam_vec, \
        offax_psf, offax_psf_offset_vec, sky_trans, contrast_convert_fac = \
          self.get_yield_input_products(pixscale_lamoD, star_diam_vec, Npts_star_diam, Nlam, norm)

        design_portrait_fig = \
          self.get_design_portrait(stellar_intens_map*contrast_convert_fac,
                                   stellar_intens_curves*contrast_convert_fac,
                                   xis, seps, star_diam_vec,
                                   second_curve_diam=second_curve_diam, use_gray_gap_zero=False)
        design_portrait_fname = os.path.join(self.fileorg['eval subdir'], 'DesignPortrait_{:s}.png'.format(self.fileorg['job name']))
        design_portrait_fig.savefig(design_portrait_fname, dpi=dpi)
        logging.info("Wrote design potrait to {:s}".format(design_portrait_fname))
        plt.close(design_portrait_fig)

        xycent_pix_coord = (float(stellar_intens_map.shape[2])/2, 
                            float(stellar_intens_map.shape[1])/2)
        obscure_ratio = (np.pi/4 - self.eval_metrics['inc energy'])/(np.pi/4)

        stellar_intens_fname = os.path.join(self.fileorg['eval subdir'], 'stellar_intens.fits')
        stellar_intens_diam_list_fname = os.path.join(self.fileorg['eval subdir'], 'stellar_intens_diam_list.fits')
        offax_psf_fname = os.path.join(self.fileorg['eval subdir'], 'offax_psf.fits')
        offax_psf_offset_list_fname = os.path.join(self.fileorg['eval subdir'], 'offax_psf_offset_list.fits')
        sky_trans_fname = os.path.join(self.fileorg['eval subdir'], 'sky_trans.fits')

        header = pyfits.Header()
        header['DESIGN'] = self.fileorg['job name']
        header['PIXSCALE'] = (pixscale_lamoD, 'pixel scale in units of lambda0/D')
        header['LAMBDA'] = (1.0, 'central wavelength of the bandpass in microns')
        header['MINLAM'] = (1.0-self.design['Image']['bw']/2, 'shortest wavelength of the bandpass in microns')
        header['MAXLAM'] = (1.0+self.design['Image']['bw']/2, 'shortest wavelength of the bandpass in microns')
        header['XCENTER'] = (xycent_pix_coord[0], 'center of the image in 00LL pixel coordinates')
        header['YCENTER'] = (xycent_pix_coord[1], 'center of the image in 00LL pixel coordinates')
        header['OBSCURED'] = (obscure_ratio, 'fraction of the aperture area that is obscured w.r.t. circular aperture')
        header['JITTER'] = (0, 'RMS jitter per axis in mas')
        header['N_LAM'] = (Nlam, 'number of wavelength samples used in evaluation')
        header['N_STAR'] = (Npts_star_diam, 'number of points across stellar diameter')

        stellar_intens_hdu = pyfits.PrimaryHDU(stellar_intens_map, header=header)
        stellar_intens_hdu.writeto(stellar_intens_fname, clobber=True)
        diam_list_hdu = pyfits.PrimaryHDU(stellar_intens_diam_vec, header=header)
        diam_list_hdu.writeto(stellar_intens_diam_list_fname, clobber=True)
        offax_psf_hdu = pyfits.PrimaryHDU(offax_psf, header=header)
        offax_psf_hdu.writeto(offax_psf_fname, clobber=True)
        offset_list_hdu = pyfits.PrimaryHDU(offax_psf_offset_vec, header=header)
        offset_list_hdu.writeto(offax_psf_offset_list_fname, clobber=True)
        sky_trans_hdu = pyfits.PrimaryHDU(sky_trans, header=header)
        sky_trans_hdu.writeto(sky_trans_fname, clobber=True)
        
        logging.info("Wrote stellar intensity map to {:s}".format(stellar_intens_fname))
        logging.info("Wrote stellar intensity diameter list to {:s}".format(stellar_intens_diam_list_fname))
        logging.info("Wrote off-axis PSF to {:s}".format(offax_psf_fname))
        logging.info("Wrote off-axis PSF offset list to {:s}".format(offax_psf_offset_list_fname))
        logging.info("Wrote sky transmission map to {:s}".format(sky_trans_fname))

class SPLC(LyotCoronagraph): # SPLC following Zimmerman et al. (2016), uses diaphragm FPM 
    _design_fields = OrderedDict([ ( 'Pupil', OrderedDict([('N',(int, 125)), ('prim',(str, 'hex3')), ('secobs',(str, 'X')), 
                                                           ('thick',(str, '025')), ('centobs',(int, 1)),
                                                           ('gap',(int, 1)), ('edge',(str, 'gray'))]) ),
                                   ( 'FPM', OrderedDict([('R0',(float, 4.)), ('R1',(float, 10.)), ('openang',(int, 180)),
                                                         ('orient',(str, 'H')), ('fpmres',(int, 10))]) ),
                                   ( 'LS', OrderedDict([('N',(int, 125)), ('shape',(str, 'ann')), ('id',(int, 25)), ('od',(int, 75)),
                                                        ('obscure',(int, 0)), ('pad',(int, 0)),
                                                        ('aligntol',(int, None)), ('aligntolcon',(float, 3.))]) ),
                                   ( 'Image', OrderedDict([('c',(float, 10.)), ('bw',(float, 0.10)), ('Nlam',(int, 1)),
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
            self.design['Image']['Nlam'] = int(np.ceil(self.design['Image']['bw+']/(0.12/4)))
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
        if self.design['Pupil']['prim'] is 'wfirst':
            self.telap_descrip = "wfirstCobs{0:02d}sthick{1:s}_N{2:04d}".format(self.design['Pupil']['centobs'], self.design['Pupil']['thick'], \
                                                                                self.design['Pupil']['N'])
        elif self.design['Pupil']['prim'] is 'wfirstCycle5':
            self.telap_descrip = "wfirstCycle5_N{0:04d}".format(self.design['Pupil']['N'])
        else:
            self.telap_descrip = "{0:s}{1:s}{2:s}cobs{3:d}gap{4:d}_N{5:04d}".format(self.design['Pupil']['prim'], self.design['Pupil']['secobs'], self.design['Pupil']['thick'], \
                                                                                    int(self.design['Pupil']['centobs']), self.design['Pupil']['gap'], self.design['Pupil']['N'])
        self.amplname_pupil = "{0:s}{1:s}".format(self.telap_descrip, self.design['Pupil']['edge'][0])

        self.amplname_fpm = "FPM{0:03d}R{1:03d}{2:s}{3:03d}res{4:02d}".format(int(10*self.design['FPM']['R0']), int(10*self.design['FPM']['R1']),
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
                  os.path.join( self.fileorg['FPM dir'], "FPM_full_diaphragm_{0:03d}R{1:03d}M{2:03d}_{3:s}{4:03d}deg.dat".format(
                                                          int(10*self.design['FPM']['fpmres']*self.design['FPM']['R0']),
                                                          int(10*self.design['FPM']['fpmres']*self.design['FPM']['R1']),
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

    def get_coron_masks(self, use_gray_gap_zero=True):
        TelAp_basename = os.path.basename(self.fileorg['TelAp fname'])
        if use_gray_gap_zero:
            gapstr_beg = TelAp_basename.find('gap')
            TelAp_nopad_basename = TelAp_basename.replace(TelAp_basename[gapstr_beg:gapstr_beg+4], 'gap0')
            TelAp_nopad_fname = os.path.join( os.path.dirname(self.fileorg['TelAp fname']), TelAp_nopad_basename )
            TelAp_p = np.loadtxt(TelAp_nopad_fname)
        elif self.design['Pupil']['edge'] == 'floor': # floor to binary
            TelAp_p = np.floor(np.loadtxt(self.fileorg['TelAp fname']))
        else:
            TelAp_p = np.round(np.loadtxt(self.fileorg['TelAp fname']))

        A_col = np.loadtxt(self.fileorg['sol fname'])[:,-1]
        FPM_p = np.loadtxt(self.fileorg['FPM fname'])
        LS_p = np.loadtxt(self.fileorg['LS fname'])
        A_p = A_col.reshape(TelAp_p.shape)
        if isinstance(self, QuarterplaneSPLC):
            TelAp = np.concatenate((np.concatenate((TelAp_p[::-1,::-1], TelAp_p[:,::-1]),axis=0),
                                    np.concatenate((TelAp_p[::-1,:], TelAp_p),axis=0)), axis=1)
            Apod = np.concatenate((np.concatenate((A_p[::-1,::-1], A_p[:,::-1]),axis=0),
                                   np.concatenate((A_p[::-1,:], A_p),axis=0)), axis=1)
            FPM = np.concatenate((np.concatenate((FPM_p[::-1,::-1], FPM_p[:,::-1]),axis=0),
                                  np.concatenate((FPM_p[::-1,:], FPM_p),axis=0)), axis=1)
            LS = np.concatenate((np.concatenate((LS_p[::-1,::-1], LS_p[:,::-1]),axis=0),
                                 np.concatenate((LS_p[::-1,:], LS_p),axis=0)), axis=1)
        elif isinstance(self, HalfplaneSPLC):
            TelAp = np.concatenate((TelAp_p[:,::-1], TelAp_p), axis=1)
            Apod = np.concatenate((A_p[:,::-1], A_p), axis=1)
            FPM = np.concatenate((np.concatenate((FPM_p[::-1,::-1], FPM_p[:,::-1]),axis=0),
                                  np.concatenate((FPM_p[::-1,:], FPM_p),axis=0)), axis=1)
            LS = np.concatenate((LS_p[:,::-1], LS_p), axis=1)
        else:
            TelAp = TelAp_p
            Apod = A_p
            FPM = FPM_p
            LS = LS_p
        return TelAp, Apod, FPM, LS

    def get_coords(self, fp2res=8, rho_out=None, Nlam=None): # for SPLC, with separate Lyot plane sampling
        D = 1.
        N_A = self.design['Pupil']['N']
        N_L = self.design['LS']['N']
        bw = self.design['Image']['bw']
        fpmres = self.design['FPM']['fpmres']
        M_fp1 = self.design['FPM']['M']
        if rho_out is None:
            rho_out = self.design['FPM']['R1'] + 0.5
        if Nlam is None:
            Nlam = self.design['Image']['Nlam']
        M_fp2 = int(np.ceil(rho_out*fp2res))
        
        # pupil plane
        dx = (D/2)/N_A
        dy = dx
        xs = np.matrix(np.linspace(-N_A+0.5,N_A-0.5,2*N_A)*dx)
        XX, YY = np.meshgrid(np.array(xs), np.array(xs))
        
        # FPM
        dmx = 1./fpmres
        mxs = np.matrix(np.linspace(-M_fp1+0.5,M_fp1-0.5,2*M_fp1)*dmx)

        # Lyot plane
        du = (D/2)/N_L
        us = np.matrix(np.linspace(-N_L+0.5,N_L-0.5,2*N_L)*du)
        
        # FP2
        dxi = 1./fp2res
        xis = np.matrix(np.linspace(-M_fp2+0.5,M_fp2-0.5,2*M_fp2)*dxi)

        # wavelength ratios
        wrs = np.linspace(1.-bw/2, 1.+bw/2, Nlam)

        return xs, dx, XX, YY, mxs, dmx, us, du, xis, dxi, wrs

    def get_onax_psf(self, fp2res=8, rho_inc=0.25, rho_out=None, Nlam=None): # for SPLC
        if self.design['Pupil']['edge'] == 'floor': # floor to binary
            TelAp_p = np.floor(np.loadtxt(self.fileorg['TelAp fname'])).astype(int)
        elif self.design['Pupil']['edge'] == 'round': # round to binary
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
            rho_out = self.design['FPM']['R1'] + 0.5
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
                theta_quad_mask = np.greater(theta_quad, self.design['FPM']['openang']/2)
            else:
                theta_quad_mask = np.less(theta_quad, self.design['FPM']['openang']/2)
            theta_rhs_mask = np.concatenate((theta_quad_mask[::-1,:], theta_quad_mask), axis=0)
            theta_mask = np.concatenate((theta_rhs_mask[:,::-1], theta_rhs_mask), axis=1)
            FoV_mask = theta_mask*rad_mask
        else:
            FoV_mask = rad_mask

        for si, sep in enumerate(seps):
            r_in = np.max([seps[0], sep-0.25])
            r_out = np.min([seps[-1], sep+0.25])
            meas_mask = (FoV_mask & (RRs >= r_in) & (RRs <= r_out))
            meas_ann_ind = np.nonzero(np.ravel(meas_mask))[0]
            if len(meas_ann_ind) > 0:
                for wi, wr in enumerate(wrs):
                    radial_intens_polychrom[wi, si] = np.mean(np.ravel(intens_polychrom[wi,:,:])[meas_ann_ind])
            else:
                radial_intens_polychrom[:, si] = np.nan

        return xis, intens_polychrom, seps, radial_intens_polychrom, FoV_mask

    def get_metrics(self, fp1res=8, fp2res=16, rho_out=None, Nlam=None, verbose=True): # for SPLC class
        TelAp_basename = os.path.basename(self.fileorg['TelAp fname'])
        gapstr_beg = TelAp_basename.find('gap')
        TelAp_nopad_basename = TelAp_basename.replace(TelAp_basename[gapstr_beg:gapstr_beg+4], 'gap0')
        TelAp_nopad_fname = os.path.join( os.path.dirname(self.fileorg['TelAp fname']), TelAp_nopad_basename )
        if os.path.exists(TelAp_nopad_fname):
            TelAp_p = np.loadtxt(TelAp_nopad_fname)
            telap_flag = 0
        else:
            TelAp_p = np.loadtxt(self.fileorg['TelAp fname'])
            telap_flag = 1
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
            rho_out = self.design['FPM']['R1']

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

        fwhm_sum_TelAp = np.sum(intens_TelAp[fwhm_ind_TelAp])*dxi*dxi
        fwhm_sum_SPLC = np.sum(intens_D_0[fwhm_ind_SPLC])*dxi*dxi
        p7ap_sum_TelAp = np.sum(intens_TelAp[p7ap_ind])*dxi*dxi
        p7ap_sum_SPLC = np.sum(intens_D_0[p7ap_ind])*dxi*dxi

        self.eval_metrics['inc energy'] = np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['tot thrupt'] = np.sum(intens_D_0*dxi*dxi)/np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['fwhm thrupt'] = fwhm_sum_SPLC/np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['fwhm circ thrupt'] = fwhm_sum_SPLC/(np.pi/4)
        self.eval_metrics['p7ap thrupt'] = p7ap_sum_SPLC/np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['p7ap circ thrupt'] = p7ap_sum_SPLC/(np.pi/4)
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
            print("Band-averaged half-max throughput, circ. ref.: {:.2f}%".format(100*self.eval_metrics['fwhm circ thrupt']))
            print("Band-averaged r=.7 lam/D throughput: {:.2f}%".format(100*self.eval_metrics['p7ap thrupt']))
            print("Band-averaged r=.7 lam/D throughput, circ. ref.: {:.2f}%".format(100*self.eval_metrics['p7ap circ thrupt']))
            print("Band-averaged relative half-max throughput: {:.2f}%".format(100*self.eval_metrics['rel fwhm thrupt']))
            print("Band-averaged relative r=0.7 lam/D throughput: {:.2f}%".format(100*self.eval_metrics['rel p7ap thrupt']))
            print("Band-averaged FWHM PSF area / (lambda0/D)^2: {:.2f}".format(self.eval_metrics['fwhm area']))
        return telap_flag

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
              os.path.join( self.fileorg['FPM dir'], "FPM_quart_diaphragm_{0:03d}R{1:03d}M{2:03d}_{3:s}{4:03d}deg.dat".format(
                                                      int(10*self.design['FPM']['fpmres']*self.design['FPM']['R0']),
                                                      int(10*self.design['FPM']['fpmres']*self.design['FPM']['R1']),
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

        if self.design['Pupil']['edge'] == 'floor': # floor to binary
            define_pupil_and_telap = """
            #---------------------

            set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] == 1} (x,y);
            param TelApProp {x in Xs, y in Ys};
            let {x in Xs, y in Ys} TelApProp[x,y] := 0;
            let {(x,y) in Pupil} TelApProp[x,y] := 1;
            """
        elif self.design['Pupil']['edge'] == 'round': # round to binary
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

        field_propagation_to_FP1 = """
        #---------------------
        var EB_real_X {mx in MXs, y in Ys, lam in Ls};
        var EB_real {mx in MXs, my in MYs, lam in Ls};
        
        subject to st_EB_real_X {mx in MXs, y in Ys, lam in Ls}:
            EB_real_X[mx,y,lam] = 2*sum {x in Xs: (x,y) in Pupil} TelApProp[x,y]*A[x,y]*cos(2*pi*x*mx/lam)*dx;
        subject to st_EB_real {(mx, my) in FPMtrans, lam in Ls}:
            EB_real[mx,my,lam] = 2/lam*sum {y in Ys} EB_real_X[mx,y,lam]*cos(2*pi*y*my/lam)*dy;
        """

        if self.solver['planeofconstr'] == 'Lyot' or self.solver['planeofconstr'] == 'FP2':
            if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon'] is not None: 
                field_propagation_to_Lyot = """
                #---------------------
                var EC_real_X {u in Us, my in MYs, lam in Ls};
                var EC_real {u in Us, v in Vs, lam in Ls};
                
                subject to st_EC_real_X {u in Us, my in MYs, lam in Ls}:
                    EC_real_X[u,my,lam] = 2*sum {mx in MXs: (mx,my) in FPMtrans} FPM[mx,my]*EB_real[mx,my,lam]*cos(2*pi*u*mx/lam)*dmx;
                subject to st_EC_real {(u,v) in Lyot union LyotDarkZone, lam in Ls}:
                    EC_real[u,v,lam] = 2/lam*sum {my in MYs} EC_real_X[u,my,lam]*cos(2*pi*v*my/lam)*dmy;

                """
            else:
                field_propagation_to_Lyot = """
                #---------------------
                var EC_real_X {u in Us, my in MYs, lam in Ls};
                var EC_real {u in Us, v in Vs, lam in Ls};
                
                subject to st_EC_real_X {u in Us, my in MYs, lam in Ls}:
                    EC_real_X[u,my,lam] = 2*sum {mx in MXs: (mx,my) in FPMtrans} FPM[mx,my]*EB_real[mx,my,lam]*cos(2*pi*u*mx/lam)*dmx;
                subject to st_EC_real {(u,v) in Lyot, lam in Ls}:
                    EC_real[u,v,lam] = 2/lam*sum {my in MYs} EC_real_X[u,my,lam]*cos(2*pi*v*my/lam)*dmy;
                
                """
        else:
            field_propagation_to_Lyot = """
	    """

        if self.solver['planeofconstr'] == 'FP2':
	    field_propagation_to_FP2 = """
            #---------------------
            var ED_real_X {xi in Xis, v in Vs, lam in Ls};
            var ED_real {xi in Xis, eta in Etas, lam in Ls};
            
            subject to st_ED_real_X {xi in Xis, v in Vs, lam in Ls}: 
                ED_real_X[xi,v,lam] = 2*sum {u in Us: (u,v) in Lyot} LS[u,v]*EC_real[u,v,lam]*cos(2*pi*u*xi/lam)*du;
            subject to st_ED_real {(xi, eta) in DarkHole, lam in Ls}: 
                ED_real[xi,eta,lam] = 2/lam*sum {v in Vs} ED_real_X[xi,v,lam]*cos(2*pi*v*eta/lam)*dv;
            
	    """
        else:
	    field_propagation_to_FP2 = """
	    """

        if self.solver['planeofconstr'] == 'Lyot' or self.solver['planeofconstr'] == 'FP2':
	    field_propagation_unocc_to_FP1 = """
            #---------------------
            var EB00_real_X {mx in MXs, y in Ys};
            var EB00_real {mx in MXs, my in MYs};
            var EC00_real_X {u in Us, my in MYs};
            var EC00_real {u in Us, v in Vs};
            var E00_ref := 0.0;

            subject to st_EB00_real_X {mx in MXs, y in Ys}:
                EB00_real_X[mx,y] = 2*sum {x in Xs: (x,y) in Pupil} TelApProp[x,y]*A[x,y]*cos(2*pi*x*mx)*dx;
            subject to st_EB00_real {(mx, my) in FPMall}: 
                EB00_real[mx,my] = 2*sum {y in Ys} EB00_real_X[mx,y]*cos(2*pi*y*my)*dy;
            """
	    field_propagation_unocc_to_Lyot = """
            subject to st_EC00_real_X {u in Us, my in MYs}:
                EC00_real_X[u,my] = 2*sum {mx in MXs: (mx,my) in FPMall} EB00_real[mx,my]*cos(2*pi*u*mx)*dmx;
            subject to st_EC00_real {(u,v) in Lyot}:
                EC00_real[u,v] = 2*sum {my in MYs} EC00_real_X[u,my]*cos(2*pi*v*my)*dmy;
            """
        else: # stop at FP1
            field_propagation_unocc_to_FP1 = """
            var E00_ref := 0.0;
            subject to st_E00_ref:
                E00_ref = 4.*sum {(x,y) in Pupil} A[x,y]*dx*dy;
            """
            field_propagation_unocc_to_Lyot = """
            """
        if self.solver['planeofconstr'] == 'FP2':
            field_propagation_unocc_to_FP2 = """
            subject to st_E00_ref:
                E00_ref = 4.*sum {u in Us, v in Vs: (u,v) in Lyot} LS[u,v]*EC00_real[u,v]*du*dv;
            """
        elif self.solver['planeofconstr'] == 'Lyot':
            field_propagation_unocc_to_FP2 = """
            subject to st_E00_ref:
                E00_ref = 1;
            """
        else:
            field_propagation_unocc_to_FP2 = """
            """

        if self.solver['planeofconstr'] == 'FP1':
            constraints = """
            #---------------------
            maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;
           
            subject to sidelobe_zero_real_pos {(mx,my) in FPMtrans}:
                EB_real[mx,my,1] <= 10^(-c/2)*E00_ref/sqrt(2.);
            subject to sidelobe_zero_real_neg {(mx,my) in FPMtrans}:
                EB_real[mx,my,1] >= -10^(-c/2)*E00_ref/sqrt(2.);
            """
        else: # for now Lyot constraint case is not defined
            if self.design['LS']['aligntol'] is not None and self.design['LS']['aligntolcon'] is not None:
                constraints = """
                #---------------------
                maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;
               
                subject to Lyot_aligntol_constr_pos {(x,y) in LyotDarkZone, lam in Ls}:
                    EC_real[x,y,lam] <= 10^-s;
                subject to Lyot_aligntol_constr_neg {(x,y) in LyotDarkZone, lam in Ls}:
                    EC_real[x,y,lam] >= -10^-s;
                subject to sidelobe_zero_real_pos {(xi,eta) in DarkHole, lam in Ls}:
                    ED_real[xi,eta,lam] <= 10^(-c/2)*E00_ref/lam/sqrt(2.);
                subject to sidelobe_zero_real_neg {(xi,eta) in DarkHole, lam in Ls}:
                    ED_real[xi,eta,lam] >= -10^(-c/2)*E00_ref/lam/sqrt(2.);
                """
            else:
                constraints = """
                #---------------------
                maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;
                
                subject to sidelobe_zero_real_pos {(xi,eta) in DarkHole, lam in Ls}:
                    ED_real[xi,eta,lam] <= 10^(-c/2)*E00_ref/lam/sqrt(2.);
                subject to sidelobe_zero_real_neg {(xi,eta) in DarkHole, lam in Ls}:
                    ED_real[xi,eta,lam] >= -10^(-c/2)*E00_ref/lam/sqrt(2.);
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

        mod_fobj.write( textwrap.dedent(field_propagation_to_FP1) )
        mod_fobj.write( textwrap.dedent(field_propagation_to_Lyot) )
        mod_fobj.write( textwrap.dedent(field_propagation_to_FP2) )

        mod_fobj.write( textwrap.dedent(field_propagation_unocc_to_FP1) )
        mod_fobj.write( textwrap.dedent(field_propagation_unocc_to_Lyot) )
        mod_fobj.write( textwrap.dedent(field_propagation_unocc_to_FP2) )

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
            #SBATCH --qos=allnccs
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

class HalfplaneSPLC(SPLC): # Zimmerman SPLC subclass for the half-plane symmetry case
    def __init__(self, **kwargs):
        super(HalfplaneSPLC, self).__init__(**kwargs)
        self.amplname_coron = "SPLC_half"
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
            self.fileorg['FPM fname'] = \
              os.path.join( self.fileorg['FPM dir'], "FPM_quart_diaphragm_{0:03d}R{1:03d}M{2:03d}_{3:s}{4:03d}deg.dat".format(
                                                      int(10*self.design['FPM']['fpmres']*self.design['FPM']['R0']),
                                                      int(10*self.design['FPM']['fpmres']*self.design['FPM']['R1']),
                                                      int(np.ceil(self.design['FPM']['fpmres']*self.design['FPM']['R1'])),
                                                      self.design['FPM']['orient'], self.design['FPM']['openang']) )

        if 'LS fname' not in self.fileorg or self.fileorg['LS fname'] is None:
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}_N{8:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['prim'], self.design['Pupil']['secobs'],
                                                         self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'],
                                                         self.design['LS']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}Pad{5:02d}_N{6:04d}.dat".format(
                                                         self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                         self.design['Pupil']['secobs'], self.design['Pupil']['thick'],
                                                         self.design['LS']['pad'], self.design['LS']['N'])) )
            else:
                self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                         "{0:s}{1:02d}D{2:02d}_clear_N{3:04d}.dat".format(self.design['LS']['shape'],
                                                         self.design['LS']['id'], self.design['LS']['od'], self.design['LS']['N'])) )

        if self.design['LS']['aligntol'] is not None and ('LDZ fname' not in self.fileorg or self.fileorg['LDZ fname'] is None):
            if self.design['LS']['obscure'] == 2:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_half_" + \
                                                          "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}_Tol{8:02d}_N{9:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['Pupil']['prim'], self.design['Pupil']['secobs'],
                                                          self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), self.design['LS']['pad'],
                                                          self.design['LS']['aligntol'], self.design['LS']['N'])) )
            elif self.design['LS']['obscure'] == 1:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_half_" + \
                                                          "{0:s}{1:02d}D{2:02d}_{3:s}{4:s}Pad{5:02d}_Tol{6:02d}_N{7:04d}.dat".format(
                                                          self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                          self.design['Pupil']['secobs'], self.design['Pupil']['thick'],
                                                          self.design['LS']['pad'], self.design['LS']['aligntol'], self.design['LS']['N'])) )
            else:
                self.fileorg['LDZ fname'] = os.path.join( self.fileorg['LS dir'], ("LDZ_half_" + \
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
        # AMPL program to optimize a half-plane symmetric SPLC
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
        set Ys := setof {j in -N_A+0.5..N_A-0.5 by 1} j*dy;
        
        set MXs := setof {i in 0.5..M-0.5 by 1} i*dmx;
        set MYs := setof {j in 0.5..M-0.5 by 1} j*dmy;

        set Us := setof {i in 0.5..N_L-0.5 by 1} i*du;
        set Vs := setof {j in -N_L+0.5..N_L-0.5 by 1} j*dv;

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

        if self.design['Pupil']['edge'] == 'floor': # floor to binary
            define_pupil_and_telap = """
            #---------------------

            set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] == 1} (x,y);
            param TelApProp {x in Xs, y in Ys};
            let {x in Xs, y in Ys} TelApProp[x,y] := 0;
            let {(x,y) in Pupil} TelApProp[x,y] := 1;
            """
        elif self.design['Pupil']['edge'] == 'round': # round to binary
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
            field_propagation_to_FP1 = """
            #---------------------
            """
        else:
            field_propagation_to_FP1 = """
            #---------------------
            var EB_part {mx in MXs, y in Ys, lam in Ls};
            var EB_real {mx in MXs, my in MYs, lam in Ls};
            var EB_imag {mx in MXs, my in MYs, lam in Ls};
            
            subject to st_EB_part {mx in MXs, y in Ys, lam in Ls}:
                EB_part[mx,y,lam] = 2*sum {x in Xs: (x,y) in Pupil} TelApProp[x,y]*A[x,y]*cos(2*pi*x*mx/lam)*dx;
            subject to st_EB_real {(mx, my) in FPMtrans, lam in Ls}:
                EB_real[mx,my,lam] = 1/lam*sum {y in Ys} EB_part[mx,y,lam]*cos(2*pi*y*my/lam)*dy;
            subject to st_EB_imag {(mx, my) in FPMtrans, lam in Ls}:
                EB_imag[mx,my,lam] = -1/lam*sum {y in Ys} EB_part[mx,y,lam]*sin(2*pi*y*my/lam)*dy;
            """

        if self.solver['planeofconstr'] == 'Lyot' or self.solver['planeofconstr'] == 'FP2':
            field_propagation_to_Lyot = """
            
            #---------------------
            var EC_part_real {u in Us, my in MYs, lam in Ls};
            var EC_part_imag {u in Us, my in MYs, lam in Ls};
            var EC_real {u in Us, v in Vs, lam in Ls};
            
            subject to st_EC_part_real {u in Us, my in MYs, lam in Ls}:
                EC_part_real[u,my,lam] = 2*sum {mx in MXs: (mx,my) in FPMtrans} FPM[mx,my]*EB_real[mx,my,lam]*cos(2*pi*u*mx/lam)*dmx;
            subject to st_EC_part_imag {u in Us, my in MYs, lam in Ls}:
                EC_part_imag[u,my,lam] = 2*sum {mx in MXs: (mx,my) in FPMtrans} FPM[mx,my]*EB_imag[mx,my,lam]*cos(2*pi*u*mx/lam)*dmx;
            subject to st_EC_real {(u,v) in Lyot, lam in Ls}:
                EC_real[u,v,lam] = 2/lam*sum {my in MYs} ( EC_part_real[u,my,lam]*cos(2*pi*v*my/lam) + EC_part_imag[u,my,lam]*sin(2*pi*v*my/lam) )*dmy;
            """
        else:
            field_propagation_to_Lyot = """
	    """
            
        if self.solver['planeofconstr'] == 'FP2':
	    field_propagation_to_FP2 = """
            #---------------------
            var ED_part {xi in Xis, v in Vs, lam in Ls};
            var ED_real {xi in Xis, eta in Etas, lam in Ls};
            var ED_imag {xi in Xis, eta in Etas, lam in Ls};
            
            subject to st_ED_part {xi in Xis, v in Vs, lam in Ls}: 
                ED_part[xi,v,lam] = 2*sum {u in Us: (u,v) in Lyot} LS[u,v]*EC_real[u,v,lam]*cos(2*pi*u*xi/lam)*du;
            subject to st_ED_real {(xi, eta) in DarkHole, lam in Ls}: 
                ED_real[xi,eta,lam] = 1/lam*sum {v in Vs} ED_part[xi,v,lam]*cos(2*pi*v*eta/lam)*dv;
            subject to st_ED_imag {(xi, eta) in DarkHole, lam in Ls}: 
                ED_imag[xi,eta,lam] = -1/lam*sum {v in Vs} ED_part[xi,v,lam]*sin(2*pi*v*eta/lam)*dv;
            
	    """
        else:
	    field_propagation_to_FP2 = """
	    """

        if self.solver['planeofconstr'] == 'Lyot' or self.solver['planeofconstr'] == 'FP2':
	    field_propagation_unocc_to_FP1 = """
            #---------------------
            var EB00_part {mx in MXs, y in Ys};
            var EB00_real {mx in MXs, my in MYs};
            var EB00_imag {mx in MXs, my in MYs};
            var EC00_part_real {u in Us, my in MYs};
            var EC00_part_imag {u in Us, my in MYs};
            var EC00_real {u in Us, v in Vs};
            var ED00_real := 0.0;

            subject to st_EB00_part {mx in MXs, y in Ys}:
                EB00_part[mx,y] = 2*sum {x in Xs: (x,y) in Pupil} TelApProp[x,y]*A[x,y]*cos(2*pi*x*mx)*dx;
            subject to st_EB00_real {(mx, my) in FPMall}: 
                EB00_real[mx,my] = sum {y in Ys} EB00_part[mx,y]*cos(2*pi*y*my)*dy;
            subject to st_EB00_imag {(mx, my) in FPMall}: 
                EB00_imag[mx,my] = sum {y in Ys} -EB00_part[mx,y]*sin(2*pi*y*my)*dy;
            """
	    field_propagation_unocc_to_Lyot = """
            subject to st_EC00_part_real {u in Us, my in MYs}:
                EC00_part_real[u,my] = 2*sum {mx in MXs: (mx,my) in FPMall} EB00_real[mx,my]*cos(2*pi*u*mx)*dmx;
            subject to st_EC00_part_imag {u in Us, my in MYs}:
                EC00_part_imag[u,my] = 2*sum {mx in MXs: (mx,my) in FPMall} EB00_imag[mx,my]*cos(2*pi*u*mx)*dmx;
            subject to st_EC00_real {(u,v) in Lyot}:
                EC00_real[u,v] = 2*sum {my in MYs} ( EC00_part_real[u,my]*cos(2*pi*v*my) + EC00_part_imag[u,my]*sin(2*pi*v*my) )*dmy;
            """
        else: # stop at FP1
            field_propagation_unocc_to_FP1 = """
            var E00_ref := 0.0;
            subject to st_E00_ref:
                E00_ref = 2.*sum {(x,y) in Pupil} A[x,y]*dx*dy;
            """
            field_propagation_unocc_to_Lyot = """
            """
        if self.solver['planeofconstr'] == 'FP2':
            field_propagation_unocc_to_FP2 = """
            subject to st_ED00_real:
                ED00_real = 2.*sum {u in Us, v in Vs: (u,v) in Lyot} LS[u,v]*EC00_real[u,v]*du*dv;
            """
        elif self.solver['planeofconstr'] == 'Lyot':
            field_propagation_unocc_to_FP2 = """
            subject to st_E00_ref:
                E00_ref = 1;
            """
        else:
            field_propagation_unocc_to_FP2 = """
            """

        if self.solver['planeofconstr'] == 'FP1':
            constraints = """
            #---------------------
            maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;
            
            subject to sidelobe_real_pos {(mx,my) in FPMtrans}:
                EB_real[mx,my,1] <= 10^(-c/2)*E00_ref/sqrt(2.);
            subject to sidelobe_real_neg {(mx,my) in FPMtrans}:
                EB_real[mx,my,1] >= -10^(-c/2)*E00_ref/sqrt(2.);
            subject to sidelobe_imag_pos {(mx,my) in FPMtrans}:
                EB_imag[mx,my,1] <= 10^(-c/2)*E00_ref/sqrt(2.);
            subject to sidelobe_imag_neg {(mx,my) in FPMtrans}:
                EB_imag[mx,my,1] >= -10^(-c/2)*E00_ref/sqrt(2.);
            """
        else: 
            constraints = """
            #---------------------
            maximize throughput: sum{(x,y) in Pupil} A[x,y]*dx*dy/TR;

            subject to sidelobe_real_pos {(xi,eta) in DarkHole, lam in Ls}: 
                ED_real[xi,eta,lam] <= 10^(-c/2)*ED00_real/lam/sqrt(2.);
            subject to sidelobe_real_neg {(xi,eta) in DarkHole, lam in Ls}: 
                ED_real[xi,eta,lam] >= -10^(-c/2)*ED00_real/lam/sqrt(2.);
            subject to sidelobe_imag_pos {(xi,eta) in DarkHole, lam in Ls}: 
                ED_imag[xi,eta,lam] <= 10^(-c/2)*ED00_real/lam/sqrt(2.);
            subject to sidelobe_imag_neg {(xi,eta) in DarkHole, lam in Ls}:
                ED_imag[xi,eta,lam] >= -10^(-c/2)*ED00_real/lam/sqrt(2.);
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
        mod_fobj.write( textwrap.dedent(sets_and_arrays_part1) )
        mod_fobj.write( textwrap.dedent(sets_and_arrays_part2) )

        mod_fobj.write( textwrap.dedent(field_propagation_to_FP1) )
        mod_fobj.write( textwrap.dedent(field_propagation_to_Lyot) )
        mod_fobj.write( textwrap.dedent(field_propagation_to_FP2) )

        mod_fobj.write( textwrap.dedent(field_propagation_unocc_to_FP1) )
        mod_fobj.write( textwrap.dedent(field_propagation_unocc_to_Lyot) )
        mod_fobj.write( textwrap.dedent(field_propagation_unocc_to_FP2) )

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
            #SBATCH --qos=allnccs
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
                                                           ('thick',(str, '025')), ('centobs',(int, 1)),
                                                           ('gap',(int, 1)), ('edge',(str, 'gray'))]) ),
                                   ( 'FPM', OrderedDict([('rad',(float, 4.)), ('M',(int, 60))]) ),
                                   ( 'LS', OrderedDict([('shape',(str, 'ann')), ('id',(int, 20)), ('od',(int, None)), ('obscure',(int, 0)),
                                                        ('pad',(int, 0)), ('aligntol',(int, None)), ('aligntolcon',(float, 3.))]) ),
                                   ( 'Image', OrderedDict([('c',(float, 10.)), ('ida',(float, -0.5)), ('oda',(float, 10.)), ('bowang',(int, 180)),
                                                           ('bw',(float, 0.10)), ('Nlam',(int, 1)), ('fpres',(int,2)),
                                                           ('wingang',(float, None)), ('incon',(float, None)), ('wingcon',(float, None))]) ) ])
    _eval_fields =   { 'Pupil': _design_fields['Pupil'], 'FPM': _design_fields['FPM'], \
                       'LS': _design_fields['LS'], 'Image': _design_fields['Image'], \
                       'Tel': {'TelAp diam':(float, 12.)}, 'Target': {}, 'Aber': {}, 'WFSC': {} }
    _LS_OD_map = {'hex1':76, 'hex2':82, 'hex3':81, 'hex4':82, 'pie08':90, 'pie12':90, 'key24':90, 'circ':90, 'wfirst':90, 'wfirstCycle5':90}
    _prim_secobs_map = {'hex1':'X', 'hex2':'X', 'hex3':'X', 'hex4':'X', 'ochex1':'X', 'ochex2':'X', 'ochex3':'X', 'ochex4':'X',
                        'pie08':'Cross', 'pie12':'Cross', 'key24':'Cross', 'circ':'Cross', 'wfirst':'WFIRST', 'wfirstCycle5':'wfirstCycle5'}

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
            self.design['Image']['Nlam'] = int(np.ceil(self.design['Image']['bw']/(0.12/4)))
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
        if self.design['LS']['shape'] is 'perim':
            self.design['LS']['shape'] = self.design['Pupil']['prim'] + 'P'
        if verbose: # Print summary of the set parameters
            logging.info("Design parameters: {}".format(self.design))
            logging.info("Optimization and solver parameters: {}".format(self.solver))
            logging.info("File organization parameters: {}".format(self.fileorg))
     
        self.amplname_coron = "APLC_full"
        if self.design['Pupil']['prim'] is 'wfirst':
            self.telap_descrip = "wfirstCobs{0:02d}sthick{1:s}_N{2:04d}".format(self.design['Pupil']['centobs'], self.design['Pupil']['thick'], \
                                                                                self.design['Pupil']['N'])
        elif self.design['Pupil']['prim'] is 'wfirstCycle5':
            self.telap_descrip = "wfirstCycle5_N{0:04d}".format(self.design['Pupil']['N'])
        else:
            self.telap_descrip = "{0:s}{1:s}{2:s}cobs{3:d}gap{4:d}_N{5:04d}".format(self.design['Pupil']['prim'], self.design['Pupil']['secobs'], \
                                                                                    self.design['Pupil']['thick'], int(self.design['Pupil']['centobs']), \
                                                                                    self.design['Pupil']['gap'], self.design['Pupil']['N'])
        self.amplname_pupil = "{0:s}{1:s}".format(self.telap_descrip, self.design['Pupil']['edge'][0])

        self.amplname_fpm = "FPM{:02}M{:03}".format(int(round(100*self.design['FPM']['rad'])), self.design['FPM']['M'])
        if self.design['LS']['obscure'] == 2: # LS includes primary and secondary aperture features
            self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}{3:s}{4:s}{5:s}cobs{6:d}Pad{7:02d}".format(self.design['LS']['shape'], self.design['LS']['id'], \
                               self.design['LS']['od'], self.design['Pupil']['prim'], self.design['Pupil']['secobs'], self.design['Pupil']['thick'], \
                               int(self.design['Pupil']['centobs']), self.design['LS']['pad'])
        elif self.design['LS']['obscure'] == 1: # LS includes secondary aperture features
            if self.design['Pupil']['prim'] is 'wfirst':
                self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}WFIRST{3:s}Pad{4:02d}".format(self.design['LS']['shape'],
                                    self.design['LS']['id'], self.design['LS']['od'],
                                    self.design['Pupil']['thick'], self.design['LS']['pad'])
            elif self.design['Pupil']['prim'] == 'wfirstCycle5':
                self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}Pad{3:02d}".format(self.design['LS']['shape'],
                                    self.design['LS']['id'], self.design['LS']['od'], self.design['LS']['pad'])
            else:
                self.amplname_ls = "LS{0:s}{1:02d}D{2:02d}{3:s}{4:s}Pad{5:02d}".format(self.design['LS']['shape'],
                                    self.design['LS']['id'], self.design['LS']['od'],
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
        elif self.design['Image']['bowang'] != 180: # bowtie image constraints
		    if self.design['Image']['bowang'] >= 0: # horizontal bowtie
				self.amplname_image = "Img{:03}C_{:02}DA{:03}HB{:03}_BW{:02}Nlam{:02}fpres{:1}".format(int(round(10*self.design['Image']['c'])),
                                       int(round(10*(self.design['FPM']['rad']+self.design['Image']['ida']))), int(round(10*self.design['Image']['oda'])),
                                       np.abs(self.design['Image']['bowang']), int(round(100*self.design['Image']['bw'])),
									   self.design['Image']['Nlam'], self.design['Image']['fpres'])
		    else:
				self.amplname_image = "Img{:03}C_{:02}DA{:03}VB{:03}_BW{:02}Nlam{:02}fpres{:1}".format(int(round(10*self.design['Image']['c'])),
                                       int(round(10*(self.design['FPM']['rad']+self.design['Image']['ida']))), int(round(10*self.design['Image']['oda'])),
                                       np.abs(self.design['Image']['bowang']), int(round(100*self.design['Image']['bw'])),
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

    def get_coron_masks(self, use_gray_gap_zero=True):
        TelAp_basename = os.path.basename(self.fileorg['TelAp fname'])
        if use_gray_gap_zero:
            gapstr_beg = TelAp_basename.find('gap')
            TelAp_nopad_basename = TelAp_basename.replace(TelAp_basename[gapstr_beg:gapstr_beg+4], 'gap0')
            TelAp_nopad_fname = os.path.join( os.path.dirname(self.fileorg['TelAp fname']), TelAp_nopad_basename )
            TelAp_p = np.loadtxt(TelAp_nopad_fname)
        elif self.design['Pupil']['edge'] == 'floor': # floor to binary
            TelAp_p = np.floor(np.loadtxt(self.fileorg['TelAp fname']))
        else:
            TelAp_p = np.round(np.loadtxt(self.fileorg['TelAp fname']))

        A_col = np.loadtxt(self.fileorg['sol fname'])[:,-1]
        FPM_p = np.loadtxt(self.fileorg['FPM fname'])
        LS_p = np.loadtxt(self.fileorg['LS fname'])
        A_p = A_col.reshape(TelAp_p.shape)
        if isinstance(self, QuarterplaneAPLC):
            TelAp = np.concatenate((np.concatenate((TelAp_p[::-1,::-1], TelAp_p[:,::-1]),axis=0),
                                    np.concatenate((TelAp_p[::-1,:], TelAp_p),axis=0)), axis=1)
            Apod = np.concatenate((np.concatenate((A_p[::-1,::-1], A_p[:,::-1]),axis=0),
                                   np.concatenate((A_p[::-1,:], A_p),axis=0)), axis=1)
            FPM = np.concatenate((np.concatenate((FPM_p[::-1,::-1], FPM_p[:,::-1]),axis=0),
                                  np.concatenate((FPM_p[::-1,:], FPM_p),axis=0)), axis=1)
            LS = np.concatenate((np.concatenate((LS_p[::-1,::-1], LS_p[:,::-1]),axis=0),
                                 np.concatenate((LS_p[::-1,:], LS_p),axis=0)), axis=1)
        elif isinstance(self, HalfplaneAPLC):
            TelAp = np.concatenate((TelAp_p[:,::-1], TelAp_p), axis=1)
            Apod = np.concatenate((A_p[:,::-1], A_p), axis=1)
            FPM = np.concatenate((np.concatenate((FPM_p[::-1,::-1], FPM_p[:,::-1]),axis=0),
                                  np.concatenate((FPM_p[::-1,:], FPM_p),axis=0)), axis=1)
            LS = np.concatenate((LS_p[:,::-1], LS_p), axis=1)
        else:
            TelAp = TelAp_p
            Apod = A_p
            FPM = FPM_p
            LS = LS_p
        return TelAp, Apod, FPM, LS

    def get_coords(self, fp2res=8, rho_out=None, Nlam=None): # for APLC
        D = 1.
        N = self.design['Pupil']['N']
        bw = self.design['Image']['bw']
        M_fp1 = self.design['FPM']['M']
        if rho_out is None:
            rho_out = self.design['Image']['oda'] + 0.5
        if Nlam is None:
            Nlam = self.design['Image']['Nlam']
        M_fp2 = int(np.ceil(rho_out*fp2res))
        
        # pupil plane
        dx = (D/2)/N
        dy = dx
        xs = np.matrix(np.linspace(-N+0.5,N-0.5,2*N)*dx)
        XX, YY = np.meshgrid(np.array(xs), np.array(xs))
        
        # FPM
        dmx = self.design['FPM']['rad']/M_fp1
        dmy = dmx
        mxs = np.matrix(np.linspace(-M_fp1+0.5,M_fp1-0.5,2*M_fp1)*dmx)
        mys = mxs

        # FP2
        dxi = 1./fp2res
        xis = np.matrix(np.linspace(-M_fp2+0.5,M_fp2-0.5,2*M_fp2)*dxi)

        # wavelength ratios
        wrs = np.linspace(1.-bw/2, 1.+bw/2, Nlam)

        return xs, dx, XX, YY, mxs, dmx, xis, dxi, wrs

    def get_onax_psf(self, fp2res=8, rho_inc=0.25, rho_out=None, Nlam=None): # for APLC class
        if self.design['Pupil']['edge'] == 'floor': # floor to binary
            TelAp_p = np.floor(np.loadtxt(self.fileorg['TelAp fname'])).astype(int)
        elif self.design['Pupil']['edge'] == 'round': # round to binary
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
            Psi_A = TelAp*A
            Psi_B = dx*dx/wr*np.exp(-1j*2*np.pi/wr*mxs.T*xs)*Psi_A*np.exp(-1j*2*np.pi/wr*xs.T*mxs)
            Psi_B_stop = np.multiply(Psi_B, FPM)
            Psi_C = Psi_A[::-1,::-1] - dmx*dmx/wr*np.exp(-1j*2*np.pi/wr*xs.T*mxs)*Psi_B_stop*np.exp(-1j*2*np.pi/wr*mxs.T*xs)
            Psi_C_stop = np.multiply(Psi_C, LS)
            Psi_D = dx*dx/wr*np.exp(-1j*2*np.pi/wr*xis.T*xs)*Psi_C_stop*np.exp(-1j*2*np.pi/wr*xs.T*xis)
            Psi_D_0_peak = np.sum(A*TelAp*LS)*dx*dx/wr
            intens_polychrom[wi,:,:] = np.power(np.absolute(Psi_D)/Psi_D_0_peak, 2)
             
        seps = np.arange(self.design['FPM']['rad']+self.design['Image']['ida'], rho_out, rho_inc)
        radial_intens_polychrom = np.zeros((len(wrs), len(seps)))
        XXs = np.asarray(np.dot(np.matrix(np.ones(xis.shape)).T, xis))
        YYs = np.asarray(np.dot(etas.T, np.matrix(np.ones(etas.shape))))
        RRs = np.sqrt(XXs**2 + YYs**2)

        if 'bowang' in self.design['Image'] and self.design['Image']['bowang'] != 180: # Define bowtie angle constraints
            if self.design['Image']['bowang'] >= 0: # horizontal dark zone
                theta_quad = np.rad2deg(np.arctan2(YYs[M_fp2:,M_fp2:], XXs[M_fp2:,M_fp2:]))
                theta_quad_mask = np.less(theta_quad, self.design['Image']['bowang']/2)
            else: # vertical dark zone
                theta_quad = np.rad2deg(np.arctan2(YYs[M_fp2:,M_fp2:], XXs[M_fp2:,M_fp2:]))
                theta_quad_mask = np.greater(theta_quad, self.design['Image']['bowang']/2)
            theta_rhs_mask = np.concatenate((theta_quad_mask[::-1,:], theta_quad_mask), axis=0)
            theta_mask = np.concatenate((theta_rhs_mask[:,::-1], theta_rhs_mask), axis=1)

        for si, sep in enumerate(seps):
            r_in = np.max([seps[0], sep-0.25])
            r_out = np.min([seps[-1], sep+0.25])
            if 'bowang' in self.design['Image'] and self.design['Image']['bowang'] != 180: # apply bowtie angle constraints
                meas_mask = (theta_mask & (RRs >= r_in) & (RRs <= r_out))
                meas_ann_ind = np.nonzero(np.ravel(meas_mask))[0]
            else: # no angle constraints
                meas_ann_ind = np.nonzero(np.logical_and(np.greater_equal(RRs, r_in).ravel(),
                                                         np.less_equal(RRs, r_out).ravel()))[0]
            for wi, wr in enumerate(wrs):
                radial_intens_polychrom[wi, si] = np.mean(np.ravel(intens_polychrom[wi,:,:])[meas_ann_ind])

        return xis, intens_polychrom, seps, radial_intens_polychrom

    def get_metrics(self, fp2res=16, rho_out=None, Nlam=None, verbose=True): # for APLC class
        TelAp_basename = os.path.basename(self.fileorg['TelAp fname'])
        gapstr_beg = TelAp_basename.find('gap')
        TelAp_nopad_basename = TelAp_basename.replace(TelAp_basename[gapstr_beg:gapstr_beg+4], 'gap0')
        TelAp_nopad_fname = os.path.join( os.path.dirname(self.fileorg['TelAp fname']), TelAp_nopad_basename )
        #if self.design['Pupil']['edge'] == 'floor': # floor to binary
        #    TelAp_p = np.floor(np.loadtxt(self.fileorg['TelAp fname'])).astype(int)
        #elif self.design['Pupil']['edge'] == 'round': # round to binary
        #    TelAp_p = np.round(np.loadtxt(self.fileorg['TelAp fname'])).astype(int)
        #else:
        #    TelAp_p = np.loadtxt(self.fileorg['TelAp fname'])
        if os.path.exists(TelAp_nopad_fname):
            TelAp_p = np.loadtxt(TelAp_nopad_fname)
            telap_flag = 0
        else:
            TelAp_p = np.loadtxt(self.fileorg['TelAp fname'])
            telap_flag = 1
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
            Psi_D_0 = dx*dy/wr*np.dot(np.dot(np.exp(-1j*2*np.pi/wr*np.dot(xis.T, xs)), TelAp*A*LS[::-1,::-1]),
                                             np.exp(-1j*2*np.pi/wr*np.dot(xs.T, xis)))
            intens_D_0_polychrom[wi] = np.power(np.absolute(Psi_D_0), 2)
            intens_D_0_peak_polychrom[wi] = (np.sum(TelAp*A*LS[::-1,::-1])*dx*dy/wr)**2
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
        self.eval_metrics['fwhm circ thrupt'] = fwhm_sum_APLC/(np.pi/4)
        self.eval_metrics['p7ap thrupt'] = p7ap_sum_APLC/np.sum(np.power(TelAp,2)*dx*dx)
        self.eval_metrics['p7ap circ thrupt'] = p7ap_sum_APLC/(np.pi/4)
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
            print("Band-averaged half-max throughput, circ. ref.: {:.2f}%".format(100*self.eval_metrics['fwhm circ thrupt']))
            print("Band-averaged r=.7 lam/D throughput: {:.2f}%".format(100*self.eval_metrics['p7ap thrupt']))
            print("Band-averaged r=.7 lam/D throughput, circ. ref.: {:.2f}%".format(100*self.eval_metrics['p7ap circ thrupt']))
            print("Band-averaged relative half-max throughput: {:.2f}%".format(100*self.eval_metrics['rel fwhm thrupt']))
            print("Band-averaged relative r=0.7 lam/D throughput: {:.2f}%".format(100*self.eval_metrics['rel p7ap thrupt']))
            print("Band-averaged FWHM PSF area / (lambda0/D)^2: {:.2f}".format(self.eval_metrics['fwhm area']))
        return telap_flag

    def get_yield_input_products(self, pixscale_lamoD=0.25, star_diam_vec=None, Npts_star_diam=7, Nlam=None, norm='aperture'):
        # Assumes quarter-plane symmetry in the final focal plane
        TelAp, Apod, FPM, LS = self.get_coron_masks(use_gray_gap_zero=True)

        if star_diam_vec is None:
            star_diam_vec = np.concatenate([np.linspace(0,0.09,10), np.linspace(0.1, 1, 10), np.array([2., 3., 4.])])
        if Nlam is None:
            Nlam = self.design['Image']['Nlam']
        bw = self.design['Image']['bw']
        wrs = np.linspace(1.-bw/2, 1.+bw/2, Nlam)
        seps = np.arange(self.design['FPM']['rad']+self.design['Image']['ida'],
                         self.design['Image']['oda']+2*pixscale_lamoD, pixscale_lamoD)
        
        N = self.design['Pupil']['N']
        M_fp1 = self.design['FPM']['M']
        fpm_rad = self.design['FPM']['rad']
        rho2 = self.design['Image']['oda'] + 0.5
        M_fp2 = int(np.ceil(rho2/pixscale_lamoD))
        M_fp2_ext = int(np.ceil((rho2+2.5)/pixscale_lamoD))
        wc = M_fp2_ext - M_fp2
        
        # pupil plane
        D = 1.
        dx = (D/2)/N
        dy = dx
        xs = np.matrix(np.linspace(-N+0.5,N-0.5,2*N)*dx)
        XX, YY = np.meshgrid(np.array(xs), np.array(xs))
        
        # FPM
        dmx = fpm_rad/M_fp1
        dmy = dmx
        mxs = np.matrix(np.linspace(-M_fp1+0.5,M_fp1-0.5,2*M_fp1)*dmx)
        
        # FP2
        dxi = pixscale_lamoD
        xis = np.matrix(np.linspace(-M_fp2+0.5,M_fp2-0.5,2*M_fp2)*dxi)
        xis_ext = np.matrix(np.linspace(-M_fp2_ext+0.5,M_fp2_ext-0.5,2*M_fp2_ext)*dxi)
        
        offax_Xis = np.array(xis[0,M_fp2:].T)
        offax_Xis_ext = np.array(xis_ext[0,M_fp2_ext:].T)
        offax_XisEtas = zip(np.ravel(np.ones_like(offax_Xis)*offax_Xis.T), np.ravel(offax_Xis*np.ones_like(offax_Xis.T)))
        offax_XisEtas_ext = zip(np.ravel(np.ones_like(offax_Xis_ext)*offax_Xis_ext.T), np.ravel(offax_Xis_ext*np.ones_like(offax_Xis_ext.T)))

        intens_2d_vs_star_diam = np.zeros((len(star_diam_vec), 2*M_fp2, 2*M_fp2))
        intens_rad_vs_star_diam = np.zeros((len(star_diam_vec), len(seps)))
        offax_psf_map_ext = np.zeros((len(offax_XisEtas_ext), 2*M_fp2, 2*M_fp2))
        offax_psf_map = np.zeros((len(offax_XisEtas), 2*M_fp2, 2*M_fp2))

        for si, star_diam in enumerate(star_diam_vec):
            intens_2d_vs_star_diam[si,:,:], \
            intens_rad_vs_star_diam[si,:] = get_finite_star_aplc_psf(TelAp, Apod, FPM, LS,
                                                                     xs, dx, XX, YY, mxs, dmx, xis, dxi,
                                                                     star_diam, Npts_star_diam, wrs=wrs,
                                                                     seps=seps, norm=norm,
                                                                     get_radial_curve=True)

        if norm is 'aperture':
            contrast_convert_fac = np.sum(np.power(TelAp, 2))*dx*dx/(dxi*dxi) / np.power(np.sum(Apod*LS)*dx*dx, 2)
        else:
            contrast_convert_fac = 1

        for oi, (delta_xi, delta_eta) in enumerate(offax_XisEtas_ext):
            offax_psf_map_ext[oi,:,:] = fast_bandavg_aplc_psf(TelAp, Apod, FPM, LS,
                                                              xs, dx, XX, YY, mxs, dmx, xis, dxi,
                                                              delta_xi, delta_eta, wrs, norm)
            if (delta_xi, delta_eta) in offax_XisEtas:
                ii = offax_XisEtas.index((delta_xi, delta_eta))
                offax_psf_map[ii,:,:] = offax_psf_map_ext[oi,:,:]
            
        sky_trans_map = np.sum(np.concatenate([offax_psf_map_ext,
                                               offax_psf_map_ext[:,::-1,:],
                                               offax_psf_map_ext[:,:,::-1],
                                               offax_psf_map_ext[:,::-1,::-1]], axis=0), axis=0)

        return intens_2d_vs_star_diam, intens_rad_vs_star_diam, np.ravel(xis), seps, star_diam_vec, \
               offax_psf_map, np.array(offax_XisEtas).T, sky_trans_map, contrast_convert_fac

def get_finite_star_aplc_psf(TelAp, Apod, FPM, LS, xs, dx, XX, YY, mxs, dmx, xis, dxi,
                             star_diam_lamoD=0.1, Npts_star_diam=7,
                             wrs=None, seps=None, get_radial_curve=False, norm='peak'):
    if wrs is None:
        wrs = np.linspace(0.95, 1.05, 5)

    disk_vec_lamoD = np.linspace(-star_diam_lamoD/2, star_diam_lamoD/2, Npts_star_diam)

    XiXi, EtaEta = np.meshgrid(disk_vec_lamoD, disk_vec_lamoD)
    star_disk = (XiXi**2 + EtaEta**2 <= (star_diam_lamoD/2)**2)
    disk_samp_XiEta = zip(XiXi[star_disk], EtaEta[star_disk])

    intens_2d_src = np.zeros((xis.shape[1], xis.shape[1]))
    
    for (delxi, deleta) in disk_samp_XiEta:
        intens_2d_bandavg = fast_bandavg_aplc_psf(TelAp, Apod, FPM, LS, xs, dx, XX, YY, mxs, dmx,
                                                  xis, dxi, delxi, deleta, wrs, norm)
        intens_2d_src += intens_2d_bandavg/len(disk_samp_XiEta)
       
    if get_radial_curve: 
        if seps is None:
            seps = np.arange(2.0, 10.25, 0.25)
        intens_radial_src = np.zeros(seps.shape)
        
        XXs = np.asarray(np.dot(np.matrix(np.ones(xis.shape)).T, xis))
        YYs = np.asarray(np.dot(xis.T, np.matrix(np.ones(xis.shape))))
        RRs = np.sqrt(XXs**2 + YYs**2)

        for si, sep in enumerate(seps):
            r_in = np.max([seps[0], sep-0.25])
            r_out = np.min([seps[-1], sep+0.25])
            meas_ann_ind = np.nonzero(np.logical_and(np.greater_equal(RRs, r_in).ravel(),
                                                     np.less_equal(RRs, r_out).ravel()))[0]
            intens_radial_src[si] = np.mean(np.ravel(intens_2d_src)[meas_ann_ind])
            
        return intens_2d_src, intens_radial_src
    else:
        return intens_2d_src
       
def fast_bandavg_aplc_psf(TelAp, A, FPM, LS, xs, dx, XX, YY, mxs, dmx, xis, dxi, delta_xi, delta_eta, wrs,
                          norm = 'peak'):
    # norm parameter is either 'aperture' for integral of illuminated aperture energy (per Stark yield input definition),
    # or 'peak' for unocculted PSF peak (contrast units)
    intens_D_polychrom = np.zeros((wrs.shape[0], xis.shape[1], xis.shape[1]))
    if norm == 'peak':
        intens_norm = np.power(np.sum(A*LS[::-1,::-1])*dx*dx/wrs, 2)
    elif norm == 'aperture':
        intens_norm = np.sum(np.power(TelAp, 2))*dx*dx
    for wi, wr in enumerate(wrs):
        Psi_A = np.exp(-1j*2*np.pi/wr*(delta_xi*XX + delta_eta*YY))
        Psi_A_stop = np.multiply(Psi_A, A)
        Psi_B = dx*dx/wr*np.exp(-1j*2*np.pi/wr*mxs.T*xs)*Psi_A_stop*np.exp(-1j*2*np.pi/wr*xs.T*mxs)
        Psi_B_stop = np.multiply(Psi_B, FPM)
        Psi_C = Psi_A_stop[::-1,::-1] - \
                dmx*dmx/wr*np.exp(-1j*2*np.pi/wr*xs.T*mxs)*Psi_B_stop*np.exp(-1j*2*np.pi/wr*mxs.T*xs)
        Psi_C_stop = np.multiply(Psi_C, LS)
        Psi_D = dx*dx/wr*np.exp(-1j*2*np.pi/wr*xis.T*xs)*Psi_C_stop*np.exp(-1j*2*np.pi/wr*xs.T*xis)
        if norm == 'peak':
            intens_D_polychrom[wi,:,:] = np.power(np.absolute(Psi_D), 2) / intens_norm[wi]
        elif norm == 'aperture':
            intens_D_polychrom[wi,:,:] = dxi*dxi*np.power(np.absolute(Psi_D), 2) / intens_norm
            
    return np.mean(intens_D_polychrom, axis=0)

def fast_bandavg_splc_psf(TelAp, A, FPM, LS, xs, dx, XX, YY, mxs, dmx, us, du, xis, dxi, delta_xi, delta_eta, wrs,
                          norm = 'peak'):
    # norm parameter is either 'aperture' for integral of illuminated aperture energy (per Stark yield input definition),
    # or 'peak' for unocculted PSF peak (contrast units)
    intens_D_polychrom = np.zeros((wrs.shape[0], xis.shape[1], xis.shape[1]))
    if norm == 'aperture':
        intens_norm = np.sum(np.power(TelAp, 2))*dx*dx
    for wi, wr in enumerate(wrs):
        Psi_A = np.exp(-1j*2*np.pi/wr*(delta_xi*XX + delta_eta*YY))
        Psi_A_stop = np.multiply(Psi_A, A)
        Psi_B = dx*dx/wr*np.exp(-1j*2*np.pi/wr*mxs.T*xs)*Psi_A_stop*np.exp(-1j*2*np.pi/wr*xs.T*mxs)
        Psi_B_0 = dx*dx/wr*np.exp(-1j*2*np.pi/wr*mxs.T*xs)*A*np.exp(-1j*2*np.pi/wr*xs.T*mxs)
        Psi_B_stop = np.multiply(Psi_B, FPM)
        Psi_C = dmx*dmx/wr*np.exp(-1j*2*np.pi/wr*us.T*mxs)*Psi_B_stop*np.exp(-1j*2*np.pi/wr*mxs.T*us)
        Psi_C_0 = dmx*dmx/wr*np.exp(-1j*2*np.pi/wr*us.T*mxs)*Psi_B_0*np.exp(-1j*2*np.pi/wr*mxs.T*us)
        Psi_C_stop = np.multiply(Psi_C, LS)
        Psi_C_0_stop = np.multiply(Psi_C_0, LS)
        Psi_D = du*du/wr*np.exp(-1j*2*np.pi/wr*xis.T*us)*Psi_C_stop*np.exp(-1j*2*np.pi/wr*us.T*xis)
        Psi_D_0_peak = du*du/wr*np.sum(Psi_C_0_stop)
        if norm == 'peak':
            intens_D_polychrom[wi,:,:] = np.power(np.absolute(Psi_D), 2) / np.power(np.absolute(Psi_D_0_peak), 2)
        elif norm == 'aperture':
            intens_D_polychrom[wi,:,:] = dxi*dxi*np.power(np.absolute(Psi_D), 2) / intens_norm
        else:
            logging.error('unrecognized value for norm parameter')
            return 1

    return np.mean(intens_D_polychrom, axis=0)
    
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
                if self.design['Pupil']['prim'] is 'wfirst':
                    self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                             "{0:s}{1:02d}D{2:02d}_WFIRST{3:s}Pad{4:02d}_N{5:04d}.dat".format(
                                                             self.design['LS']['shape'], self.design['LS']['id'], self.design['LS']['od'],
                                                             self.design['Pupil']['thick'], self.design['LS']['pad'], self.design['Pupil']['N'])) )
                elif self.design['Pupil']['prim'] is 'wfirstCycle5':
                    self.fileorg['LS fname'] = os.path.join( self.fileorg['LS dir'], ("LS_half_" + \
                                                             "ann{0:02d}D{1:02d}_{2:s}Pad{3:02d}_N{4:04d}.dat".format(
                                                             self.design['LS']['id'], self.design['LS']['od'], self.design['LS']['shape'], 
                                                             self.design['LS']['pad'], self.design['Pupil']['N'])) )
                else:
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
        # load amplgsl.dll;
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
            param ang := {4:0.1f};      # opening angle of dark bowtie constraint region
            
            #---------------------
            param N := {5:d};				# discretization parameter (pupil)
            param M := {6:d};				# discretization parameter (mask)
            param Nimg := {7:d};           # discretization parameter (image)
                                  
            #---------------------
            param bw := {8:0.2f};
            param Nlam := {9:d};
            """.format(self.design['Image']['c'], self.design['FPM']['rad'], self.design['FPM']['rad']+self.design['Image']['ida'],
                       self.design['Image']['oda'], np.abs(self.design['Image']['bowang']),
                       self.design['Pupil']['N'], self.design['FPM']['M'], self.design['Image']['Nimg'],
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

        if self.design['Pupil']['edge'] == 'floor': # floor to binary
            define_pupil_and_telap = """
            #---------------------

            set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] == 1} (x,y);
            param TelApProp {x in Xs, y in Ys};
            let {x in Xs, y in Ys} TelApProp[x,y] := 0;
            let {(x,y) in Pupil} TelApProp[x,y] := 1;
            """
        elif self.design['Pupil']['edge'] == 'round': # round to binary
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
            set Lyot := setof {x in Xs, y in Ys: LS[x,y] >= 0.5} (x,y);
            set LyotDarkZone := setof {x in Xs, y in Ys: LDZ[x,y] == 1 && TelApProp[x,y] > 0} (x,y);

            param TR := sum {(x,y) in Pupil} TelApProp[x,y]*dx*dy; # Transmission of the Pupil. Used for calibration.
            
            var A {x in Xs, y in Ys} >= 0, <= 1, := 0.5;
            """
        else:
            sets_and_arrays = """
            set Mask := setof {mx in MXs, my in MYs: FPM[mx,my] > 0} (mx,my);
            set Lyot := setof {x in Xs, y in Ys: LS[x,y] >= 0.5} (x,y);
            set LyotFlip := setof {x in Xs, y in Ys: LS[x,-y] >= 0.5} (x,y);
            set TransArea := Pupil union LyotFlip;

            param TransAreaNorm := sum {(x,y) in Pupil} TelApProp[x,y];
            
            var A {x in Xs, y in Ys} >= 0, <= 1, := 0.5;
            """

        if self.design['Image']['bowang'] == 180:
            dark_hole = """
            #---------------------

            set DarkHole := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= rho0 && sqrt(xi^2+eta^2) <= rho1} (xi,eta);
            """
        elif self.design['Image']['bowang'] < 0: # vertical bowtie region
            dark_hole = """
            #---------------------

            set DarkHole := setof {xi in Xis, eta in Etas:
                sqrt(xi^2+eta^2) >= rho0 && 
                sqrt(xi^2+eta^2) <= rho1 &&
                eta >= xi*tan(ang/2*pi/180)} (xi,eta);
            """
        else: # horizontal bowtie region
            dark_hole = """
            #---------------------

            set DarkHole := setof {xi in Xis, eta in Etas: 
                sqrt(xi^2+eta^2) >= rho0 && 
                sqrt(xi^2+eta^2) <= rho1 &&
                eta <= xi*tan(ang/2*pi/180)} (xi,eta);
            """

        field_propagation = """
        #---------------------
        var EBm_part {mx in MXs, y in Ys, lam in Ls};
        var EBm_real {mx in MXs, my in MYs, lam in Ls};
        var EBm_imag {mx in MXs, my in MYs, lam in Ls};

        subject to st_EBm_part {mx in MXs, y in Ys, lam in Ls}:
            EBm_part[mx,y,lam] = 2*sum {x in Xs: (x,y) in Pupil} TelApProp[x,y]*A[x,y]*cos(2*pi*x*mx/lam)*dx;
        subject to st_EBm_real {(mx,my) in Mask, lam in Ls}:
            EBm_real[mx,my,lam] = 1/lam*sum {y in Ys} EBm_part[mx,y,lam]*cos(2*pi*y*my/lam)*dy;
        subject to st_EBm_imag {(mx,my) in Mask, lam in Ls}:
            EBm_imag[mx,my,lam] = -1/lam*sum {y in Ys} EBm_part[mx,y,lam]*sin(2*pi*y*my/lam)*dy;
        
        #---------------------
        var EC_part_real {x in Xs, my in MYs, lam in Ls};
        var EC_part_imag {x in Xs, my in MYs, lam in Ls};
        var ECm {x in Xs, y in Ys, lam in Ls};
        
        subject to st_EC_part_real {x in Xs, my in MYs, lam in Ls}:
            EC_part_real[x,my,lam] = 2*sum {mx in MXs: (mx,my) in Mask} FPM[mx,my]*EBm_real[mx,my,lam]*cos(2*pi*x*mx/lam)*dmx;
        subject to st_EC_part_imag {x in Xs, my in MYs, lam in Ls}:
            EC_part_imag[x,my,lam] = 2*sum {mx in MXs: (mx,my) in Mask} FPM[mx,my]*EBm_imag[mx,my,lam]*cos(2*pi*x*mx/lam)*dmx;
        subject to st_EC {(x,y) in Lyot, lam in Ls}:
            ECm[x,y,lam] = 2/lam*sum {my in MYs} ( EC_part_real[x,my,lam]*cos(2*pi*my*y/lam) + EC_part_imag[x,my,lam]*sin(2*pi*my*y/lam) )*dmy;
        
        #---------------------
        var ED_part {xi in Xis, y in Ys, lam in Ls};
        var ED_real {xi in Xis, eta in Etas, lam in Ls};
        var ED_imag {xi in Xis, eta in Etas, lam in Ls};
        
        subject to st_ED_part {xi in Xis, y in Ys, lam in Ls}:
            ED_part[xi,y,lam] = 2*sum {x in Xs: (x,y) in Lyot} (TelApProp[x,-y]*A[x,-y] - ECm[x,y,lam])*cos(2*pi*x*xi/lam)*dx;
        subject to st_ED_real {xi in Xis, eta in Etas, lam in Ls}:
            ED_real[xi,eta,lam] = 1/lam*sum {y in Ys} ED_part[xi,y,lam]*cos(2*pi*y*eta/lam)*dy;
        subject to st_ED_imag {xi in Xis, eta in Etas, lam in Ls}:
            ED_imag[xi,eta,lam] = -1/lam*sum {y in Ys} ED_part[xi,y,lam]*sin(2*pi*y*eta/lam)*dy;
        
        #---------------------
        var ED00_real := 0.0;
        subject to st_ED00_real: ED00_real = 2*sum {x in Xs, y in Ys: (x,y) in LyotFlip} (A[x,y]*TelApProp[x,y])*dx*dy;
        """

        constraints = """
        #---------------------
        maximize throughput: sum{(x,y) in TransArea} A[x,y]/TransAreaNorm;

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
        mod_fobj.write( textwrap.dedent(dark_hole) )
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
            param ang := {5:0.1f};      # opening angle of dark bowtie constraint region
            
            #---------------------
            param N := {6:d};				# discretization parameter (pupil)
            param M := {7:d};				# discretization parameter (mask)
            param Nimg := {8:d};           # discretization parameter (image)
                                  
            #---------------------
            param bw := {9:0.2f};
            param Nlam := {10:d};
            """.format(self.design['Image']['c'], self.design['LS']['s'], self.design['FPM']['rad'],
                       self.design['FPM']['rad']+self.design['Image']['ida'], self.design['Image']['oda'],
					   np.abs(self.design['Image']['bowang']), self.design['Pupil']['N'], self.design['FPM']['M'],
					   self.design['Image']['Nimg'], self.design['Image']['bw'], self.design['Image']['Nlam'])
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
            param bowang := {4:0.1f};
            
            #---------------------
            param N := {5:d};				# discretization parameter (pupil)
            param M := {6:d};				# discretization parameter (mask)
            param Nimg := {7:d};           # discretization parameter (image)
                                  
            #---------------------
            param bw := {8:0.2f};
            param Nlam := {9:d};
            """.format(self.design['Image']['c'], self.design['FPM']['rad'], self.design['FPM']['rad']+self.design['Image']['ida'],
                       self.design['Image']['oda'], np.abs(self.design['Image']['bowang']),
					   self.design['Pupil']['N'], self.design['FPM']['M'], self.design['Image']['Nimg'],
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

        if self.design['Pupil']['edge'] == 'floor': # floor to binary
            define_pupil_and_telap = """
            #---------------------

            set Pupil := setof {x in Xs, y in Ys: TelAp[x,y] == 1} (x,y);
            param TelApProp {x in Xs, y in Ys};
            let {x in Xs, y in Ys} TelApProp[x,y] := 0;
            let {(x,y) in Pupil} TelApProp[x,y] := 1;
            """
        elif self.design['Pupil']['edge'] == 'round': # round to binary
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
            """
        else:
            sets_and_arrays = """
            set Mask := setof {mx in MXs, my in MYs: FPM[mx,my] > 0} (mx,my);
            set Lyot := setof {x in Xs, y in Ys: LS[x,y] > 0} (x,y);

            param TR := sum {(x,y) in Pupil} TelApProp[x,y]*dx*dy; # Transmission of the Pupil. Used for calibration.
            
            var A {x in Xs, y in Ys} >= 0, <= 1, := 0.5;
            """
        if self.design['Image']['bowang'] == 180:
            dark_hole = """
            #---------------------

            set DarkHole := setof {xi in Xis, eta in Etas: sqrt(xi^2+eta^2) >= rho0 && sqrt(xi^2+eta^2) <= rho1} (xi,eta);
            """
        elif self.design['Image']['bowang'] < 0: # vertical bowtie region
            dark_hole = """
            #---------------------

            set DarkHole := setof {xi in Xis, eta in Etas:
                sqrt(xi^2+eta^2) >= rho0 && 
                sqrt(xi^2+eta^2) <= rho1 &&
                eta >= xi*tan(ang/2*pi/180)} (xi,eta);
            """
        else: # horizontal bowtie region
            dark_hole = """
            #---------------------

            set DarkHole := setof {xi in Xis, eta in Etas: 
                sqrt(xi^2+eta^2) >= rho0 && 
                sqrt(xi^2+eta^2) <= rho1 &&
                eta <= xi*tan(ang/2*pi/180)} (xi,eta);
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
        mod_fobj.write( textwrap.dedent(dark_hole) )
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
