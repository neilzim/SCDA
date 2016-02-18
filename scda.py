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

def configure_log(log_fname=None):
    logger = logging.getLogger("scda.logger")
    if not len(logger.handlers):
        logger.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout) # console dump
        ch.setLevel(logging.INFO)
        logger.addHandler(ch)
        if log_fname is not None: # optional file log in addition to the console dump
            fh = logging.FileHandler(log_fname, mode="w")
            fh.setLevel(logging.DEBUG)
            logger.addHandler(fh)

class LyotCoronagraph(object): # Lyot coronagraph base class
    _key_fields = dict([('fileorg', ['ampl src dir', 'TelAp dir', 'FPM dir', 'LS dir', 'sol dir', 'eval dir']),\
                        ('solver', ['constr', 'method', 'presolve'])])
    _solver_menu = dict([('constr',['lin', 'quad']), ('method', ['bar', 'barhom', 'dualsimp']), \
                         ('presolve',[True, False]), ('Nthreads', range(1,33))])
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
                        self.logger.warning("Warning: The specified location for \"{0}\", \"{1}\" does not exist".format(dirkey, location))
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
        if 'Nthreads' not in self.solver or self.solver['Nthreads'] is None: self.solver['Nthreads'] = None

        if 'TelAp_fname' in kwargs:
            self.telap_fname = os.path.join(self.fileorg['TelAp dir'], kwargs['TelAp_fname'])
            if not os.path.exists(self.telap_fname):
                self.logger.warning("Warning: The specified telescope aperture file \"{0}\" does not exist".format(self.telap_fname))
        if 'FPM_fname' in kwargs:
            self.fpm_fname = os.path.join(self.fileorg['FPM dir'], kwargs['FPM_fname'])
            if not os.path.exists(self.fpm_fname):
                self.logger.warning("Warning: The specified focal plane mask file \"{0}\" does not exist".format(self.fpm_fname))
        if 'LS_fname' in kwargs:
            self.ls_fname = os.path.join(self.fileorg['LS dir'], kwargs['LS_fname'])
            if not os.path.exists(self.ls_fname):
                self.logger.warning("Warning: The specified Lyot stop file \"{0}\" does not exist".format(self.ls_fname))

class NdiayeAPLC(LyotCoronagraph): # Image-constrained APLC following N'Diaye et al. (2015, 2016)
    _design_fields = { 'Pupil': {'N':(int, 1000), 'pm':(str, 'hex1'), 'so':(bool, True), 'ss':(str, 'x'), 'sst':(str, '100'), 'tel diam':(float, 12.)}, \
                       'FPM':   {'M':(int, 50), 'rad':(float, 4.)}, \
                       'LS':    {'id':(int, 20), 'od':(int, 90), 'ovsz':(int, 0), 'altol':(int, None)}, \
                       'Image': {'c':(float, 10.), 'ci':(float, None), 'co':(float, None), 'iwa':(float, 4.), \
                                 'owa':(float, 10.), 'oca':(float, 10.), 'fpres':(int,2), 'bw':(float, 0.1), 'Nlam':(int, 3)} }
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
        # Finally, set a private attribute for the number of image plane samples between the center and the OCA
        self.design['Image']['_Nimg'] = int( np.ceil( self.design['Image']['fpres']*self.design['Image']['oca']/(1. - self.design['Image']['bw']/2) ) )
     
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
            self.amplname_solver += "threads{:02d}".format(self.solver['Nthreads'])

        self.ampl_src_fname = self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                              self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".mod" 

        if not hasattr(self, 'telap_fname'):
            self.telap_fname = os.path.join( self.fileorg['TelAp dir'], ("TelAp_full_" + self.amplname_pupil + ".dat") )
        if not hasattr(self, 'fpm_fname'):
            self.fpm_fname = os.path.join( self.fileorg['FPM dir'], "FPM_occspot_M{:03}.dat".format(self.design['FPM']['M']) )
        if not hasattr(self, 'ls_fname'):
            self.ls_fname = os.path.join( self.fileorg['LS dir'], ("LS_full_" + self.amplname_pupil + \
                                          "_{0:02d}D{1:02d}ovsz{2:02d}.dat".format(self.design['LS']['id'], \
                                          self.design['LS']['od'], self.design['LS']['ovsz'])) ) 
                                                              

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
        self.amplname_coron = "APLC_half" 
        self.ampl_src_fname = self.amplname_coron + "_" + self.amplname_pupil + "_" + self.amplname_fpm + "_" + \
                              self.amplname_ls + "_" + self.amplname_image + "_" + self.amplname_solver + ".mod" 
        if not hasattr(self, 'telap_fname') or os.path.split(self.telap_fname)[1].startswith("TelAp_full"):
            self.telap_fname = os.path.join(os.path.split(self.telap_fname)[0], \
                                            os.path.split(self.telap_fname)[1].replace('TelAp_full', 'TelAp_half'))
        if not hasattr(self, 'fpm_fname'):
            self.fpm_fname = os.path.join( self.fileorg['FPM dir'], "FPM_occspot_M{:03}.dat".format(self.design['FPM']['M']) )
        if not hasattr(self, 'ls_fname') or os.path.split(self.ls_fname)[1].startswith("LS_full"):
            self.ls_fname = os.path.join(os.path.split(self.ls_fname)[0], \
                                         os.path.split(self.ls_fname)[1].replace('LS_full', 'LS_half'))

    def write_ampl(self): 
        self.logger.info("Writing the AMPL program for this APLC with halfplane symmetry")
        self.abs_ampl_src_fname = os.path.join(self.fileorg['ampl src dir'], self.ampl_src_fname)
        if not os.path.exists(self.fileorg['ampl src dir']):
           os.mkdir(self.fileorg['ampl src dir'])
        if os.path.exists(self.abs_ampl_src_fname):
            self.logger.warning("Warning: Overwriting the existing copy of {0}".format(self.abs_ampl_src_fname))
        mod_fobj = open(self.abs_ampl_src_fname, "w")

        header = """\
        # AMPL program to optimize a half-plane symmetric APLC
        # Created by {0} with {1} at {2}
        load amplgsl.dll;


        #---------------------

        param pi:= 4*atan(1);

        #---------------------
        param c := {3:.2f};

        #---------------------
        param Rmask := {4:0.3f};
        param rho0 := {5:0.2f};
        param rho1 := {6:0.2f};
        
        #---------------------
        param N := {7};				# discretization parameter (pupil)
        param M := {8};				# discretization parameter (mask)
        
        param Nimg := {9};			# discretization parameter (image)
        param Fmax := {10:0.2f};    # NTZ: If we paramaterize our image plane resolution by fpres = sampling rate at 
                                    # the shortest wavelength, then Nimg should be an integer function of fpres, oca,
        #---------------------      # and bw. This integer is not specified directly by the user, but computed "privately"
        param bw := {11:0.2f};      # by the APLC class constructor.
        param lam0 := 1.;
        param dl := bw*lam0;
        param Nlam := {12};
        
        #---------------------
        param obs := 20;             # NTZ: We could eliminate this section of parameter definitions, 
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
        """.format(getpass.getuser(), os.path.basename(__file__), datetime.datetime.now(), \
                   self.design['Image']['c'], self.design['FPM']['rad'], self.design['Image']['iwa'], self.design['Image']['owa'], \
                   self.design['Pupil']['N'], self.design['FPM']['M'], self.design['Image']['_Nimg'], \
                   self.design['Image']['oca'], self.design['Image']['bw'], self.design['Image']['Nlam'])
        
        load_masks = """\
        #---------------------
        # Loading Pupil
        param PupilFile {{1..2*N,1..N}};
        
        read {{i in 1..2*N,j in 1..N}} PupilFile[i,j] < {0:s};
        close {1:s};
        
        # Loading FPM
        param MaskFile {{1..2*M,1..M}};
        
        read {{i in 1..2*M,j in 1..M}} MaskFile[i,j] < {2:s}; 
        close {3:s};
        
        # Loading Lyot stop
        param LyotFile {{1..2*N,1..N}};
        
        read {{i in 1..2*N,j in 1..N}} LyotFile[i,j] < {4:s};
        close {5:s};
        """.format(self.telap_fname, self.telap_fname, self.fpm_fname, self.fpm_fname, self.ls_fname, self.ls_fname)


        mod_fobj.write( textwrap.dedent(header) )
        mod_fobj.write( textwrap.dedent(load_masks) )

        mod_fobj.close()
        self.logger.info("Wrote %s"%self.abs_ampl_src_fname)
