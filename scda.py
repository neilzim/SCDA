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

class LyotCoronagraph(object): # Lyot coronagraph base class
    _key_fields = { 'fileorg': ['work dir', 'ampl src dir', 'TelAp dir', 'FPM dir', 'LS dir', 'sol dir', 'eval dir', \
                                'ampl src fname', 'TelAp fname', 'FPM fname', 'LS fname', 'sol fname'], \
                    'solver': ['constr', 'method', 'presolve', 'Nthreads'] }

#    _key_fields = dict([('fileorg', ['ampl src dir', 'TelAp dir', 'FPM dir', 'LS dir', 'sol dir', 'eval dir']),\
#                        ('solver', ['constr', 'method', 'presolve', 'Nthreads'])])

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
        #   intermediate masks, co-eval AMPL programs, solutilons, logs, etc.
        #////////////////////////////////////////////////////////////////////////////////////////////////////
        setattr(self, 'fileorg', {})
        if 'fileorg' in kwargs:
            for namekey, location in kwargs['fileorg'].items():
                if namekey in self._key_fields['fileorg']:
                    self.fileorg[namekey] = location
                    if location is not None:
                        if not os.path.exists(location):
                            self.logger.warning("Warning: The specified location of \"{0}\", \"{1}\" does not exist".format(namekey, location))
                else:
                    self.logger.warning("Warning: Unrecognized field {0} in fileorg argument".format(dirkey))
        # Handle missing directory values
        if 'work dir' not in self.fileorg or self.fileorg['work dir'] is None:
            self.fileorg['work dir'] = os.getcwd()
        for namekey in self._key_fields['fileorg']: # Set other missing directory locations to 'work dir'
            if namekey.endswith('dir') and ( namekey not in self.fileorg or self.fileorg[namekey] is None ):
                self.fileorg[namekey] = self.fileorg['work dir']

        if 'TelAp fname' in self.fileorg and self.fileorg['TelAp fname'] is not None and \
        not os.path.exists(self.fileorg['TelAp fname']) and os.path.exists(self.fileorg['TelAp dir']) and \
        not os.path.isdir(self.fileorg['TelAp fname']):
            try_fname = os.path.join(self.fileorg['TelAp dir'], self.fileorg['TelAp fname']) 
            if os.path.exists(try_fname):
                self.fileorg['TelAp fname'] = try_fname
            else:
                self.logger.warning("Warning: Could not find the specified telescope aperture file \"{0}\" in {1}".format(self.fileorg['TelAp fname'], \
                                    self.fileorg['TelAp dir']))
#        if 'FPM fname' in kwargs:
#            self.fpm_fname = os.path.join(self.fileorg['FPM dir'], kwargs['FPM_fname'])
#            if not os.path.exists(self.fpm_fname):
#                self.logger.warning("Warning: The specified focal plane mask file \"{0}\" does not exist".format(self.fpm_fname))
#        if 'LS fname' in kwargs:
#            self.ls_fname = os.path.join(self.fileorg['LS dir'], kwargs['LS_fname'])
#            if not os.path.exists(self.ls_fname):
#                self.logger.warning("Warning: The specified Lyot stop file \"{0}\" does not exist".format(self.ls_fname))
#        if 'ampl src fname' in kwargs:
#            if os.path.isabs(kwargs['ampl_src_fname']):
#                self.ampl_src_fname = kwargs['ampl_src_fname']
#            else:
#                self.ampl_src_fname = os.path.join(self.fileorg['ampl src dir'], kwargs['ampl_src_fname']) 
                
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
        # Print summary of the set parameters
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

    def write_ampl(self): 
        self.logger.info("Writing the AMPL program for the specified half-plane APLC")
        if not os.path.exists(self.fileorg['ampl src dir']):
           os.mkdir(self.fileorg['ampl src dir'])
        if os.path.exists(self.fileorg['ampl src fname']):
            self.logger.warning("Warning: Overwriting the existing copy of {0}".format(self.fileorg['ampl src fname']))
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
        param Fmax := {11:0.2f};        # NTZ: We paramaterize our image plane resolution by fpres = sampling rate at 
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

        printf {{x in Xs, y in Ys}}: "%15g %15g %15g \\n", x, y, A[x,y] > "{0:s}"
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
        self.logger.info("Wrote %s"%self.fileorg['ampl src fname'])
