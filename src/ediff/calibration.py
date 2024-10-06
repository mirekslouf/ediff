'''
Module: ediff.calibration
-------------------------

The module performs calibration of electron diffraction patterns:

* The original diffraction pattern
  shows intensities as a function of *distance-in-pixels*.
* The calibrated diffraction pattern
  shows intensities as a function of *distance-in-q-vectors*.
* The term *distance-in-q-vectors*
  means the distance equals to the magnitude of q-vector in [1/A].
* This module gives the *calibration constant*,
  which converts *distance-in-pixels* to *distance-in-q-vectors*.
* The final conversion is very simple:
  `distance_in_q = distance_in_pixels * calibration_constant`

Various ways how to determine the calibration constant:

* The calibration constant can be determined
  by means of various functions in ediff.calibration module.
* The calibration is not difficult,
  but we need at least one of the following:
    - ELD and XRD profiles
      (which can be calculated by ediff)
    - some microscope constants
      (voltage, camera_lenght, camera_pix_size)
    - a pre-calibrated microscope
      (stored in the sub-dataclass ediff.calibration.Microscopes)
* If the constants are not known,
  you can always use the EDIFF-calculated ELD and XRD profiles.
* The short examples are below,
  more details are in the *worked examples* coming with this documentation.

>>> # Standard import of ediff
>>> import ediff as ed
>>>
>>> # (1) Calibration constant from the whole ELD and XRD profiles.
>>> # (If the max.peak in ELD corresponds to the max.peak in XRD
>>> calibration_constant = ed.calibration.Calculate.from_max_peaks(ELD, XRD)
>>> 
>>> # (2) Calibration constant from selected parts of ELD and XRD profiles.
>>> # (If the max.peak in ELD corresponds to low peak in XRD or vice versa
>>> calibration_constant = ed.calibration.Calculate.from_max_peaks_in_range(
>>>     ELD, XRD, eld_range=(50,120), xrd_range=(2.0,2.5))
>>>
>>> # (3) Calibration constant from known microscope parameters.
>>> # (The three parameters below are often known or you can find them.
>>> calibration_constant = ed.calibration.Calculate.from_microscope_constants(
>>>     voltage_kV = 120,
>>>     camera_length_mm = 170,
>>>     camera_pixel_size_um = 26.2)
>>>
>>> # (4) Calibration constant from known/calibrated microscope.
>>> # (The calibrated microscopes are in: ediff.calibration.Microscopes
>>> calibration_constant = ed.calibration.Microscopes.TecnaiVeleta(D=750)
'''



import numpy as np
import scipy.constants
from dataclasses import dataclass



class Calculate:
    '''
    A class with functions to *calculate* the SAED calibration constant.
    '''
    
    
    def from_max_peaks(eld_profile, xrd_profile, messages=True):
        '''
        Calibration constant from the *maximal* peaks on ED and PXRD profiles.

        Parameters
        ----------
        eld_profile : str or numpy.array
            The *eld_profile* is
            an electron diffraction profile in EDIFF format.
            It can come as file (if *eld_profile* = str = filename)
            or array (*eld_profile* = numpy.array).
            Both file or array must be in EDIFF format,
            i.e. containing the following 3 columns:
            pixels, intensity, bkgr-corrected intensity.
        xrd_profile : str or numpy.array
            The *xrd_profile* is
            an X-rayd diffraction profile in EDIFF format.
            It can come as file (if *xrd_profile* = str = filename)
            or array (*xrd_profile* = numpy.array).
            Both file or array must be in EDIFF format,
            i.e. containing the following 4 columns:
            2theta[deg], S[1/A], q[1/A], normalized-intensity.
        messages : bool, optional, default is True
            If *messages* = True,
            print some information
            and the final calibration constant to stdout.
        
        Returns
        -------
        calibration constant : float
            The multiplicative constant that converts
            ED-profile X-coordinate-in-pixels
            to X-coordinate-in-q-vectors [1/A].
        '''
        
        # Function {Calculate.from_max_peak} is a special case
        # of more general function {Calculate.from_max_peaks_in_range},
        # which we call without specifying the ranges => to use whole ranges.
        calibration_constant = Calculate.from_max_peaks_in_range(
            eld_profile, xrd_profile, messages=messages)
        
        # Return calibration constant
        return(calibration_constant)        

    
    def from_max_peaks_in_range(
            eld_profile, xrd_profile,
            eld_range=None, xrd_range=None, messages=True):
        '''
        Calibration constant from the *selected* peaks on ED and PXRD profile.
        
        The peaks are selected using arguments *eld_range* and *xrd_range*.
        Both arguments are tuples of two floats = x-ranges.
        Only the maximal peaks in given ranges are considered.
        ED range is defined in [pixels] and PXRD range is given in [q-vectors].

        Parameters
        ----------
        eld_profile : str or numpy.array
            The *eld_profile* is
            an electron diffraction profile in EDIFF format.
            It can come as file (if *eld_profile* = str = filename)
            or array (*eld_profile* = numpy.array).
            Both file or array must be in EDIFF format,
            i.e. containing the following 3 columns:
            pixels, intensity, bkgr-corrected intensity.
        xrd_profile : str or numpy.array
            The *xrd_profile* is
            an X-rayd diffraction profile in EDIFF format.
            It can come as file (if *xrd_profile* = str = filename)
            or array (*xrd_profile* = numpy.array).
            Both file or array must be in EDIFF format,
            i.e. containing the following 4 columns:
            2theta[deg], S[1/A], q[1/A], normalized-intensity.
        eld_range : tuple of two floats, optional, default is None
            The x-range in 1D ED profile,
            in which we should search for the maximal peak.
            The ED x-range is given in [pixels].
        xrd_range : tuple of two floats, optional, default is None
            The x-range in 1D XRD profile,
            in which we should search for the maximal peak.
            The XRD x-range is given in [q-vectors].
        messages : bool, optional, default is True
            If *messages* = True,
            print some information
            and the final calibration constant to stdout.
            
        Returns
        -------
        calibration constant : float
            The multiplicative constant that converts
            ED-profile X-coordinate-in-pixels
            to X-coordinate-in-q-vectors [1/A].
        '''
        
        # (1) Read ED and XRD diffraction profiles.
        # * The profiles are supposed to be either filenames or numpy.arrays
        # * In any case, the filenames or arrays should be in EDIFF format:
        #   ED  = 3 columns: pixels, intensity, bkgr-corrected intensity
        #   XRD = 4 columns: 2theta[deg], S[1/A], q[1/A], normalized-intensity
        eld = Utils.read_profile(eld_profile)
        xrd = Utils.read_profile(xrd_profile)
            
        # (2) Determine ranges, in which we search for peaks.
        # (We search peaks either in the whole x-ranges
        # (or in sub-ranges, if the args eld_range and xrd_range are given.
        # ! X-ranges for ED and XRD are given in [pixel] and [q], respectively.
        if eld_range is None:
            e_range = eld 
        else:
            emin,emax = eld_range
            e_range = eld[:,(emin<=eld[0])&(eld[0]<=emax)]
        if xrd_range is None:
            x_range = xrd
        else:
            xmin,xmax = xrd_range
            x_range = xrd[:,(xmin<=xrd[2])&(xrd[2]<=xmax)]
        
        # (3) Get the peak/maximum value for ED and XRD.
        # (The ranges where to search for peaks were determined in prev.step.
        max_eld = float(e_range[0,(e_range[2]==np.max(e_range[2]))])
        max_xrd = float(x_range[2,(x_range[3]==np.max(x_range[3]))])
        print(f'Position of max.peak in ELD-profile: {max_eld:6.4f} in d[pix]')
        print(f'Position of max.peak in XRD-profile: {max_xrd:6.4f} in q[1/A]')
        
        # (4) Calculate calibration constant
        # (the constant converts d[pixels] to q[1/A]
        calibration_constant = max_xrd/max_eld
        print(f'Calibration constant: {calibration_constant:.6f}')
        
        # (5) Return the calculated calibration constant
        return(calibration_constant)
    
    
    def from_microscope_constants(
            voltage_kV, camera_length_mm, camera_pixel_size_um, messages=True):
        '''
        Calibration constant from microscope-specific constants.
        
        * The calibration constant can be estimated from certain constants,
          which are typical of given microscope + camera system.
        * The constants we need to know are: (i) accelerating_voltage,
          (ii) camera_length, and (iii) camera_pixel_size.
        * Warning: All three 'constants' may change for given experiment,
          as explained below in Technical notes section.

        Parameters
        ----------
        voltage_kV : int
            Accelerating voltage in [kV].
        camera_length_mm : float
            Camera lenght in [mm].
        camera_pixel_size_um : TYPE
            Camera pixel size in [um].
        messages : bool, optional, default is True
            If {True}, the function prints outputs to stdnout.

        Returns
        -------
        calibration constant : float
            The multiplicative constant that converts
            ED-profile X-coordinate-in-pixels
            to X-coordinate-in-q-vectors [1/A].
        
        Technical notes
        ---------------
        * This calculation is based on three microscope "constants",
          namely: (i) accelerating_voltage, (ii) camera_length,
          and (iii) camera_pixel_size.
        * All three "constants" may change from experiment to experiment.
            - *accelerating_voltage_kV* for given TEM
              can be changed (although it is not done frequently).
            - *camera_length_mm* is usually adjusted
              to see the desired range of difractions;
              the TEM software usually displays an info on the camera length.
            - *camera_pixel_size_um* does not change physically,
              but it increases 2x,4x ... for binning=2,4 ...
        * IMPORTANT: The camera length shown by the microscope
          may not be the correct camera lenght needed for the calculation.
            - Possible reason: Different real camera position,
              for example the difference between bottom x upper camera.
            - Real-life solution: we need to know real camera length
              for each theoretical/TEM-software-displayed cammera length!
        '''
        
        # (1) Calculate relativistic wavelength of electrons
        # (this is needed to convert CL to CC in the next step
        Lambda = Utils.electron_wavelenght(voltage_kV)
        
        # (2) Determine CL and CC
        # * CL = CameraLenght
        #   CL should be known from exp => argument camera_length_mm
        # * CC = Camera Constant
        #   CC = CL * Lamda => CC[mmA] = CL[mm]*Lamda[A]
        #   The CC-CL => Camera equation: R*d = Lambda*CL = CC
        CL = camera_length_mm
        CC = camera_length_mm * Lambda
        
        # (3) Calculate final calibration constants.
        # * The final calibration constants are:
        #   R_calibration = pixel size in [mm]
        #   S_calibration = k_calibration = pixel size in S [1/A]
        #   q_calibration ............... = pixel size in q = 2*pi*S [1/A]
        # * To calculate all three constants, we need just two parameters:
        #   camera_pixel_size  = real size of the camera pixels
        #   CC = camera_length = camera length, determined above
        # * The calculation is easy
        #   but it is performed in a separate function,
        #   because it is employed also in calibration.Microscopes class.
        # * The function contains
        #   more details, including full justification of the calculations.
        R_calibration, S_calibration, q_calibration = \
            Utils.final_calibration_constants(
                camera_pixel_size_um = camera_pixel_size_um,
                camera_constant_mmA = CC)
        
        # (4) Print the calculated values if requested
        if messages:
            Utils.print_calibration_constants(
                Lambda, CL, CC,
                R_calibration, S_calibration, q_calibration)

        # (5) Return the calculated calibration constant
        # (We want the calibration constant in q [1/A]
        return(q_calibration)
    


@dataclass
class Microscopes:
    '''
    A dataclass with calibrated microscopes to determine calibration constants.
    '''
    # Dataclasses
    #   https://docs.python.org/3/library/dataclasses.html
    # Dataclasses in PDoc:
    #   https://mirekslouf.github.io/myimg/docs/pdoc.html/myimg/settings.html
    # Technical notes/tricks:
    #   * @dataclass cannot use self during the definition/initialization.
    #     Reason: During a dataclass initialization, self is undefined yet.
    #     Solution: Define additional properties later, using @property.
    #   * @property decorator converts methods to properties.
    #     More precisely, a method can be used as a property ...
    #     ... i.e: class.method() => class.method (without parentheses).
    
    
    @dataclass
    class TecnaiVeleta:
        '''
        Microscope:Tecnai + camera:Veleta3G.
        '''
        D                    : float = 1000
        voltage_kV           : float = 120
        camera_pixel_size_um : float = 26.2
        D_to_CL_coefficient  : float = 1 / 4.5
        messages             : bool  = True

    
        @property 
        def calibration_constant(self):
            
            # The following calculation is almost analogous to
            #  calibration.Calculate.from_microscope_constants - see above.
            # There is just one small difference:
            #  Step(2) - CL is calculated by means of D_to_CC_coefficient.
            #  Reason: D given by the microscope may not be the real CL[mm].
            
            # (1) Calculate relativistic wavelength
            Lambda = Utils.electron_wavelenght(self.voltage_kV)
            
            # (2) Determine CL and CC
            # * CL = CameraLenght
            #   CL is estimated from known distance D from the microscope
            #   CL = D * D_to_CL_coefficient (coeff. typical of given device)
            # * CC = Camera Constant
            #   CC = CL * Lamda => CC[mmA] = CL[mm]*Lamda[A]
            #   The CC-CL => Camera equation: R*d = Lambda*CL = CC
            CL = self.D * self.D_to_CL_coefficient
            CC = CL * Lambda
            
            # (3) Calculate the final calibration constants
            # (see  Utils.final_calibration_constants for the justification
            R_calibration, S_calibration, q_calibration = \
                Utils.final_calibration_constants(
                    camera_pixel_size_um = self.camera_pixel_size_um, 
                    camera_constant_mmA = CC)
            
            # (4) Print calibrations if requested
            if self.messages:
                Utils.print_calibration_constants(
                    Lambda, CL, CC,
                    R_calibration, S_calibration, q_calibration)
                
            # (5) Return the final q_calibration constant
            return(q_calibration)


    @dataclass 
    class TecnaiMorada:
        '''
        Microscope:TecnaiSpirit + camera:Morada.
        '''
        pass


        @property 
        def calibration_constant(self):
            pass



class Utils:
    '''
    Utilities for the calculation of calibration constants.
    '''

    
    def read_profile(profile):
        '''
        Read the ED or XRD profile in EDIFF format.

        Parameters
        ----------
        profile : str or numpy.array
            
            * If profile = str,
              we assume that it is the filename
              of the file with ELD or XRD profile in EDIFF format.
            * If profile = numpy.array,
              we assume that it is the 2D-array
              containing ELD or XRD profile in EDIFF format.
            * See section *Technical notes* below
              for explanation of the EDIFF format of the ELD and XRD profiles.

        Returns
        -------
        profile : 2D numpy.array
            The array representing ELD or XRD profile in EDIFF format.
            See section *Technical notes* below
            for explanation of the EDIFF format of the ELD and XRD profiles.  

        Technical notes
        ---------------
        * EDIFF format of ELD and XRD profiles.
            * EDIFF format = format employed in EDIFF package.
            * ELD and XRD profiles can come in the form of files or arrays.
            * ELD profile = 3 cols = pixels, intensity, bkgr-corrected-intsty
            * XRD profile = 4 cols = 2theta[deg], S[1/A], q[1/A], norm-intsty
        * Note1: ELD = electron diffraction, XRD = X-ray diffraction.
        * Note2: Columns in files <=> rows in arrays (unpack=True).
        '''
        if type(profile)==np.ndarray:
            return(profile)
        else:
            profile = np.loadtxt(profile, unpack=True)
            return(profile)
        
    def electron_wavelenght(U):
        '''
        Calculate relativistic wavelenght of accelerated electrons.
    
        Parameters
        ----------
        U : float
            Accelerating voltage [kV].
    
        Returns
        -------
        Lambda : float
            Relativistic wavelenght of electrons,
            which were accelerated with the voltage *U*[kV]
            in an electron microscope.
        '''
        # The formula below gives Lambda = wavelength in [A]
        # (justification: textbooks of physics
        # (my source: M:/MIREK.PRE/0_EM.PY/1_ELN-IN-EM/1eln_v-wl_v2.nb.html
        
        # Convert U from [kV] to [V]
        U = U * 1000
        
        # Collect constants from scipy
        h = scipy.constants.h
        c = scipy.constants.c
        e = scipy.constants.e 
        m = scipy.constants.m_e
        
        # Calculate and return relativistic wavelenght of electrons
        Lambda = h/np.sqrt(2*m*e*U) * (1/(np.sqrt(1+(e*U)/(2*m*c**2)))) * 1e10
        return(Lambda)


    def final_calibration_constants(camera_pixel_size_um, camera_constant_mmA):
        '''
        Final calibration constants from known CC and camera_pixel_size.

        * Once we know *CC* and *camera_pixel_size*,
          the final calibration constants are calculated surprisingly easy.
        * The full justification of the calculation
          is given in the comments of the source code below.
        
        Parameters
        ----------
        camera_constant_mmA : float
            The camera constant (from R*d = CL*Lambda = CC) in [mmA].
        camera_pixel_size_um : float
            The real dimension of one pixel of the camera (detector) in [um].

        Returns
        -------
        R_calibration, S_calibration, q_calibration : three floats
            * R_calibration = the camera pixel size in [mm].
            * S_calibration = the camera pixel size in S-vector units [1/A].
            * q_calibration = the camaera pixel size in q-vector units [1/A].
        '''

        # What do we want to calculate?
        #  => three calibration constants = pixel size in [um],S-units,q-units        
        # Brief justification:
        #  (a) Camera Equation: R*d = CL*Lambda = CC
        #  (b) Bragg's Law: 2*d*sin(theta) = Lamda => S*d = 1 => S = 1/d
        #  (c) Combine (a)+(b) for a known R[um] => here: R = size of 1 pixel
        #      in general : R*d = CC => S = 1/d = R/CC
        #      here/below : S_calibration[1/A] = R_calibration[mm] / CC [mm*A]
        R_calibration = camera_pixel_size_um / 1000
        S_calibration = R_calibration / camera_constant_mmA
        q_calibration = 2*np.pi * S_calibration

        # Return the three calibration constants
        return(R_calibration, S_calibration, q_calibration)

    
    def print_calibration_constants(
            Lambda, CL, CC, R_calibration, S_calibration, q_calibration):
        '''
        Print all constants employed in calibraton to stdout in nice form.

        Parameters
        ----------
        Lambda : float
            The relativistic wavelenght of electrons in [A].
        CL : float
            The camera lenght (from R*d = CL*Lambda) in [mm].
        CC : float
            The camera constant (from R*d = CL*Lambda = CC) in [mmA].
        R_calibration : float, 
            The camera pixel size in [mm].
        S_calibration : float
            The camera pixel size in S-vector units [1/A].
        q_calibration : TYPE
            The camera pixel size in q-vector units [1/A].

        Returns
        -------
        None
            The arguments are just nicely printed in stdout.

        '''
        print(f'Lambda(relativistic)         : {Lambda:.5f} [A]')
        print(f'CL = CameraLength            : {CL:.1f} [mm]')
        print(f'CC = CameraConstant          : {CC:.5f} [mmA]')
        print(f'R_calibration = pixel_size   : {R_calibration:.5f} [mm]')
        print(f'S_calibration = pixel_size_S : {S_calibration:.5f} [1/A]')
        print(f'q_calibration = pixel_size_q : {q_calibration:.5f} [1/A]')
        
    
    def calibrate_and_normalize_eld_profile(eld_profile, calibration_constant):
        '''
        Calibrate and normalize ELD profile (in EDIFF format)

        Parameters
        ----------
        eld_profile : numpy.array
            The original ELD profile in EDIFF format.
            
        calibration_constant : float
            DESCRIPTION
            
        Returns
        -------
        eld_profile : numpy.array
            The calibrated and normalized ELD profile (in EDIFF format).
        '''
        
        # X-data = calibrate = convert [pixels] to [q-vectors]. 
        eld_profile[0] = eld_profile[0] * calibration_constant
        
        # Y-data = normalize to 1 (for both raw and bkgr-correctet intensities)
        eld_profile[1] = eld_profile[1]/np.max(eld_profile[1])
        eld_profile[2] = eld_profile[2]/np.max(eld_profile[2])
        
        # Return the calibrated profile
        return(eld_profile)
