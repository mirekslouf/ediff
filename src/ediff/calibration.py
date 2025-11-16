'''
Module: ediff.calibration
-------------------------

The calibration of electron diffraction patterns:

* The original diffraction pattern
  shows intensities as a function of *distance-in-pixels*.
* The calibrated diffraction pattern
  shows intensities as a function of *distance-in-q-vectors*.
* The term *distance-in-q-vectors*
  means that the distances = magnitudes of q-vector in [1/A].
* This module gives the *calibration constant*,
  which converts *distance-in-pixels* to *distance-in-q-vectors*.
* The final conversion is very simple:
  `distance_in_q = distance_in_pixels * calibration_constant`

How to determine the `calibration_constant`:
    
* In EDIFF, the calibration constant is calculated in
  ediff.calibration.Calculate
* The calculation is simplified for known microscopes in
  ediff.calibration.Microscopes 
* Links to real-life examples are given at GitHub pages:
  https://mirekslouf.github.io/ediff/docs
'''      

import sys
import numpy as np
import scipy.constants
import ediff.io
from dataclasses import dataclass


class Calculate:
    '''
    A class with functions to calculate the *calibration constant*.
    
    * The meaning of *calibration constant* is explained
      above on the top of ediff.calibration module.
    * The calculation of *calibration constant* by means of
      ediff.calibration.Calculate functions is shown below.
     
    >>> # EDIFF :: calculation of calibration constant
    >>> import ediff as ed
    >>>
    >>> # Note: two functions below use ELD and XRD profiles as arguments.
    >>> # - both ELD and XRD profiles are 1D diffraction patterns
    >>> # - ELD = experimental electron diffractogram, radially averaged
    >>> # - XRD = theoretical X-ray diffractogram, calculated in ediff.pcryst 
    >>>
    >>> # (1) Calibration constant from the whole ELD and XRD profiles.
    >>> # (If the max.peak in ELD corresponds to the max.peak in XRD.
    >>> calibration_constant = \\
    >>>     ed.calibration.Calculate.from_max_peaks(ELD, XRD)
    >>> 
    >>> # (2) Calibration constant from selected parts of ELD and XRD profiles.
    >>> # (If the max.peak in ELD does not correspond to the max.peak in XRD,
    >>> # (we define the peak for calibration in specific ranges/regions. 
    >>> calibration_constant = \\
    >>>     ed.calibration.Calculate.from_max_peaks_in_range(
    >>>         ELD, XRD, eld_range=(50,120), xrd_range=(2.0,2.5))
    >>>
    >>> # (3) Calibration constant from known microscope parameters.
    >>> # (3a) Simple case - known microscope with default parameters
    >>> # (more details and examples => ediff.calibration.Microscopes
    >>> my_TEM = ed.calibration.Microscopes.TecnaiVeleta
    >>> calibration_constant = \\
    >>>     ed.calibration.Calculate.from_microscope_parameters(my_TEM)
    >>> # (3b) General case - arbitrary microscope
    >>> # (the four parameters below must be known for given microscope
    >>> calibration_constant = \\
    >>>     ed.calibration.Calculate.from_microscope_parameters(
    >>>         voltage_kV = 120,
    >>>         camera_length_mm = 170,
    >>>         camera_pixel_size_um = 13.2, binning = 2)
    '''
    
    def from_max_peaks(eld_profile, xrd_profile, messages=True):
        '''
        Calibration constant from the *maximal* peaks on ELD and PXRD profiles.

        * ELD/XRD profiles represent 1D electron/X-ray diffraction patterns.
        * More info about ELD/XRD profiles in EDIFF
          => see the docs of ediff.io.read_profile function.

        Parameters
        ----------
        eld_profile : str or numpy.array
            The *eld_profile* (ELD) is
            an electron diffraction profile in EDIFF format.
            It can come as file (if *eld_profile* = str = filename)
            or array (if *eld_profile* = numpy.array).
        xrd_profile : str or numpy.array
            The *xrd_profile* (XRD) is
            an X-rayd diffraction profile in EDIFF format.
            It can come as file (if *xrd_profile* = str = filename)
            or array (if *xrd_profile* = numpy.array).
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
        Calibration constant from the *selected* peaks on ED and PXRD profiles.
        
        * ELD/XRD profiles represent 1D electron/X-ray diffraction patterns.
        * More info about ELD/XRD profiles in EDIFF
          => see the docs of the ediff.io.Profile.read function.
        * The peaks are selected using arguments *eld_range* and *xrd_range*.
            - Both arguments are tuples of two floats = x-ranges.
            - Only the maximal peaks in given ranges are considered.
            - ED range is in [pixels] and PXRD range is given in [q-vectors].

        Parameters
        ----------
        eld_profile : str or numpy.array
            The *eld_profile* is
            an electron diffraction profile in EDIFF format.
            It can come as file (if *eld_profile* = str = filename)
            or array (if *eld_profile* = numpy.array).
        xrd_profile : str or numpy.array
            The *xrd_profile* is
            an X-rayd diffraction profile in EDIFF format.
            It can come as file (if *xrd_profile* = str = filename)
            or array (if *xrd_profile* = numpy.array).
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

        Technical notes
        ---------------
        * ELD profile = 1D radially averaged
          powder electron diffraction pattern
            - in EDIFF, it is obtained from an experimental 2D difractogram
        * XRD profile = 1D powder X-ray diffraction pattern
            - in EDIFF, it is calculated from a standard CIF file
              = Crystallographic Information File
        * EDIFF format of ELD and XRD profiles employed in EDIFF package
            - ELD and XRD profiles can come in the form of files or np.arrays
            - Columns in files <=> rows in np.arrays (we use: *unpack=True*)
            - ELD profile = 3 cols = pixels, intensity, bkgr-corrected-intsty
            - XRD profile = 4 cols = 2theta[deg], S[1/A], q[1/A], norm-intsty
        * EDIFF calculation of ELD and XRD profiles is best seen from examples:
            - https://mirekslouf.github.io/ediff/docs -> worked example
           '''
        
        # (1) Read ED and XRD diffraction profiles.
        # * The profiles are supposed to be either filenames or numpy.arrays
        # * In any case, the filenames or arrays should be in EDIFF format:
        #   ED  = 3 columns: pixels, intensity, bkgr-corrected intensity
        #   XRD = 4 columns: 2theta[deg], S[1/A], q[1/A], normalized-intensity
        eld = ediff.io.Profile.read(eld_profile)
        xrd = ediff.io.Profile.read(xrd_profile)
            
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
    
    
    def from_microscope_parameters(
        microscope=None, voltage_kV=None, camera_length_mm=None, 
        camera_pixel_size_um=None, binning=None, verbose=False):
        '''
        Calibration constant from microscope-specific parameters.
        
        * The calibration constant can be estimated from parameters,
          which are typical of given microscope + camera system.
        * The parameters we need to know are: (i) accelerating_voltage,
          (ii) camera_length, (iii) camera_pixel_size, and (iv) binning.
        * The parameters can be inserted directly OR
          in the form of microscope object, which contains them.
        
        Parameters
        ----------
        microscope : microscope object, optional, default is None
            One of the microscopes defined in ediff.calibration.Microscopes.
        voltage_kV : int, optional, default is None
            Accelerating voltage in [kV].
            If {None} it is taken from {microscope} argument.
        camera_length_mm : float, optional, default is None
            Camera lenght in [mm].
            If {None} it is taken from {microscope} argument.
            WARNING: this must be the *real_camera_length*;
            see the *Note on camera lenght* section below for more details.
        camera_pixel_size_um : float, optional, default is None
            Camera pixel size in [um].
            If {None} it is taken from {microscope} argument.
        binning : int
            Binning=1,2,4,... increases camera_pixel_size (1x,2x,4x,...).
            If {None} it is taken from {microscope} argument.
        verbose : bool, optional, default is False
            If {True}, the function prints outputs to stdnout.

        Returns
        -------
        calibration constant : float
            The multiplicative constant that converts
            ED-profile X-coordinate-in-pixels
            to X-coordinate-in-q-vectors [1/A].
        
        Technical notes
        ---------------
        * This calculation is based on four microscope parameters,
          namely: (i) accelerating_voltage, (ii) camera_length_mm,
          (iii) camera_pixel_size, and (iv) binning.
        * Nevertheless, three of the parameters may change
          from experiment to experiment:
            - *accelerating_voltage_kV* for given TEM
              can be changed by user (although it is not typical).
            - *camera_length_mm* is routinely adjusted
              to see the desired range of difractions;
              the TEM software usually displays some information
              about the nominal camera length.
            - *binning* averages the neighboring pixels of the camera;
              For example, the binning=1,2,4,...
              increases the pixel size 1x,2x,4x,..
        
        Note on camera length
        ---------------------
        * The *nominal camera length* shown by the microscope (*D*)
          may be different from *real camera length*,
          which is the argument (camera_lenght_um) used in the calculation.
            - Possible reasons: Exact amera position (bottom x upper camera),
              projective lenses (which may change the magnification
              of the diffraction pattern.
            - Real-life solution: we need to know real *camera_length_mm*
              for each theoretical/nominal/TEM-software-displayed
              cammera length *D*!
        * Note: For known/calibrated microscopes
          (that are defined in ediff.calibration.Microscopes)
          we insert the *nominal camera length*, because we know
          the conversion coefficient from nominal to real camera length.
          Here, for the arbitrary microscope, we have to insert
          the *real camera length* in order to get the correct result.
        '''
        
        # (0) If {microscope} argument was given,
        # take MISSING parameters for the calculation from {microscope}.
        if (voltage_kV is None) and (microscope is not None):
            voltage_kV = microscope.voltage_kV
        else:
            sys.exit('Microscope params - missing accelerating voltage!')
        if (camera_length_mm is None) and (microscope is not None): 
            camera_length_mm = microscope.D * microscope.D_to_CL_coefficient
        else:
            sys.exit('Microscope params - camera length is missing!')
        if (camera_pixel_size_um is None) and (microscope is not None):
            camera_pixel_size_um = microscope.camera_pixel_size_um
        else:
            sys.exit('Microscope params - camera pixel size is missing!')
        if (binning is None) and (microscope is not None):
            binning = microscope.binning
        else:
            sys.exit('Microscope params - no info about binning!')
        
        # (1) Calculate relativistic wavelength of electrons
        # (this is needed to convert CL to CC in the next step
        Lambda = Utils.electron_wavelength(voltage_kV)
        
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
            Utils.calc_final_calibration_constants(
                camera_pixel_size_um = camera_pixel_size_um,
                binning = binning,
                camera_constant_mmA = CC)
        
        # (4) Print the calculated values if requested
        if verbose:
            Utils.print_final_calibration_constants(
                Lambda, CL, CC,
                R_calibration, S_calibration, q_calibration)

        # (5) Return the calculated calibration constant
        # (We want the calibration constant in q [1/A]
        return(q_calibration)


class Microscopes:
    '''
    A dataclass containing know/calibrated microscopes.

    * The microscopes are defined as subclasses.
    * More precisely, the microscopes are sub-dataclasses.
    * The calibrated microscopes can be used to determine calibration constant.
    
    >>> How to calculate calibration constant for known/calibrated microscope
    >>> import ediff as ed
    >>>
    >>> # (1) Define the microscope + its parameters
    >>> # (a) Simple case - all parameters at their default/standard values
    >>> my_TEM = ed.calibration.Microscopes.TecnaiVeleta()
    >>> # (b) Intermediate case - standard values + camera lenght changed
    >>> my_TEM = ed.calibration.Microscopes.TecnaiVeleta(D=750)
    >>> # (c) The most comlex case - multiple parameters have changed
    >>> my_TEM = ed.calibration.Microscopes.TecnaiVeleta(
    >>>     D=750, voltage_kV=80, binning=4)
    >>>
    >>> # (2) Calculate the calibration constant
    >>> const = ed.calibration.Calculate.from_microscope_parameters(my_TEM)
    
    Technical notes:
    
    * Each calibrated {microscope + camera} system
      is defined by one ediff.calibration.Microscopes subclass.
    * Example: {Tecnai microscope + Veleta camera}
      is represented by ediff.calibration.Microscopes.TecnaiVeleta subclass.
    * The final calculation of the calibration constant is made by the
      general function: ediff.calibration.Calculate.from_microscope_parameters,
      but in this case the arguments of the function are pre-defined
      within each subclass (such as TecnaiVeleta, TecnaiMorada...).
    '''
    
    # TECHNICAL NOTE
    # ==============
    # * the previous version used dataclasses with a property
    # * the property yielded the calibration constant within the dataclass
    # * the problem was that dataclasses CANNOT use self during initialization
    # * the solution how to calculate within dataclasses => properties
    # * the sources/hints/notes how this can be done are given/kept below
    #
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
        The default parameters for
        the microscope/camera = Tecnai/Veleta3G.
        
        Parameters
        ----------
        voltage_kV : float, default is 120
            The accelerating voltage, typical of given system.
        D : float, default is 1000
            The *nominal camera length* from the control software.
        D_to_CL_coefficient : float, default is 1/4.5
            The coefficient to convert
            *nominal camera lenght* to *real_camera_length* in mm.
            This coefficient is specific for given system
            as the *nominal* and *real* camera length may be different.
        camera_pixel_size_um : float, default is 13.2
            The camera pixel size = a constant for given camera.
        binning : int, default is 2
            The binning, which is typical of given system.
        '''
        
        voltage_kV           : float = 120
        D                    : float = 1000
        D_to_CL_coefficient  : float = 1 / 4.5
        camera_pixel_size_um : float = 13.2
        binning              : int   = 2


    @dataclass 
    class TecnaiMorada:
        '''
        The default parameters for
        the microscope/camera = Tecnai/Morada.
        
        Parameters
        ----------
        voltage_kV : float, default is 120
            The accelerating voltage, typical of given system.
        D : float, default is 660
            The *nominal camera length* from the control software.
        D_to_CL_coefficient : float, default is 1/4
            The coefficient to convert
            *nominal camera lenght* to *real_camera_length*.
            This coefficient is specific for given system
            as the *nominal* and *real* camera length may be different.
        camera_pixel_size_um : float, default is 9.0
            The camera pixel size = a constant for given camera.
        binning : int, default is 4
            The binning, which is typical of given system.
        '''
        voltage_kV           : float = 120
        D                    : float = 660
        D_to_CL_coefficient  : float = 1 / 4
        camera_pixel_size_um : float = 9.0
        binning              : int   = 4
        

class Utils:
    '''
    Utilities for the calculation of calibration constants.
    '''
    
    
    def electron_wavelength(U):
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
        return(float(Lambda))


    def calc_final_calibration_constants(
            camera_pixel_size_um, binning, camera_constant_mmA):
        '''
        Calculate final calibration constants
        from known CC and camera_pixel_size.

        * Once we know *CC*, *camera_pixel_size*, and *binning*,
          the calculation of the final calibration constants
          is surprisingly easy.
        * The full justification of the calculation
          is given in the comments of the source code below.
        
        Parameters
        ----------
        camera_constant_mmA : float
            The camera constant (from R*d = CL*Lambda = CC) in [mmA].
        camera_pixel_size_um : float
            The real dimension of one pixel of the camera (detector) in [um].
        binning : int
            Binning = 1,2,4,... increases pixel size (1x,2x,4x,...)

        Returns
        -------
        R_calibration, S_calibration, q_calibration : three floats
            * R_calibration = the camera pixel size in [mm].
            * S_calibration = the camera pixel size in S-vector units [1/A].
            * q_calibration = the camaera pixel size in q-vector units [1/A].
        '''

        # What do we want to calculate?
        #  => three calibration constants = pix_size, S-units, q-units        
        # Brief justification:
        #  (a) Camera Equation: R*d = CL*Lambda = CC
        #  (b) Bragg's Law: 2*d*sin(theta) = Lamda => S*d = 1 => S = 1/d
        #  (c) Combine (a)+(b) for a known R[um] => here: R = size of 1 pixel
        #      in general : R*d = CC => S = 1/d = R/CC
        #      here/below : S_calibration[1/A] = R_calibration[mm] / CC [mm*A]
        # -----
        # Note: camera_pixel_size is influenced by binning
        
        # Calculate the calibration constants acc.to justification above
        R_calibration = camera_pixel_size_um * binning / 1000
        S_calibration = R_calibration / camera_constant_mmA
        q_calibration = 2*np.pi * S_calibration

        # Return the three calibration constants
        return(R_calibration, S_calibration, q_calibration)

    
    def print_final_calibration_constants(
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
        Calibrate and normalize ELD profile.
        
        * Assumption: The ELD profile is in EDIFF format as described below.
        
        Parameters
        ----------
        eld_profile : numpy.array
            The original/non-calibrated ELD profile in EDIFF format.
            The details are in the *Technical notes* section below.
        calibration_constant : float
            The calibration constant,
            which converts *distance-in-pixels* to *distance-in-q-vectors*.
            The details are in the *Technical notes* section below.
            
            
        Returns
        -------
        eld_profile : numpy.array
            The calibrated and normalized ELD profile (in EDIFF format).

        Technical notes
        ---------------
        * ELD profile in EDIFF format
          = text file with 3 columns + comments (starting with #).
            - EDIFF file columns
              = {pixel}, {intensity}, {bkg-corrected-intensity}.
            - The {pixel} column
              = {distance-in-pixels-from-the-diffractogram-center}.
            - More details concerning ELD profiles
              = the documentation of ediff.io.Profile.read function.
        * The *calibration constant*
          = the final constant for the conversion/calibration.
            - It is explained in more detail at the initial desription
              of ediff.calibration module.
        '''
        
        # X-data = calibrate = convert [pixels] to [q-vectors]. 
        eld_profile[0] = eld_profile[0] * calibration_constant
        
        # Y-data = normalize to 1 (for both raw and bkgr-correctet intensities)
        eld_profile[1] = eld_profile[1]/np.max(eld_profile[1])
        eld_profile[2] = eld_profile[2]/np.max(eld_profile[2])
        
        # Return the calibrated profile
        return(eld_profile)
