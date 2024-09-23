'''
Module: ediff.calibration
-------------------------
Semi-automated calibration of electron diffraction patterns.    
'''


# calibration_constant = \
#   ed.calibration.Calculate.from_max_peak(
#       eld_file, xrd_file)
# calibration_constant = \
#   ed.calibration.Calculate.from_max_peak_in_range(
#       eld_file, xrd_file, eld_range_pix=(100,150), xrd_range_q=(3.5,4.5))
# calibration_constant = \
#   ed.calibration.Calculate.from_microscope_constants(
#       voltage_kV=120, camera_lenght_mm=170, camera_pixel_size_um=26.2)
#
# calibration_constant = \
#   ed.calibration.Microscopes.TecnaiVeleta(d=750)
# calibration_constant = \
#   ed.calibration.Microscopes.TecnaiVeleta(d=1000)
 


import numpy as np
from dataclasses import dataclass



class Calculate:
    '''
    A class with functions to *calculate* the SAED calibration constant.
    '''
    
    
    def from_max_peaks(eld_file, xrd_file):
        
        # Function {Calculate.from_max_peak} is a special case
        # We just call {Calculate.from_max_peaks_in_range} for whole ranges,
        # i.e. we call the function without specifying the ranges
        calibration_constant = Calculate.from_max_peaks_in_range(
            eld_file, xrd_file)
        
        # Return calibration constant
        return(calibration_constant)        

    
    def from_max_peaks_in_range(
            eld_file, xrd_file, eld_range_pix=None, xrd_range_q=None):
                        
        # (1) Read ED and XRD diffraction files.
        # (the files are supposed to be in EDIFF format
        # (ED  = 3 columns: pixels, intensity, bkgr-corrected intensity
        # (XRD = 4 columns: 2theta[deg], S[1/A], q[1/A], normalized intensity
        eld = np.loadtxt(eld_file, unpack=True)
        xrd = np.loadtxt(xrd_file, unpack=True)
        
        # (2) Determine ranges, in which we search for peaks.
        # (We search peaks either in whole ranges (all pixels/q for eld/xrd)
        # (or in sub-ranges, if args eld_range_pix and xrd_range_q are given.
        if eld_range_pix is None:
            e_range = eld 
        else:
            emin,emax = eld_range_pix
            e_range = eld[:,(emin<=eld[0])&(eld[0]<=emax)]
        if xrd_range_q is None:
            x_range = xrd
        else:
            xmin,xmax = xrd_range_q
            x_range = xrd[:,(xmin<=xrd[2])&(xrd[2]<=xmax)]
        
        # (3) Get the peak/maximum value for ED and XRD.
        # (The ranges where to search for peaks were determined in prev.step.
        max_eld = float(e_range[0,(e_range[2]==np.max(e_range[2]))])
        max_xrd = float(x_range[2,(x_range[3]==np.max(x_range[3]))])
        print(f'Position of max.peak in ED-profile:  {max_eld:6.4f} in d[pix]')
        print(f'Position of max.peak in XRD-profile: {max_xrd:6.4f} in q[1/A]')
        
        # (4) Calculate calibration constant
        # (the constant converts d[pixels] to q[1/A]
        calibration_constant = max_xrd/max_eld
        print(f'Calibration constant: {calibration_constant:.6f}')
        
        # (5) Return the calculated calibration constant
        return(calibration_constant)
    
    
    def from_microscope_constants(
            voltage_kV, camera_length_mm, camera_pixel_size_um):
        pass
    


@dataclass
class Microscopes:
    '''
    A dataclass with microscopes with the *known* SAED calibration constants.
    '''

    @dataclass
    class TecnaiVeleta:
        pass
    
    
    @dataclass 
    class TecnaiMorada:
        pass
