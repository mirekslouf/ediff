'''
Module ediff.radial
-------------------
The conversion of a 2D powder diffraction pattern
to a 1D powder diffraction pattern = radially averaged intensity distribution.
'''

import numpy as np
import ediff.center

def calc_radial_distribution(arr, center=None, output_file=None):
    """
    Calculate 1D radially averaged distrubution profile
    from 2D diffraction pattern.

    Parameters
    ----------
    arr : 2D-numpy array
        The numpy array which contains the 2D-diffractogram.
    center : tuple/list of two floats, optional, default is None
        The accurate coordinates of the 2D-diffractogram.
        This argument should be determined by ediff.center.CenterLocator
        to get the best results.
        If not given (= if it defaults to None), the center is determined
        an approximate procedure using intensity center,
        without any refinement;
        this is imprecise, especially in case of
        diffraction patterns with a beamstopper.
    output_file : str, optional, default is None
        Name of the output file.
        If given, the calculated 1D profile is saved to *output_file*.

    Returns
    -------
    profile : 2D numpy array containing two rows [R,I]
        * R = radial_distance = dist. from the diffractogram center [pixels]
        * I = intensity = intensities at given distances [arbitrary units]
    """
    
    # (1) Get the center of 2D-diffractogram.
    if center is None:
        # If center argument was not given,
        # we ESTIMATE center as mass/intensity center of the array.
        # We use simple functions from skimage.measure,
        # do not consider central square and do not use refinement.
        # This is usually sufficient for 4D-STEM/PNBD diffractograms,
        # but it fails if the diffraction image contains a beamstopper.
        xc,yc = ediff.center.IntensityCenter.center_of_intensity(
           arr, csquare=50, cintensity=0.8)
    else:
        # If center accurate center coordinates were given,
        # it is much better (and more accurate) and we use them.
        (xc,yc) = center
    
    # (2) Get image dimensions
    (width,height) = arr.shape
    
    # (3) Calculate radial distribution
    # --- (3a) Prepare 2D-array/meshgrid with calculated radial distances
    # (trick 1: the array/meshgrid will be employed for mask
    # (...the meshgrid size = the same as the original array size
    # (...note: when creating meshgrid, xc,yc must be switched
    # (   it is not so clear why, but it was verified repeatedly
    [X,Y] = np.meshgrid(np.arange(width)-yc, np.arange(height)-xc)
    R = np.sqrt(np.square(X) + np.square(Y))
    # --- (3b) Initialize variables
    radial_distance = np.arange(1,np.max(R),1)
    intensity       = np.zeros(len(radial_distance))
    index           = 0
    bin_size        = 1
    # --- (3c) Calcualte radial profile
    # (Gradual calculation of average intenzity
    # (in circles with increasing distance from the center 
    # (trick 2: to create the circles, we will employ mask from trick 1
    for i in radial_distance:
        mask = np.greater(R, i - bin_size) & np.less(R, i + bin_size)
        values = arr[mask]
        intensity[index] = np.mean(values)
        index += 1 
    
    # (4) Save profile to array, save it to file if requested, and return it
    profile = np.array([radial_distance, intensity])
    if output_file: save_radial_distribution(profile, output_file)
    return(profile)

def save_radial_distribution(profile, output_file):
    """
    Save 1D radially averaged distrubution profile to output_file.

    Parameters
    ----------
    profile : 2D numpy array containing two rows [R,I]
        * R = radial_distance = dist. from the diffractogram center [pixels]
        * I = intensity = intensities at given distances [arbitrary units]
    filename : str
        Name of the output file.

    Returns
    -------
    None.
        The output is the radial distribution saved in a file with *filename*. 
    """
    np.savetxt(output_file, np.transpose(profile), fmt='%3d %8.1f')

def read_radial_distribution(filename):
    """
    Read 1D-radially averaged distrubution profile from a TXT-file.

    Parameters
    ----------
    filename : str
        Name of the input file;
        the file is expected to contain two columns [distance, intensity].

    Returns
    -------
    arr : 2D-numpy array
        The array containing two columns: distance, intensity.
    """
    arr = np.loadtxt(filename, unpack=True)
    return(arr)
