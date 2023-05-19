'''
Module ediff.radial
-------------------
The conversion of a 2D powder diffraction pattern
to a 1D powder diffraction pattern = radially averaged intensity distribution.
'''

import numpy as np
from skimage import measure

def calc_radial_distribution(arr, output_file=None):
    """
    Calculate 1D radially averaged distrubution profile
    from 2D diffraction pattern.

    Parameters
    ----------
    arr : 2D-numpy array
        The numpy array which contains the 2D-PNBD pattern.
    output_file : str, optional, default is None
        Name of the output file.
        If given, the calculated 1D profile is saved to *output_file*.

    Returns
    -------
    profile : 2D numpy array containing two rows [R,I]
        * R = radial_distance = dist. from the diffractogram center [pixels]
        * I = intensity = intensities at given distances [arbitrary units]
    """
    # 1) Find center
    # (We employ function from skimage.measure (not from stemdiff.io),
    # (because we want float/non-integer values from the whole array.
    M =  measure.moments(arr,1)
    (xc,yc) = (M[1,0]/M[0,0], M[0,1]/M[0,0])
    # 2) Get image dimensions
    (width,height) = arr.shape
    # 3) 2D-pole/meshgrid with calculated radial distances
    # (trick 1: the array/meshgrid will be employed for mask
    # (it has the same size as the original array for rad.distr.calculation
    [X,Y] = np.meshgrid(np.arange(width)-yc, np.arange(height)-xc)
    R = np.sqrt(np.square(X) + np.square(Y))
    # 4) Initialize variables
    radial_distance = np.arange(1,np.max(R),1)
    intensity       = np.zeros(len(radial_distance))
    index           = 0
    bin_size        = 2
    # 5) Calcualte radial profile
    # (Gradual calculation of average intenzity
    # (in circles with increasing distance from the center 
    # (trick 2: to create the circles, we will employ mask from trick 1
    for i in radial_distance:
        mask = np.greater(R, i - bin_size/2) & np.less(R, i + bin_size/2)
        values = arr[mask]
        intensity[index] = np.mean(values)
        index += 1 
    # 6) Save profile to array, save it to file if requested, and return it
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
