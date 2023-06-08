'''
Module ediff.center
-------------------
Find center of 2D diffraction pattern.    
'''

import numpy as np
import skimage as sk

def central_square(arr, csquare, xcenter=None, ycenter=None):
    '''
    Return central square from an array
    '''
    xsize,ysize = arr.shape
    # If center of was not given, take geometrical center
    # (for array selections/slicing, we need integers => round, //
    xc = round(xcenter) or xsize // 2
    yc = round(ycenter) or ysize // 2
    # Half of the central square
    # (for array selections/slicing, we need integers => //
    half_csquare = csquare // 2
    # Create sub-array = just central square around xc,yc
    arr2 = arr[
        xc-half_csquare:xc+half_csquare,
        yc-half_csquare:yc+half_csquare].copy()
    return(arr2)

def center_of_intensity(arr, csquare=20, cintensity=0.8):
    '''
    Find center of intensity/mass of an array.
    
    Parameters
    ----------
    arr : 2D-numpy array
        The array, whose intensity center will be determined.
    csquare : int, optional, default is 20
        The size/edge of the square in the (geometrical) center.
        The intensity center will be searched only within the central square.
        Reasons: To avoid other spots/diffractions and
        to minimize the effect of possible intensity assymetry around center. 
    cintensity : float, optional, default is 0.8
        The intensity fraction.
        When searching the intensity center, we will consider only
        pixels with intensity > max.intensity.
        
    Returns
    -------
    xc,yc : float,float
        XY-coordinates of the intensity/mass center of the array.
        Round XY-coordinates if you use them for image/array calculations.
    '''
    # Get image/array size
    xsize,ysize = arr.shape
    # Calculate borders around the central square
    xborder = (xsize - csquare) // 2
    yborder = (ysize - csquare) // 2
    # Create central square = cut off the borders
    arr2 = arr[xborder:-xborder,yborder:-yborder].copy()
    # In the central square, set all values below cintenstity to zero
    arr2 = np.where(arr2>np.max(arr2)*cintensity, arr2, 0)
    # Calculate 1st central moments of the image
    M = sk.measure.moments(arr2,1)
    # Calculate the intensity center = centroid according to www-help
    (xc,yc) = (M[1,0]/M[0,0], M[0,1]/M[0,0])
    # We have centroid of the central square => recalculate to whole image
    (xc,yc) = (xc+xborder,yc+yborder)
    # Return the final center
    return(xc,yc)

def center_of_powder_diffractogram(arr):
    '''
    TODO - Pavlina.

    Parameters
    ----------
    arr : 2D-numpy array
        The array containing 2D-diffractogram, whose center will be determined.

    Returns
    -------
    xc,yc : float,float
        XY-coordinates of the intensity/mass center of the array.
        Round XY-coordinates if you use them for image/array calculations.
    '''
    pass