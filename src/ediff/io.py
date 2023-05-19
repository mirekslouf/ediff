'''
Module ediff.io
---------------
Input/output functions for package ediff.    
'''

from PIL import Image

import numpy as np
import matplotlib.pyplot as plt

import ediff.radial

def read_image(image_name, itype=None):
    '''
    Read grayscale image into 2D numpy array.
    
    Parameters
    ----------
    image_name : string or pathlib object
        Name of image that should read into numpy 2D array.
    itype: string ('8bit'  or '16bit')
        type of the image: 8 or 16 bit grayscale    
        
    Returns
    -------
    2D numpy array
    '''
    img = Image.open(image_name)
    if itype=='8bit':
        arr = np.asarray(img, dtype=np.uint8)
    else:
        arr = np.asarray(img, dtype=np.uint16)
    return(arr)

def plot_radial_distributions(
        data_to_plot, xlimit, ylimit, output_file=None):
    """
    Plot one or more 1D-radial distrubution files in one graph.

    Parameters
    ----------
    data_to_plot : 2D-list 
        * list with several rows containing [data, linestyle, label]
        * data = data for plotting - they can be one of the following:
            - PNG filename = str, a PNG-file = 2D diffraction pattern
            - TXT filename = str, a text file = 1D diffraction profile
            - 2D-array = a numpy array, containg 2D diffraction pattern
            - 1D-array = a numpy array, containing 1D diffraction profile
            - Note1: 2D-pattern = a square image/array with intensities
            - Note2: 1D-profile = a text file/array with two cols/rows = [R,I],
              where R = distance from center, I = diffraction intensity
        * linestyle = matplotlib.pyplot format, such as 'r-' (red line)
        * label = name of the data, which will appear in the plot legend
    xlimit : int
        maximum of the X-axis
    ylimit : int
        maximum of the Y-axis
    output_file : int, optional, default=None
        Name of the output file;
        if the *output* argument is given,
        the plot is not only shown on screen, but also saved in *output* file. 

    Returns
    -------
    Nothing
        The output is the plot on screen
        (and in *output file* if the *output* argument was given).
    
    Technical note
    --------------
    This function is quite flexible.
    It can plot one radial distribution or more.
    It can take data from PNG-files, TXT-files, 2D-arrays and 1D-arrays.
    This makes the code a bit more complex, but it is convenient for the user.
    A fast comparison of three 1D-distributions from three 2D-diffractograms:
    
    >>> ediff.io.plot_radial_distributions(
    >>>     data_to_plot = [
    >>>         ['sum_all_16bit.png', 'k:',  'All data'],
    >>>         ['sum_f_16bit.png',   'b--', 'F data'],
    >>>         ['sum_fd_16bit.png',  'r-',  'FD data']]
    >>>     xlimit=200, ylimit=300,
    >>>     output_file='sums_final_1d.png')
    """
    # Read radial distribution files
    n = len(data_to_plot)
    rdist = data_to_plot
    # Plot radial distribution files
    for i in range(n):
        # Read data
        data = rdist[i][0]
        if type(data) == str:  # Datafile
            if data.lower().endswith('.png'):  # ....PNG file, 2D-diffractogram
                arr = read_image(data)
                profile = ediff.radial.calc_radial_distribution(arr)
            else:  # ......................................TXT file, 1D-profile
                profile = ediff.radial.read_radial_distribution(data)
        elif type(data) == np.ndarray:  # Numpy array
            if data.shape[0] == data.shape[1]:  # sqare array, 2D-diffractogram
                profile = ediff.radial.calc_radial_distribution(data)
            else:  # ..............................non-square rrray, 1D-profile
                profile = data
        # Read plot parameters        
        my_format = rdist[i][1]
        my_label  = rdist[i][2]
        # Plot data
        R,I = profile
        plt.plot(R,I, my_format, label=my_label)
    # ...adjust plot
    plt.xlabel('Radial distance [pixel]')
    plt.ylabel('Intensity [grayscale]')
    plt.xlim(0,xlimit)
    plt.ylim(0,ylimit)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    # ...save plot as PNG (only if argument [output] was given)
    if output_file: plt.savefig(output_file, dpi=300)
    # ...show plot
    plt.show()
