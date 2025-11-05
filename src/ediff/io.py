'''
Module: ediff.io
----------------
Input/output functions for package EDIFF.    
'''

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import ediff.radial

from pymatgen.core import Lattice as pmLattice
from pymatgen.core import Structure as pmStructure


class Lattice(pmLattice):
    '''
    Lattice object = crystal lattice.

    * Lattice object is identical to pymatgen.core.Lattice: <br>
      https://pymatgen.org/pymatgen.core.html#module-pymatgen.core.lattice
    * In EDIFF:
        - Lattice can be defined by all methods of the original object.
        - Lattice can help to define a crystal structure = ediff.io.Structure.
        
    >>> # How to define Lattice = crystall lattice using ediff.io
    >>> import ediff as ed
    >>> lattice1 = ed.io.Lattice.cubic(a=5.47)
    >>> lattice2 = ed.io.Lattice.hexagonal(a=5.91, c=3.50)
    >>> lattice3 = ed.io.Lattice.from_parameters(3, 4, 5, 90, 110, 90) 
    '''
    pass


class Structure(pmStructure):
    '''
    Structure object = crystal structure.

    * Structure object is identical to pymatgen.core.Structure: <br>
      https://pymatgen.org/pymatgen.core.html#module-pymatgen.core.structure
    * In EDIFF:
        - Structure can be defined by all methods of the original object.
        - Structure is employed in calculation of theoretical diffractograms.
    
    >>> # How to define Structure = crystal structure using ediff.io
    >>>
    >>> # (0) Standard import of ediff
    >>> import ediff as ed
    >>>
    >>> # (1) Structure from the very beginning
    >>> sg = 'Fm-3m'
    >>> lat = ed.io.Lattice.cubic(a=4.08) 
    >>> atoms = ['Au']
    >>> coords = [[0,0,0]]
    >>> structure1 = ed.io.Structure.from_spacegroup(sg, lat, atoms, coords)
    >>>
    >>> # (2) Structure from CIF-file
    >>> struct2 = ed.io.Structure.from_file('au.cif')
    '''
    pass


class Diffractogram:
    '''
    Read and show 2D diffractograms in a simple and reproducible way.
    
    * A class with two funcs (read, show); show func can be used for saving.
    * Assumption: the diffractogram is a 2D numpy array or grayscale image.
    '''
    
    
    def read(diffractogram, itype=None):
        '''
        Read 2D diffraction pattern (=grayscale image) into 2D numpy array.
        
        Parameters
        ----------
        diffractogram : string or path-like object or numpy array
            Name of image that should read into numpy 2D array
            or directly the 2D numpy array representing the diffractogram.
        itype : string ('8bit'  or '16bit')
            Type of the image: 8 or 16 bit grayscale    
            
        Returns
        -------
        arr : 2D numpy array
            The *arr* is the input image read to an array by means of numpy.
        '''
        
        # Test the input type
        if type(diffractogram) == np.ndarray:
            # Diffractogram is an array - just assign to arr variable.
            arr = diffractogram
        else:
            # Diffractogram is not an array - we expect a grayscale image.
            img = Image.open(diffractogram)
            if itype=='8bit':
                arr = np.asarray(img, dtype=np.uint8)
            else:
                arr = np.asarray(img, dtype=np.uint16)
        # Return the diffractogram in the form of 2D array.
        return(arr)


    def show(diffractogram, icut=None, origin=None, 
             title=None, output_file=None, output_file_dpi=300):
        '''
        Show/plot 2D diffraction pattern.

        Parameters
        ----------
        diffractogram : numpy.array
            A numpy.array object representing a 2D diffractogram image.
            In EDIFF,
            this array is usually obtained by ediff.ioi.read_image function.
        icut : integer, optional, default is None
            Upper limit of intensity shown in the diffractogram.
            The argument *icut* is used as *vmax* in plt.imshow function.
            Example: If *icut*=300, then all intensities >300 are set to 300.
        origin : 'upper' or 'lower' or None, optional, default is None
            Orientation of the image during final rendering.
            If the argument is None, we follow the Matplotlib default,
            which is *origin*='upper' = [0,0] in the upper left corner.
            Alternative: *origin*='lower' = [0,0] is in the lower left corner.
        title : str, optional, default is None
            If given, then it is the title of the plot.
        output_file : str, optional, default is None
            Name of the output file.
            If the argument is not None, the plot is saved to *output_file*.
        output_file_dpi : int, optional, default is 300
            Resolution of the output file.

        Returns
        -------
        None
            The plot is shown in the stdout
            and saved to *output_file* if requested.
        '''
        
        # Plot the 2D-diffractogram
        # (quite simple, we employ plt.imshow function with a few arguments
        # (the function is defined in order to simplify user's input even more
        
        # (1) Plot title if requested
        if title is not None: plt.title(title)
        
        # (2) The plot itself
        plt.imshow(diffractogram, origin=origin, vmax=icut)
        plt.tight_layout()
        
        # (3) Save the plot if requested
        if output_file is not None:
            plt.savefig(output_file, dpi=output_file_dpi)
        
        # (4) Show the plot
        plt.show()

    
class Profile:
    '''
    Read and show 1D diffraction profiles in a simple and reproducible way.
    
    * A class with two funcs (read, show); show func can be used for saving.
    * Assumption: the profile is a numpy array or text file in EDIFF format. 
    '''
    
    def read(profile):
        '''
        Read the ELD or XRD profile in EDIFF format.
    
        * More info about ELD/XRD profiles in EDIFF
          => see the section *Technical notes* below.
    
        Parameters
        ----------
        profile : str or numpy.array
            
            * If profile = str,
              we assume a filename
              of the file with ELD or XRD profile in EDIFF format.
            * If profile = numpy.array,
              we assume a 2D-array
              containing ELD or XRD profile in EDIFF format.
    
        Returns
        -------
        profile : 2D numpy.array
            The array representing ELD or XRD profile in EDIFF format.
            See section *Technical notes* below
            for explanation of the EDIFF format of the ELD and XRD profiles.
        
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
        if type(profile)==np.ndarray:
            return(profile)
        else:
            profile = np.loadtxt(profile, unpack=True)
            return(profile)

    def show(Xvalues, Yvalues, Xlabel, Ylabel, Xrange, Yrange,
        title=None, output_file=None, output_file_dpi=300):
        '''
        Plot a 1D profile in a simple and stadnard way.

        Parameters
        ----------
        Xvalues : array or list-like object
            X values for plotting.
        Yvalues : array or list-like object
            Y values for plotting.
        Xlabel : str
            Label of the X-axis.
        Ylabel : str
            Label of the Y-axis.
        Xrange : list/tuple of two floats
            X range = minimum and maximu for Xvalues to plot.
        Yrange : list/tuple of two floats
            Y range = minimum and maximu for Yvalues to plot.
        title : str, optional, default is None
            The title of the plot.
        output_file : str, optional, default is None
            Name of the output file.
            If the argument is not None, the plot is saved to *output_file*.
        output_file_dpi : int, optional, default is 300
            Resolution of the output file.

        Returns
        -------
        None
            The plot is shown in the stdout
            and saved to *output_file* if requested.
       '''
        
        # (1) Plot title if requested
        if title is not None: plt.title(title)
        
        # (2) The plot itself
        plt.plot(Xvalues, Yvalues)
        plt.xlabel(Xlabel)
        plt.ylabel(Ylabel)
        plt.xlim(Xrange)
        plt.ylim(Yrange)
        plt.grid()
        plt.tight_layout()
        
        # (3) Save the plot if requested
        if output_file is not None:
            plt.savefig(output_file, dpi=output_file_dpi)
        
        # (4) Show the plot
        plt.show()


def set_plot_parameters(
        size=(8,6), dpi=100, fontsize=8, my_defaults=True, my_rcParams=None):
    '''
    Set global plot parameters (mostly for plotting in Jupyter).

    Parameters
    ----------
    size : tuple of two floats, optional, the default is (8,6)
        Size of the figure (width, height) in [cm].
    dpi : int, optional, the defalut is 100
        DPI of the figure.
    fontsize : int, optional, the default is 8
        Size of the font used in figure labels etc.
    my_defaults : bool, optional, default is True
        If True, some reasonable additional defaults are set,
        namely line widths and formats.
    my_rcParams : dict, optional, default is None
        Dictionary in plt.rcParams format
        containing any other allowed global plot parameters.

    Returns
    -------
    None
        The result is a modification of the global plt.rcParams variable.
    '''
    # (1) Basic arguments -----------------------------------------------------
    if size:  # Figure size
        # Convert size in [cm] to required size in [inch]
        size = (size[0]/2.54, size[1]/2.54)
        plt.rcParams.update({'figure.figsize' : size})
    if dpi:  # Figure dpi
        plt.rcParams.update({'figure.dpi' : dpi})
    if fontsize:  # Global font size
        plt.rcParams.update({'font.size' : fontsize})
    # (2) Additional default parameters ---------------------------------------
    if my_defaults:  # Default rcParams if not forbidden by my_defaults=False
        plt.rcParams.update({
            'lines.linewidth'    : 0.8,
            'axes.linewidth'     : 0.6,
            'xtick.major.width'  : 0.6,
            'ytick.major.width'  : 0.6,
            'grid.linewidth'     : 0.6,
            'grid.linestyle'     : ':'})
    # (3) Further user-defined parameter in rcParams format -------------------
    if my_rcParams:  # Other possible rcParams in the form of dictionary
        plt.rcParams.update(my_rcParams)

    
def plot_final_eld_and_xrd(eld_profile, xrd_profile, fine_tuning, x_range,
        eld_data_label='ED experiment', xrd_data_label='XRD calculation',
        x_axis_label='$q$ [1/\u212B]', y_axis_label='Intensity',
        xticks=None, yticks=None, mxticks=None, myticks=None,
        output_file=None, output_file_dpi=300, transparent=False, CLI=False):
    '''
    Final plot/comparison of ELD and XRD profiles.

    * During the final plotting, we fine-tune the ELD calibration.
    * This is done by iterative modification of fine_tuning constant.
    
    Parameters
    ----------
    eld_profile : str or numpy.array
        The *eld_profile* (ELD) is
        an electron diffraction profile in EDIFF format.
        It can come as file (if *eld_profile* = str = filename)
        or array (if *eld_profile* = numpy.array).
        More info about ELD profiles in EDIFF
        => see docs of ediff.calibration.Calculate.from_max_peaks function.
    xrd_profile : str or numpy.array
        The *xrd_profile* (XRD) is
        an X-rayd diffraction profile in EDIFF format.
        It can come as file (if *xrd_profile* = str = filename)
        or array (if *xrd_profile* = numpy.array).
        More info about XRD profiles in EDIFF
        => see docs of ediff.calibration.Calculate.from_max_peaks function.
    fine_tuning : float
        The constant for the final fine-tuning of peak position.
        The *fine_tuning* constant has a starting value of 1.000.
        If ELD and XRD peaks are shifted, the constant should be adjusted.
        The constant multiplies the X-values of ELD profile.
    x_range : tuple of two floats
        The limits for X-axis (minimum and maximu q-vectors on X-axis).
    eld_data_label : str, optional, default is 'ED experiment'
        The label of ELD data (= name of the electron diffraction data).
    xrd_data_label : str, optional, the default is 'XRD calculation'
        The label of XRD data (= name of the X-ray diffraction data).
    x_axis_label : str, optional, the default is '$q$ [1/\u212B]' ~ q [1/A]
        The label of X-axis.
    y_axis_label : str, optiona, the default is 'Intensity'.
        The label of Y-axis.
    xticks : float, optional, default is None
        The X-axis ticks (if not omitted, use the default).
    yticks : float, optional, default is None
        The Y-axis ticks (if not omitted, use the default).
    mxticks : float, optional, default is None
        The Y-axis minor ticks (if not omitted, use the default).
    myticks : float, optional, default is None
        The Y-axis minor ticks (if not omitted, use the default).
    output_file : str, optional, default is None
        The filename, to which the final graph should be saved.
        If *output_file* is not None (= the default),
        the plot is not only shown in stdout,
        but also saved in the *output_file*.
    output_file_dpi : int, optional, default is 300
        The DPI of the output graph.
    transparent : bool, optional, default is False
        If *transparent* = True, then the image background is transparent.
    CLI : bool, optional, default is False
        If *CLI* = True, we assume command line interface
        and the plot is not shown, just saved.

    Returns
    -------
    None
        The plot is shown in the stdout
        and saved to *output_file* if requested.
    '''
    
    # (1) Read ED and XRD diffraction profiles.
    # * The profiles are supposed to be either filenames or numpy.arrays
    # * In any case, the filenames or arrays should be in EDIFF format:
    #   ED  = 3 columns: pixels, intensity, bkgr-corrected intensity
    #   XRD = 4 columns: 2theta[deg], S[1/A], q[1/A], normalized-intensity
    eld = Profile.read(eld_profile)
    xrd = Profile.read(xrd_profile)

    # Plot the data
    plt.plot(xrd[2], xrd[3], label=xrd_data_label)
    plt.plot(eld[0]*fine_tuning, eld[2], color='red', label=eld_data_label)
    # Define axis labels
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)
    # Define xlim = x-limits = x-range
    plt.xlim(x_range)
    # Define ticks and minor ticks if requested
    if xticks is not None:
        plt.gca().xaxis.set_major_locator(plt.MultipleLocator(xticks))
    if yticks is not None:
        plt.gca().yaxis.set_major_locator(plt.MultipleLocator(yticks))
    if mxticks is not None:
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(mxticks))
    if myticks is not None:
        plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(myticks))
    # Add legent (considering transparency)
    if transparent == True:
        plt.legend(framealpha=0.0)
    else: plt.legend()
    # Additional parameters
    plt.grid()
    plt.tight_layout()
    # Save plot if requested
    if output_file is not None:
        if transparent == True:
            plt.savefig(output_file, dpi=output_file_dpi, transparent=True)
        else:
            plt.savefig(output_file, dpi=output_file_dpi, facecolor='white')
    # Show the plot
    if CLI == False: plt.show()


def plot_radial_distributions(
        data_to_plot, xlimit, ylimit, output_file=None):
    """
    Plot one or more 1D-radial distrubution files in one graph.
    
    * This is a rather specific function.
    * It is employed mostly when we combine STEMDIFF and EDIFF.

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
    None
        The plot is shown in the stdout
        and saved to *output_file* if requested.
    
    Technical notes
    ---------------
    * This function is quite flexible.
    * It can plot one radial distribution or more.
    * It can take data from PNG-files, TXT-files, 2D-arrays and 1D-arrays.
    * If the input is a PNG-file or2D-array,
      the center is just *estimated* as as the center of intensity;
      therefore, this works only for good diffractograms with a central spot.
    * This makes the code a more complex, but it is convenient for the user.
    * An example of fast comparison of three 1D-distributions
      taken from three 2D-diffractograms in the form of 16-bit PNG images:
    
    >>> ediff.io.plot_radial_distributions(
    >>>     data_to_plot = [
    >>>         ['sum_all_16bit.png', 'k:',  'All data'],
    >>>         ['sum_f_16bit.png',   'b--', 'F data'],
    >>>         ['sum_fd_16bit.png',  'r-',  'FD data']]
    >>>     xlimit=200, ylimit=300,
    >>>     output_file='sums_final_1d.png')
    """
    # Initialize
    n = len(data_to_plot)
    rdist = data_to_plot
    # Plot radial distribution files
    for i in range(n):
        # Read data
        data = rdist[i][0]
        if type(data) == str:  # Datafile
            if data.lower().endswith('.png'):  # ....PNG file, 2D-diffractogram
                arr = Diffractogram.read(data)
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
    # ...change Jupyter default transparent color to white
    plt.gcf().patch.set_facecolor('white')
    # ...save plot as PNG (only if argument [output] was given)
    if output_file: plt.savefig(output_file, dpi=300)
    # ...show plot
    plt.show()
