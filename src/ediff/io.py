'''
Module: ediff.io
----------------
Input/output functions for package EDIFF.    
'''

import numpy as np
import pandas as pd
import skimage as ski
import matplotlib.pyplot as plt
from pathlib import Path

import ediff.radial

from pymatgen.core import Lattice as pmLattice
from pymatgen.core import Structure as pmStructure


class Lattice(pmLattice):
    '''
    Lattice object = crystal lattice.

    * Lattice object is identical to pymatgen.core.Lattice:
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

    * Structure object is identical to pymatgen.core.Structure:
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


class Data1D:
    '''
    Read, show and save 1D-data (such as 1D diffraction profiles).
    
    * PP (procedural programming) class - intentionally.
    * It is kept for simple data treatment and backward compatibility.
    * Note #1: the profiles are in pandas-saved files (with column headers).
    * Note #2: the show function can be used for saving (using out_file arg).
    
    TODO :: this function is finished, but not fully incorporated.
        
    * In the next step, the whole package/examples should be adjusted
      to read/save XRD and ELD profiles as pandas.DataFrames reproducibly.
    '''

    
    def read(XYdata, **kwargs):
        '''
        Read XY data as dataframe with several columns.
    
        Parameters
        ----------
        XYdata : str or pandas.DataFrame or numpy.ndarray
            
            * If XYdata = str,
              we assume that it is a filename
              of the file with several columns of data.
              The optional argument {dtype} determines how to read the data.
              The optional argument {dcols} determines how to name the columns.
            * If XYdata = pandas.Dataframe,
              we simply read the dataframe and return it.
              Optionally, we can change column names if {dcols} is given.
            * If XYdata = numpy.ndarray,
              we convert the array to DataFrame
              and define column names if {dcols} is given.
              
    
        Returns
        -------
        XYdata : pandas.Dataframe
            The DataFrame representing XYdata.
            Typically, the XYdata represent a XRD or ELD diffraction profile.        
        '''
        
        # (1) Check if the column names were defined.
        #  - the column names can be defined/changed when calling this func
        #  - typical column names for XRD profiles: [TwoTheta, S, q, Intensity]
        #  - typical column names for ELD profiles: [q, Iraw, Ibkg, I] 
        my_columns  = kwargs.get('columns') or None
        
        # (2) Read XYdata depending on the data type.
        #  - the data type can be filename/pd.DataFrame/np.ndarray
        if isinstance(XYdata, str):
            # XYdata is a str = file
            # ... comments may be present
            my_comment  = kwargs.get('comment') or '#'
            # ... header should be None if there are no column headers!
            my_header   = kwargs.get('header') or None
            # ... skiprows can be defined for completness and flexibility
            my_skiprows = kwargs.get('skiprows') or 0
            # ... separator can be defined, defalt is any whitespace
            my_sep      = kwargs.get('sep') or r'\s+'
            # ... usecols may be defined, default is the first two columns
            my_usecols  = kwargs.get('usecols') or None
            # (c) read the file with pd.read_csv and the **kwargs from above
            # (Trick: we are not sure, if the file has header
            # (Step1: read with header=None
            # (Step2: try to convert df to float
            # (Step3: if it fails, 1st line was header => re-read with header=0
            # (Note: this works also with comment='#' => 1st non-commented line
            # -----
            # Step1: reading header given as argument OR header=None (default)
            df = pd.read_csv(XYdata, 
                comment=my_comment, header=my_header, skiprows=my_skiprows, 
                sep=my_sep, usecols=my_usecols)
            # Step2: try if the first line is header or not
            try:
                df.astype(float)
            except ValueError:
                # Step3: if conversion failed, re-read with header=0
                df = pd.read_csv(XYdata, 
                    comment=my_comment, header=0, skiprows=my_skiprows, 
                    sep=my_sep, usecols=my_usecols)
            # (d) Set column names if columns was given as argument
            if my_columns is not None: df.columns = my_columns
            return(df)
        elif isinstance(XYdata, pd.DataFrame):
            # XYdata is a DataFrame => add column names if required and return.
            if my_columns is not None: XYdata.columns = my_columns
            return(XYdata)
        elif isinstance(XYdata, np.ndarray):
            # XYdata is an array => convert to Dataframe.
            df = pd.DataFrame(np.transpose(XYdata))
            if my_columns is not None: df.columns =  my_columns
            return(df)
        elif XYdata is None:
            # XYdata is None => create an empty dataframe with {dcols}
            # * XYdata = None is a rare, but real case
            #   It may happen if we call: ELD = ediff.io.Profile(),
            #   which calls ediff.io.Data1D.read internally (and XYdata=None).
            # * The following command simply creates
            #   empty DataFrame with predefined column names,
            #   and it works even without column names => for {dcols} = None.
            df = pd.DataFrame(columns=my_columns)
            return(df)
        else:
            raise TypeError('Unknown format of XYdata!')


    def show(df, X, Y, ptype='plot',
             Xlabel=None, Ylabel=None, Xlim=None, Ylim=None,
             title=None, out_file=None, out_dpi=300,
             **kwargs):
        '''
        Show/plot a 1D profile in a simple and stadnard way.
        
        Parameters
        ----------
        df : pandas.Dataframe
            A DataFrame, from which we select columns for plotting.
        X : str
            X values for plotting, column name in {df}.
        Y : str
            Y values for plotting, column name in {df}.
        ptype : 'plot' or 'vlines', optional, default is 'plot'
            Plot type. The values 'plot' and 'vlines' are suitable for
            the diffraction profile and individual diffractions, respectively.
        Xlabel : str, optional, default is None
            Label of the X-axis.
        Ylabel : str, optional, default is None
            Label of the Y-axis.
        Xlim : float ar list/tuple of two floats, optional, default is None
            X-range, minimum and maximum of X-axis.    
            If Xlim = 300, then we plot with plt.xlim(0,300).
            If Xlim = (50,350), the we plot with plt.xlim(50,350).
        Ylim : float or list/tuple of two floats, optional, default is None
            Y-range, minimum and maximum of Y-axis.    
            If Ylim = 300, then we plot with plt.ylim(0,300).
            If Ylim = (50,350), the we plot with plt.ylim(50,350).
        title : str, optional, default is None
            The title of the plot.
        out_file : str, optional, default is None
            Name of the output file.
            If the argument is not None, the plot is saved to *output_file*.
        out_dpi : int, optional, default is 300
            Resolution of the output file.
        **kwargs : arguments supplied to plt.plot command
            Any arguments that are transferred to plt.plot command.
            For example: color='red', linestyle=':', linewidth='1', label='ED'
            
        Returns
        -------
        None
            The plot is shown in the stdout
            and saved to {out_file} if the arqument is given.
        '''
        
        # (0) Check if the X,Y are existing column names.
        if X not in df.columns:
            raise ValueError(f"Column '{X}' does not exist!")
        if Y not in df.columns:
            raise ValueError(f"Column '{Y}' does not exist!")
        
        # (1) Basic plot
        if ptype == 'plot':
            plt.plot(df[X], df[Y], **kwargs)
        elif ptype == 'vlines':
            plt.vlines(df[X], 0, df[Y], **kwargs)
        else:
            raise TypeError('Unknown plot type!')
            
        # (2) Optional arguments
        # Plot title :: if requested
        if title is not None: plt.title(title)
        # Axes labels :: if Xlabel,Ylabel not given, use column names instead
        if Xlabel is None: Xlabel = X
        if Ylabel is None: Ylabel = Y
        plt.xlabel(Xlabel)
        plt.ylabel(Ylabel)
        # Axes limits :: if not defined, use defaults
        if Xlim is not None: 
            if isinstance(Xlim, (int,float)): 
                plt.xlim(0,Xlim)
            elif isinstance(Xlim, (tuple,list)): 
                plt.xlim(Xlim)
            else:
                raise TypeError('Wrong type of Xlim argument!')
        if Ylim is not None: 
            if isinstance(Ylim, (int,float)): 
                plt.ylim(0,Ylim)
            elif isinstance(Ylim, (tuple,list)): 
                plt.ylim(Ylim)
            else:
                raise TypeError('Wrong type of Ylim argument!')
        
        # (3) Additional plot parameters
        plt.grid()
        plt.tight_layout()
        
        # (4) Save the plot if requested
        if out_file is not None: plt.savefig(out_file, dpi=out_dpi)
        
        # (5) Show the plot
        plt.show()
        
    
    def save(df, out_file, pretty=True, **kwargs):
        '''
        Save the {df} containing XYdata to TXT-file in EDIFF format.
        
        * EDIFF format: DataFrame saved as TXT, the 1st line = column names.
        
        Parameters
        ----------
        df : pandas.DataFrame
            XYdata = several (named) columns of data
        out_file : str or PathLike object
            Name of the output file (filename or complete path).
        pretty : bool, optional, default is True
            If {pretty} is True (default)
            use df.to_string to get nicely fomatted columns in the {out_file}.
        **kwargs : optional arguments for pandas
            The {**kwargs} are passed to df.to_string or df.to_csv methods.
            Typically, it can be something optimizing the output format,
            such as: `float_format='%.2f'
            or `formatters=({'pixels':'{:d}'.format,'I0:'{:.2f}'.format})`

        Returns
        -------
        None
            The DataFrame containing XYdata is just saved to file.
            The file contains column names and the data.
            Later it can be read by ediff.io.Data1D.read function.
        '''
        
        # Convert out_file to Path object if it is just a string/filename.
        if isinstance(out_file, str):
            out_file = Path(out_file)
        
        # Save df to out_file.
        # (**kwargs can be anything optimizing the output format
        # (Ex 1: foat_format='%.2f'
        # (Ex 2: formatters=({'pixels':':d'.format,'I0':':.2f'.format})
        if pretty is True:
            # If {pretty} = True => save df with the cols nicely aligned.
            out_file.write_text(df.to_string(index=False, **kwargs))
        else:
            # If {pretty} != True => save df as CSV with tab as separator.
            df.to_csv(out_file, index=False, sep='\t', **kwargs)
            

class Data2D:
    '''
    Read, show and save 2D-data (such as 2D diffraction patterns).
    
    * PP (procedural programming) class - intentionally.
    * It is kept for simple data treatment and backward compatibility.
    * Note 1: the diffractogram is a 2D numpy array or grayscale image.
    * Note 2: only '8bit' and '16bit' images/arrays are allowed.
    * Note 3: {show} function can save the plot - using {out_file} argument.
    * Note 4: {save} function can save the array - as an 8bit/16bit image.
    '''
       
    def read(image):
        '''
        Read 2D diffraction pattern (image or array) into 2D numpy array.
        
        Parameters
        ----------
        image : str or PathLike object or numpy.ndarray
            Name of image that should read into numpy 2D array
            or directly the 2D numpy array representing the diffractogram.
          
        Returns
        -------
        arr : 2D numpy array
            The *arr* is the input image read to an array by means of numpy.
        '''
        
        # Test the input type
        if isinstance(image, np.ndarray):
            # Diffractogram is an array - just assign to arr variable.
            arr = image
        else:
            # Diffractogram is not an array - we expect a grayscale image.
            # (a) Read the image
            # (skimage should read both 8bit and 16bit, with correct img.dtype
            img = ski.io.imread(image)
            # (b) Grayscale images can be 1-channel (G) or 3-channel (R=G=B)
            # (the following code converts the image to 1-channel grayscale
            img = img[..., 0] if img.ndim == 3 else img
            # (c) Just to be sure, check if the image is 8bit or 16bit
            # (if it is 8bit in 16bit containers, we leave it like this
            # (ImageJ sees it as 8bit even if dtype is uint16
            if img.dtype not in (np.uint8, np.uint16):
                raise TypeError("Only 8bit and 16bit images allowed!")
            # (d) Assign img to arr if everything was Ok so far
            arr = img
        
        # Return the diffractogram in the form of 2D array.
        return(arr)


    def show(image, title=None, icut=None, cmap='gray',
             origin=None, out_file=None, out_dpi=300):
        '''
        Show/plot 2D data (usually a 2D diffraction pattern).

        Parameters
        ----------
        image : numpy.array
            A numpy.array object corresponding to an image.
            Typically, the array represents a 2D diffraction pattern,
            which had been read by means by means of ediff.io funcs.
        title : str, optional, default is None
            If given, then it is the title of the plot.
        icut : integer, optional, default is None
            Upper limit of intensity shown in the diffractogram.
            The argument {icut} is used as *vmax* in plt.imshow function.
            Example: If {icut}=300, then all intensities >300 are set to 300.
        origin : 'upper' or 'lower' or None, optional, default is None
            Orientation of the image during final rendering.
            If the argument is None, we follow the Matplotlib default,
            which is {origin}='upper' = [0,0] in the upper left corner.
            Alternative: {origin}='lower' = [0,0] is in the lower left corner.
        out_file : str, optional, default is None
            Name of the output file.
            If the argument is not None, the plot is saved to {output_file}.
        out_dpi : int, optional, default is 300
            Resolution of the output file.

        Returns
        -------
        None
            The plot is shown in the stdout
            and saved to {out_file} if requested.
        '''
        
        # Plot the 2D-diffractogram
        # (quite simple, we employ plt.imshow function with a few arguments
        # (the function is defined in order to simplify user's input even more
        
        # (1) Plot title if requested
        if title is not None: plt.title(title)
        
        # (2) The plot itself
        plt.imshow(image, cmap=cmap, origin=origin, vmax=icut)
        plt.tight_layout()
        
        # (3) Save the plot if requested
        if out_file is not None:
            plt.savefig(out_file, dpi=out_dpi)
        
        # (4) Show the plot
        plt.show()
    
        
    def save(image, out_file, itype=None, icut=None):
        '''
        Save 2D data (usually a 2D diffraction pattern) as an image file.

        Parameters
        ----------
        image : numpy.array
            A numpy.array object corresponding to an image.
            Typically, the array represents a 2D diffraction pattern,
            which had been read by means by means of ediff.io funcs.
        out_file : str or PathLike object
            The name of the file into which the image/array will be saved.
            The image format is deduced from the extension of the {out_file}.
        itype : '8bit or '16bit', optional, default is '16bit'
            Bit depth of the saved image (8bit or 16bit grayscale).
            The bit depth is auto-estimated by skimage.io.save
            from {image}.dtype but it can be also set here explicitly.
            The explicit setting is maintained mostly for historical reasons.
        icut : int or float, optional, default is None
            Intensity cut.
            If {icut} = 200, then all intensity values >200 are set to 200.

        Returns
        -------
        None
            The input image/array is saved in {out_file}.
        '''
        # Apply {icut} if requested.
        if icut is not None:
            image = np.where(image > icut, icut, image)
        
        # Check the image type and hard-convert to required type if needed
        if itype == '8bit':
            if image.dtype != np.uint8: image = ski.img_as_ubyte(image)
        elif itype == '16bit':
            if image.dtype != np.uint16: image = ski.img_as_uint(image)
        else:
            raise ValueError("Only '8bit' or '16bit' images are allowed!")
        
        # Save the image.
        ski.io.imsave(out_file, image)


class DiffData1D(pd.DataFrame):
    '''
    Class defining a general 1D diffraction data object.
    
    * The general object is a pandasDataFrame with additional methods.
    * The object is used to defie two specific objects = subclasses:
      ediff.io.Diffractions1D or ediff.io.Diffractogram1D.
    * The object has three additional methods (has, show, save);
      the show method is adjusted in the subclasses.
    '''
    
    # (1) DataFrame subclassing, item 1 :: _metadata
    # Optional class property, useful to keep properties with selected names.
    # * pandas convention used only by DataFrame/Series objects
    # * it keeps selected properties in derived objects
    #   >>> profile.name = 'Au powder'
    #   >>> profile2 = profile[profile.q > 1]
    #   >>> profile2.name  # works only if _metadata=['name'], otherwise error
    
    _metadata = ['name']
    
    # (2) DataFrame subclassing, item 2 :: _constructor
    # Necessary class property to keep our object after any pandas operation.
    # * pandas operation return a plan DataFrame (if _constructor not defined)
    #   >>> type(profile)        # Profile
    #   >>> type(profile.copy()  # pandas.DataFrame => wrong!
    
    @property
    def _constructor(self):
        return DiffData1D
    
    # (3) Dataframe subclassing, item 3 :: __init__
    # Necessary for correct object initialization.
    
    def __init__(self, data=None, name=None, **kwargs):
        
        # (a) EDIFF-specific pre-processing:
        # * In a general subclass of pandas.DataFrame, this would be skipped
        # * In this specific subclass, we pre-read the profile with Data1D.read
        # * Reason: the file can come from numpy and/or come without col names
        # * Solution: we pre-read data to {df} and pass {df} to super.__init__
        df = Data1D.read(data, **kwargs)        
        
        # (b) STANDARD SUBCLASSING
        # (the columns are NOT fixed/pre-defined intentionally
        # (the later/lazy initialization of arbitrary columns is more flexible
        super().__init__(df)
        # Note: if it were not for the previous EDIFF pre-processing,
        #  the correct command would be: super().__init__(df, *args, **kwargs)
        
        # (c) Object name
        # (optional, it can be something like: 'Au powder'
        # (we define it only if {name} argument is not None
        # (reason: if we auto-construct object somehow, old name can be kept
        # (see also the comments above, around the {_metadata} class variable
        if name is not None: self.name = name
    
    # (4) DataFrame subclassing, item 4 :: additonal functions
    
    def has(self, column):
        '''
        Test if the {Profile} object contains specific column.

        Parameters
        ----------
        column : str
            Column name.

        Returns
        -------
        Bool
            True if the {Profile} object contains the {column},
            false atherwise.
        '''
        return column in self.columns

    
    def show(self, X, Y, 
             Xlabel=None, Ylabel=None, Xlim=None, Ylim=None,
             title=None, out_file=None, out_dpi=300,
             **kwargs):
        '''
        Show/plot Profile in a simple and stadnard way.

        * This function is a wrapper of ediff.io.Data1D.show.
        
        Parameters
        ----------
        X : str
            X values for plotting, column name in {df}.
        Y : array or list-like object
            Y values for plotting.
        Xlabel : str, optional, default is None
            Label of the X-axis.
        Ylabel : str, optional, default is None
            Label of the Y-axis.
        Xlim : list/tuple of two floats, optional, default is None
            X range = minimum and maximu for Xvalues to plot.
        Ylim : list/tuple of two floats, optional, default is None
            Y range = minimum and maximu for Yvalues to plot.
        title : str, optional, default is None
            The title of the plot.
        out_file : str, optional, default is None
            Name of the output file.
            If the argument is not None, the plot is saved to *output_file*.
        out_dpi : int, optional, default is 300
            Resolution of the output file.
        **kwargs : arguments supplied to plt.plot command
            Any arguments that are transferred to plt.plot command.
            For example: `color='red', linestyle=':', linewidth='1'`
            
        Returns
        -------
        None
            The plot is shown in the stdout
            and saved to {out_file} if requested.
        '''
        Data1D.show(self, X, Y,
                    Xlabel, Ylabel, Xlim, Ylim,
                    title, out_file, out_dpi, **kwargs)
    
    
    def save(self, out_file, pretty=True, fine_tune=None, **kwargs):
        '''
        Save {Profile} object to TXT-file in EDIFF format.

        * This method is a wrapper of ediff.io.Data1D.save.
        * EDIFF format: DataFrame saved as TXT, the 1st line = column names.
        
        Parameters
        ----------
        out_file : str or PathLike object
            Name of the output file (filename or complete path).
        pretty : bool, optional, default is True
            If {pretty} is True (default)
            use df.to_string to get nicely fomatted columns in the {out_file}.
        kwargs : optional arguments for pandas
            The {kwargs} are passed to df.to_string or df.to_csv methods.
            Typically, it can be something optimizing the output format,
            such as: `float_format='%.2f'
            or `formatters=({'pixels':'{:d}'.format,'I0:'{:.2f}'.format})`

        Returns
        -------
        None
            The DataFrame containing XYdata is just saved to file.
            The file contains column names and the data.
            Later it can be read by ediff.io.Data1D.read function.
        '''
        
        # Special case: multiplicative fine_tune coefficient for q-vectors
        if (fine_tune is not None) and ('q' in self.columns):
            diffractogram = self.copy()
            diffractogram.q = diffractogram.q * fine_tune
            Data1D.save(diffractogram, out_file, pretty, **kwargs)
        else:
            Data1D.save(self, out_file, pretty, **kwargs)


class Diffractions1D(DiffData1D):
    
    _metadata = ['name']
    
    @property
    def _constructor(self):
        return Diffractogram1D
    
    def show(self, X, Y, 
             Xlabel=None, Ylabel=None, Xlim=None, Ylim=None,
             title=None, out_file=None, out_dpi=300,
             **kwargs):
        
        Data1D.show(self, X, Y, 'vlines',
                    Xlabel, Ylabel, Xlim, Ylim,
                    title, out_file, out_dpi, **kwargs)


class Diffractogram1D(DiffData1D):
    
    _metadata = ['name']
    
    @property
    def _constructor(self):
        return Diffractogram1D
    
    def show(self, X, Y, 
             Xlabel=None, Ylabel=None, Xlim=None, Ylim=None,
             title=None, out_file=None, out_dpi=300,
             **kwargs):
        
        Data1D.show(self, X, Y, 'plot',
                    Xlabel, Ylabel, Xlim, Ylim,
                    title, out_file, out_dpi, **kwargs)
        

class Diffractogram2D(np.ndarray):
    '''
    Read, show, work with, and save 2D diffractograms.
    '''
    
    def __new__(cls, image):
        '''
        Initialize the {Diffractogram} object.
        
        * {Diffractogram} object is a numpy.ndarray with additional functions.
        * If we subclass numpy.ndarray objects, we use __new__ method.
        * Reason: specific numpy feature, recommended way.

        Parameters
        ----------
        image : str or array or PIL-image
            The image representing the 2D diffraction pattern.

        Returns
        -------
        Diffractogran object
            Diffractogram object is a numpy.array with additional methods.
        '''
                
        # Create the object = np.ndarray.
        # * The function np.asarray is surprisingly versatile!
        # * It can transform many common objects to np.ndarray, such as:
        #   arrays, list-of-lists, PIL images, and other array-like objects...
        # * The pattern: obj = np.asarray(input).view(cls)
        #   is the official NumPy pattern how to create objects from arrays.
        
        # (1) Read array or image
        arr = ediff.io.Data2D.read(image)
        
        # (2) Use the above-mentioned numpy initialization trick
        # (obj = np.asarray(input).view(cls)
        obj = np.asarray(arr).view(cls)
        
        # (3) Return the created object/array.
        return obj

    
    def show(self, title=None, icut=None, cmap='gray',
             origin=None, out_file=None, out_dpi=300):
        '''
        Show/plot the {Diffractogram} object.
        
        * This function is a wrapper of ediff.io.Data2D.show.
        
        Parameters
        ----------
        title : str, optional, default is None
            If given, then it is the title of the plot.
        icut : integer, optional, default is None
            Upper limit of intensity shown in the diffractogram.
            The argument {icut} is used as *vmax* in plt.imshow function.
            Example: If {icut}=300, then all intensities >300 are set to 300.
        origin : 'upper' or 'lower' or None, optional, default is None
            Orientation of the image during final rendering.
            If the argument is None, we follow the Matplotlib default,
            which is {origin}='upper' = [0,0] in the upper left corner.
            Alternative: {origin}='lower' = [0,0] is in the lower left corner.
        out_file : str, optional, default is None
            Name of the output file.
            If the argument is not None, the plot is saved to {output_file}.
        out_dpi : int, optional, default is 300
            Resolution of the output file.

        Returns
        -------
        None
            The plot is shown in the stdout
            and saved to {output_file} if requested.
        '''       
        Data2D.show(self, 
            title=title, icut=icut, cmap=cmap, 
            origin=origin, out_file=out_file, out_dpi=out_dpi)
        
        
    def save(self, out_file, itype=None, icut=None):
        '''
        Save the {Diffractogram} object as an image file.

        * This function is a wrapper of ediff.io.Data2D.save
        
        Parameters
        ----------
        out_file : str or PathLike object
            The name of the file into which the image/array will be saved.
            The image format is deduced from the extension of the {out_file}.
        itype : '8bit or '16bit', optional, default is '16bit'
            Bit depth of the saved image (8bit or 16bit grayscale).
            The bit depth is auto-estimated by skimage.io.save
            from {image}.dtype but it can be also set here explicitly.
            The explicit setting is maintained mostly for historical reasons.
        icut : int or float, optional, default is None
            Intensity cut.
            If {icut} = 200, then all intensity values >200 are set to 200.

        Returns
        -------
        None
            The input image/array is saved in {out_file}.
        '''
        Data2D.save(self, out_file, itype, icut)
        

class Plots:

    
    def set_plot_params(size=(8,6), dpi=100, fontsize=8, 
            my_defaults=True, my_rcParams=None):
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
        # (1) Basic arguments -------------------------------------------------
        if size:  # Figure size
            # Convert size in [cm] to required size in [inch]
            size = (size[0]/2.54, size[1]/2.54)
            plt.rcParams.update({'figure.figsize' : size})
        if dpi:  # Figure dpi
            plt.rcParams.update({'figure.dpi' : dpi})
        if fontsize:  # Global font size
            plt.rcParams.update({'font.size' : fontsize})
        # (2) Additional default parameters -----------------------------------
        if my_defaults:  # Defaults => not forbidden by my_defaults=False
            plt.rcParams.update({
                'lines.linewidth'    : 0.8,
                'axes.linewidth'     : 0.6,
                'xtick.major.width'  : 0.6,
                'ytick.major.width'  : 0.6,
                'grid.linewidth'     : 0.6,
                'grid.linestyle'     : ':'})
        # (3) Further user-defined parameter in rcParams format ---------------
        if my_rcParams:  # Other possible rcParams in the form of dictionary
            plt.rcParams.update(my_rcParams)
    
        
    def plot_final_eld_and_xrd(eld_profile, xrd_profile, fine_tune, Xlim,
            eld_data_label='ED experiment', xrd_data_label='XRD calculation',
            Xlabel='$q$ [1/\u212B]', Ylabel='Intensity',
            Xticks=None, Yticks=None, mXticks=None, mYticks=None,
            out_file=None, out_dpi=300, transparent=False, CLI=False):
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
        fine_tune : float
            The constant for the final fine-tuning of peak position.
            The *fine_tuning* constant has a starting value of 1.000.
            If ELD and XRD peaks are shifted, the constant should be adjusted.
            The constant multiplies the X-values of ELD profile.
        Xlim : tuple of two floats
            The limits for X-axis (minimum and maximu q-vectors on X-axis).
        eld_data_label : str, optional, default is 'ED experiment'
            The label of ELD data (= name of the electron diffraction data).
        xrd_data_label : str, optional, the default is 'XRD calculation'
            The label of XRD data (= name of the X-ray diffraction data).
        Xlabel : str, optional, the default is '$q$ [1/\u212B]' ~ q [1/A]
            The label of X-axis.
        Ylabel : str, optional, the default is 'Intensity'.
            The label of Y-axis.
        Xticks : float, optional, default is None
            The X-axis ticks (if not omitted, use the default).
        Yticks : float, optional, default is None
            The Y-axis ticks (if not omitted, use the default).
        mXticks : float, optional, default is None
            The Y-axis minor ticks (if not omitted, use the default).
        mYticks : float, optional, default is None
            The Y-axis minor ticks (if not omitted, use the default).
        out_file : str or PathLike object, optional, default is None
            The name of the (optional) output file (preferably a PNG-file).
            If {out_file} = None (default), the plot is just shown, not saved.
        out_dpi : int, optional, default is 300
            The DPI of the output file = plot = image of ELD and XRD profiles.
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
        # TODO: more general reading of profiles?
        # => not necessary at the moment
        eld = eld_profile
        xrd = xrd_profile
    
        # Plot the data
        plt.plot(xrd.q, xrd.I, label=xrd_data_label)
        plt.plot(eld.q * fine_tune, eld.I, color='red', label=eld_data_label)
        # Define axis labels
        plt.xlabel(Xlabel)
        plt.ylabel(Ylabel)
        # Define xlim = x-limits = x-range
        plt.xlim(Xlim)
        # Define ticks and minor ticks if requested
        if Xticks is not None:
            plt.gca().xaxis.set_major_locator(plt.MultipleLocator(Xticks))
        if Yticks is not None:
            plt.gca().yaxis.set_major_locator(plt.MultipleLocator(Yticks))
        if mXticks is not None:
            plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(mXticks))
        if mYticks is not None:
            plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(mYticks))
        # Add legent (considering transparency)
        if transparent == True:
            plt.legend(framealpha=0.0)
        else: plt.legend()
        # Additional parameters
        plt.grid()
        plt.tight_layout()
        # Save plot if requested
        if out_file is not None:
            if transparent == True:
                plt.savefig(out_file, dpi=out_dpi, transparent=True)
            else:
                plt.savefig(out_file, dpi=out_dpi, facecolor='white')
        # Show the plot
        if CLI == False: plt.show()
    
    
    def plot_radial_distributions(
            data_to_plot, xlimit, ylimit, output_file=None):
        """
        Plot one or more 1D-radial distrubution files in one graph.
        
        * This is a special/specific function.
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
                - Note2: 1D-profile = a text file/array with 2 cols/rows =
                  [R,I], where R = distance from center, I = diff.intensity.
            * linestyle = matplotlib.pyplot format, such as 'r-' (red line)
            * label = name of the data, which will appear in the plot legend
        xlimit : int
            maximum of the X-axis
        ylimit : int
            maximum of the Y-axis
        output_file : int, optional, default=None
            Name of the output file;
            is given, the plot is shown on screen AND saved in {out_file}. 
    
        Returns
        -------
        None
            The plot is shown in the stdout
            and saved to {out_file} if requested.
        
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
        >>>     out_file='sums_final_1d.png')
        """
        # Initialize
        n = len(data_to_plot)
        rdist = data_to_plot
        # Plot radial distribution files
        for i in range(n):
            # Read data
            data = rdist[i][0]
            if type(data) == str:  # Datafile
                if data.lower().endswith('.png'):  # PNG file, 2D-diffractogram
                    arr = Data2D.read(data)
                    profile = ediff.radial.calc_radial_distribution(arr)
                else:  # TXT file, 1D-profile
                    profile = ediff.radial.read_radial_distribution(data)
            elif type(data) == np.ndarray:  # Numpy array
                if data.shape[0] == data.shape[1]:  # .....sqare array, 2D-diff
                    profile = ediff.radial.calc_radial_distribution(data)
                else:  # .............................non-square array, 1D-diff
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
