'''
Module: ediff.center
--------------------
Find the center of a 2D diffraction pattern.

* The center determination may be surprisingly tricky in certain cases.
* Nevertheless, usually it is enough to call `CenterLocator` as shown below.
* More details and examples at GitHub: https://mirekslouf.github.io/ediff/docs

>>> # Example: How to use CenterLocator and get results?
>>> import ediff as ed
>>>
>>> center = ed.center.CenterLocator(
>>>    input_image='some_diffractogram.png',
>>>    determination='intensity',
>>>    refinement='sum',
>>>    verbose=False, final_replot=True)
>>>
>>> print('Determined center coordinates:', center.x1, center.y1)
>>> print('Refined center coordinates   :', center.x2, center.y2)
'''

# CenterDet, major update and improvements by MS->PS
# PS 2023-10-06: CentDet update, methods compatibility
# MS 2023-11-26: Improved code formatting and docs + TODO notes for PS
# MS 2024-23-09: Re-desing/draft of CenterLocator => CenterLocator_new
# PS 2025-XX-XX: Updates + imporovements in branch {center_determination}
# MS 2025-08-01: Import of {center_determination} branch to {main) 
    

import numpy as np
import skimage as sk
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.legend_handler import HandlerBase
from math import floor
from scipy.optimize import curve_fit
from scipy.ndimage import map_coordinates

import ediff.io 
import ediff.radial 
import os
import cv2

from skimage.measure import moments
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.registration import phase_cross_correlation
import skimage.feature as skf
from scipy.ndimage import gaussian_filter

from scipy.signal import find_peaks
from scipy.fft import _pocketfft, set_backend
from textwrap import dedent

import sys
import warnings
warnings.filterwarnings("ignore")

from functools import wraps
CPU_COUNT = os.cpu_count()


class CenterLocator:
    '''
    CenterLocator object : determine and refine the center of a 2D electron 
    diffraction pattern (diffractogram)
    
    This class detects the center of diffraction patterns. It offers several 
    automatic or manual methods and optionally refines the center position.


    Parameters
    ----------
    input_image : str, path, or numpy.ndarray
        The input 2D diffraction image, either a file path or a NumPy array.

    determination : str or None, optional, default=None
        Method used for the initial center determination.
        Options include:
        - 'manual'   : Manual selection of 3 points on a diffraction ring.
        - 'hough'    : Automatic detection using Hough transform.
        - 'intensity': Automatic detection based on the intensity center 
                       of the central region.
        - 'phase'    : Automatic detection using phase correlation to find 
                       the symmetry center.
        - 'ccorr'    : Automatic detection using cross-correlation to find 
                       the symmetry center.
        - 'curvefit' : Automatic detecting using pseudo-Voigt profile fitting                      

    refinement : str or None, optional, default=None
        Method used for refining the initially detected center.
        Options include:
        - 'manual': Manual fine-tuning of the center along the selected 
                    diffraction ring.
        - 'sum'   : Automatic refinement by maximizing the intensity sum along 
        - 'var'   : Automatic refinement by minimizing intensity variance along 

    rtype : int, default=1
        Radius estimation method.
        Options:
        - 0 : Peak matching method (historical).
        - 1 : Radial distribution method (more universal, new default).
        
    in_file : str or None, optional, default=None
        Path to the file from which previously saved coordinates will be loaded.

    out_file : str or None, optional, default=None
        Path to the file where the determined coordinates will be saved.
        
    sobel : bool, optional, default=False
        If True, apply Sobel filter to the image before processing.
        
    heq : bool or None, optional, default=False
        If True, apply histogram equalization internally (display unchanged).

    icut : float or None, optional, default=None
        Cut-off intensity level for processing the image.
        
    cmap : str, optional, default='gray'
        Matplotlib colormap name used for displaying the image.
        
    csquare : int, optional, defaul=50
        Size (in pixels) of the central square used for initial center search.
        
    cintensity : float, optional, default=0.8
        Threshold intensity for finding the intensity center. Pixels with a 
        (relative) intensity lower than cintensity are ignored.
        
    verbose : int, optional, default=0
        Verbosity level for printing messages during processing:
        - 0 : Silent mode — no messages are printed.
        - 1 : Show help messages during manual refinement only.
        - 2 : Full verbose mode — print all messages and debug information.
        - 3 : Progress mode — print brief status updates during time-consuming 
              processing.
        
    print_sums : bool, optional, default=False
        If True, print intensity sums after each adjustment (in manual modes).
        
    final_print : bool, optional, default=True
        If True, print final coordinates to stdout.
    
    live_plot : bool, optional, default=True
        If True, show live visualization of some processes if available.
 
    Returns
    -------
    None
        The result is stored in the instance variables:
        - `(x1, y1)` for the initially determined center
        - `(x2, y2)` for the refined center (if refinement is applied)
            
    Technical notes
    ---------------
    * The class automatically initializes and uses two internal components:
      `CenterDetermination` and `CenterRefinement`.
    * These internal processes are hidden from the user but can be accessed 
      directly if needed.

    Examples
    --------
    >>> import ediff as ed
    >>>
    >>> locator = ed.center.CenterLocator(
    >>>    input_image   = '.r/path/to/your/data.png',   
    >>>    determination = 'manual',           
    >>>    refinement    = 'sum',
    >>>    final_print   = True)   
                              
    '''
    
    
    def __init__(self,
                 input_image, 
                 determination = None, 
                 refinement = None,
                 rtype = 1,
                 in_file = None,
                 out_file = None,
                 sobel = False,
                 masking = True,
                 ellipse = False,
                 mcorrect = None,
                 heq = False, 
                 icut = None,
                 cmap = 'gray',
                 csquare=50,
                 cintensity=0.8,
                 verbose = 0,
                 print_sums = False,
                 final_print = True,
                 live_plot = False):
        
        ######################################################################
        # PRIVATE FUNCTION: Initialize CenterLocator object.
        # The parameters are described above in class definition.
        ######################################################################

        ## Initialize input attributes ---------------------------------------
        
        # Input image = diffraction pattern
        # (it can be either image file or numpy array)
        self.input_image = input_image
        
        # Methods of center determination, refinement, and radius estimation
        self.determination = determination
        self.refinement = refinement
        self.rtype = rtype
        
        # For reproducibility, convert the method names to lowercase
        if determination is not None:
            self.determination = determination.lower()
        if refinement is not None:
            self.refinement = refinement.lower()
        
        # Text files for reading/saving center coordinates
        self.in_file = in_file
        self.out_file = out_file
        
        # Image processing parameters
        self.ellipse = ellipse
        self.mcorrect = mcorrect
        self.sobel = sobel
        self.masking = masking
        self.heq = heq
        self.icut = icut
        self.cmap = cmap
        self.csquare = csquare
        self.cintensity = cintensity
        
        # Process-Information parameters
        self.verbose = verbose
        self.print_sums = print_sums
        self.final_print = final_print
        self.live_plot = live_plot
        
        ## Initialize new attributes ------------------------------------------
        self.to_refine = []
        self.dText = []
        self.rText = []
        
        # The value of the following attribute can be adjusted
        self.marker_size = 100          

        ## (1) Read input image -----------------------------------------------
        if self.verbose==3:
            print("[INFO] Loading image.")
        # The input image can be numpy.array (ndarray) or image file (path/str)
        if isinstance(input_image, np.ndarray):
            self.image = input_image
        else:
            self.image = ediff.io.read_image(self.input_image)
        
        ## (2) Correct ellipse ------------------------------------------------
        if self.ellipse:
            if self.verbose==3:
                print("[INFO] Correcting ellipse distortion... ", 
                      end="", flush=True)
            self.image = self.ellipse_distortion(self.image, show=True,
                                                 method=self.mcorrect)
            if self.verbose == 3:
                print(" [DONE]")
                
        ## (3) Read center coordinates ----------------------------------------
        # The center coordinates may have been determined in a previous run 
        # of the program and saved to a text file.
        # We can Load coordinates from an input file if specified, this is done 
        # by the following command. As the saved coordinates can be used INSTEAD 
        # of CenterDetermination we will read them here.
        if self.in_file is not None: 
            self.load_results()  
            if self.verbose==3:
                print("[INFO] Loading saved center coordinates.")
        
        ## (4) Run CenterDetermination ----------------------------------------
        if self.verbose == 3:
            print("[INFO] Detecting center...", end="", flush=True)
            
        #  (4a) Initialize/run CenterDetermination
        self.center1 = CenterDetermination(self,
                self.input_image,
                self.determination,
                self.heq,
                self.icut,
                self.cmap,
                self.csquare,
                self.cintensity,
                self.verbose,
                self.print_sums)
        
        if self.verbose == 3:
            print(" [DONE]")
    
        #  (4b) Find radius of a circle/diffraction ring, which is needed 
        #       for the next step = CenterRefimenement.
        if self.verbose==3:
            print("[INFO] Estimating radius of the diffraction pattern...",
                  end="", flush=True)
    
        if self.determination != "manual":
            self.center1.r = self.center1.get_radius(
                self.rtype,
                self.image, 
                self.center1.x, 
                self.center1.y, 
                disp=False)
        
        if self.verbose == 3:
            print(" [DONE]")
            
        ## (5) Initialize/run CenterRefinement --------------------------------
        if self.verbose==3:
            print("[INFO] Refining center coordinates...",
                  end="", flush=True)
                
        self.center2 = CenterRefinement(self,
            self.input_image, 
            self.refinement,
            self.in_file,
            self.out_file,
            self.heq, 
            self.icut,
            self.cmap,
            self.verbose,
            self.print_sums)
        
        if self.verbose == 3:
            print(" [DONE]")
            
            
        ## (6) Collect results ------------------------------------------------
        if self.verbose==3:
            print("[INFO] Collecting results.")
        self.x1 = self.center1.x
        self.y1 = self.center1.y
        
        # Center only detected and not refined
        if self.refinement is not None:
            self.x2 = self.center2.xx
            self.y2 = self.center2.yy
        else: 
            # Center detected and refined
            self.x2 = self.center1.x
            self.y2 = self.center1.y
            self.center2.xx = self.center1.x
            self.center2.yy = self.center1.y
            self.center2.rr = self.center1.r
        
        ## (7) Switch coordinates if necessary --------------------------------
        #      This step is important, as center detection/refinement methods
        #      have various coordinate system origins. The conversion is 
        #      necessary for some combinations of methods
        
        if (self.determination == "intensity"):
            self.x1, self.y1 = self.convert_coords(self.x1, self.y1)
            self.x2, self.y2 = self.convert_coords(self.x2, self.y2)  
        
        if (self.determination =="hough" and self.refinement == "manual"):
            self.x2, self.y2 = self.convert_coords(self.x2, self.y2)  

        if (self.determination =="phase" and self.refinement == "manual"):
            self.x2, self.y2 = self.convert_coords(self.x2, self.y2)  
            
        if (self.determination =="ccorr" and self.refinement == "manual"):
            self.x2, self.y2 = self.convert_coords(self.x2, self.y2) 
                        
        if (self.determination == "curvefit" and self.refinement == "manual"):
            self.x2, self.y2 = self.convert_coords(self.x2, self.y2) 
                            
        ## (8) Print the coordinates if required ------------------------------
        if self.final_print:
            self.dText=str(self.dText)
            self.rText=str(self.rText)
            
            print("----------------------------------------------------------")
            print(self.dText.format(float(self.x1),float(self.y1)))
            
            if self.refinement is not None:
                print(self.rText.format(float(self.x2),float(self.y2)))
                        
        ## (9) Save results to a .txt file if specified -----------------------
        
        if out_file is not None:
            self.save_results()
            if self.verbose==3:
                print("[INFO] Saving results.")

    
    def output(self):
        """
        Return the final results of center determination and refinement. 
        
        Returns
        -------
        x1 : float
            x-coordinate of the center determined by the *determination* method.
            
        y1 : float
            y-coordinate of the center determined by the *determination* method.
        
        x2 : float
            x-coordinate of the center after *refinement*, or same as x1 if no 
            refinement was used.
        
        y2 : float
            y-coordinate of the center after *refinement*, or same as y1 if no 
            refinement was used.

            
        Notes
        -----
        - The function always returns four values: (x1, y1, x2, y2).
        - If no refinement method was used (`refinement=None`), then (x2, y2) 
          will be equal to (x1, y1).
        - All values are rounded to one decimal place before returning.
        - Coordinates are internally converted to float to ensure consistency 
          in output.
        """
        
        # (1) refinement=None -------------------------------------------------
        if self.center2.ret == 1:
            # Convert to float
            if type(self.x1) != float:
                self.x1, self.y1, self.x2, self.y2 = \
                    [float(value) for value in (self.x1, self.y1, 
                                                self.x2, self.y2)]
                    
        # (2) refinement!=None ------------------------------------------------
        else:
            # Convert to float
            if type(self.x1) != float:
                self.x1, self.y1 = \
                    [float(value) for value in (self.x1,self.y1)]
            # Define x2,y2 = x1,y1
            self.x2 = self.x1
            self.y2 = self.y1
            
        # (3) Return (rounded) values of center coordinates -------------------
        return (np.round(self.x1,1), np.round(self.y1,1), 
                np.round(self.x2,1), np.round(self.y2,1))  

    
    def save_results(self):
        '''
        Save the determined and refined center coordinates to a file.
    
        This method writes the coordinates (x1, y1) from the *determination*
        step and (x2, y2) from the *refinement* step to a file specified by the 
        `out_file` attribute.
    
        - If the file already exists, the results are appended to the end.
        - If the file does not exist, it is created and the results are written.
        - Coordinates are formatted to four decimal places for clarity.
        
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        Notes
        -----
        - If `self.out_file` is None, the method does nothing.
        - The results are written in the format:
            x1: <value>, y1: <value>
            x2: <value>, y2: <value>
        - This method does not raise an error if the path is invalid;
          ensure `out_file` is a valid writable path.
        '''
        
        if self.out_file is not None:
            # Check if the specified file exists
            if os.path.isfile(self.out_file):
                # Append results to in_file
                with open(self.out_file, 'a') as f:  # Open in append mode
                    f.write(f"x1: {self.x1:.4f}, y1: {self.y1:.4f}\n")
                    f.write(f"x2: {self.x2:.4f}, y2: {self.y2:.4f}\n")
            else:
                # If the file does not exist, create it and write results
                with open(self.out_file, 'w') as f:  # Open in write mode
                    f.write(f"x1: {self.x1:.4f}, y1: {self.y1:.4f}\n")
                    f.write(f"x2: {self.x2:.4f}, y2: {self.y2:.4f}\n")


    def load_results(self, path=None):
        """
        Load center coordinates from a results file.
    
        This method reads a text file (defaulting to `self.in_file` or
        a provided `path`) containing previously saved center coordinates.
        It extracts coordinate pairs (x1, y1) from the determination step and
        (x2, y2) from the refinement step. All sets of values are stored
        internally, with the most recent values assigned to instance variables 
        `self.x1`, `self.y1`, `self.x2`, and `self.y2`.
    
        Parameters
        ----------
        path : str, optional
            Path to the results file. 
            If not provided, the method uses `self.in_file`.
    
        Returns
        -------
        tuple or None
            If `path` is provided and successfully loaded:
                (x1_values, y1_values) : Lists of all loaded x1 and y1 values.
            If loading fails:
                None
    
        Notes
        -----
        - The file is expected to contain lines in the following format:
            x1: <value>, y1: <value>
            x2: <value>, y2: <value>
        - Multiple coordinate entries can be loaded; the last set is stored in
          the instance variables for easy access.
        - If the file doesn't exist or is unreadable, a warning is printed and
          the function returns None.
        """      
        
        if path is not None:
            self.in_file = path
            
        if self.in_file is not None and os.path.isfile(self.in_file):
            with open(self.in_file, 'r') as f:
                lines = f.readlines()
                
                # Initialize lists to hold multiple values
                x1_values = []
                y1_values = []
                x2_values = []
                y2_values = []
    
                for line in lines:
                    # Strip and split each line to get the coordinates
                    parts = line.strip().split(',')
                    if len(parts) == 2:
                        # Extract x1 and y1 or x2 and y2 based on line content
                        if 'x1' in parts[0]:
                            x1 = float(parts[0].split(':')[1].strip())
                            y1 = float(parts[1].split(':')[1].strip())
                            x1_values.append(x1)
                            y1_values.append(y1)
                        elif 'x2' in parts[0]:
                            x2 = float(parts[0].split(':')[1].strip())
                            y2 = float(parts[1].split(':')[1].strip())
                            x2_values.append(x2)
                            y2_values.append(y2)
    
                # Store the last set of values as instance variables 
                self.x1 = x1_values
                self.y1 = y1_values
                self.x2 = x2_values
                self.y2 = y2_values
    
                # Store the last loaded values:
                if x1_values and y1_values and x2_values and y2_values:
                    self.x1 = x1_values[-1]
                    self.y1 = y1_values[-1]
                    self.x2 = x2_values[-1]
                    self.y2 = y2_values[-1]
                    
                if path is not None:
                    return x1_values, y1_values
        else:
            print("Error reading text file with center coordinates!")
            return


    def get_circle_pixels(self, image, cx, cy, r, num_points=360):
        """
        Extract pixel values along a circular path in the image.
    
        Parameters
        ----------
        image : np.ndarray
            2D grayscale image from which pixel values are extracted.
            
        cx : float
            X-coordinate (row) of the circle center.
        
        cy : float
            Y-coordinate (column) of the circle center.
        
        r : float
            Radius of the circle.
        
        num_points : int, optional
            Number of points to sample along the circular path. Default is 360.
    
        Returns
        -------
        pixels : np.ndarray
            1D array of pixel intensity values sampled along the circular path.
        """
        points = np.linspace(0, 2 * np.pi, num_points)
        
        # Generate coordinates for the circle
        x = cx + r * np.cos(points)
        y = cy + r * np.sin(points)
        
        # Use map_coordinates for accurate sub-pixel interpolation
        coords = np.vstack((y, x))
        
        # 'order=1' is linear interpolation
        pixels = map_coordinates(image, 
                                 coords, 
                                 order=1, 
                                 mode='constant', 
                                 cval=0.0)
        
        # Return pixels
        return pixels
     
     
    def intensity_sum(self, im, cx, cy, r):
        """
        Compute the sum of pixel intensities along a circular path.
    
        This function extracts pixel values along a circle defined by the 
        given center coordinates (cx, cy) and radius r, then returns 
        the sum of these intensities.
    
        Parameters
        ----------
        im : np.ndarray
            2D grayscale image from which the intensities are extracted.
        
        cx : float
            X-coordinate (row) of the circle center.
        
        cy : float
            Y-coordinate (column) of the circle center.
        
        r : float
            Radius of the circle.
    
        Returns
        -------
        float
            Sum of pixel intensities along the circular path.
        """
        # Get pixel intensities
        pixels = self.get_circle_pixels(im, cx, cy, r)
        
        # Sum pixel intensities
        return np.sum(pixels)
    
    
    def intensity_variance(self, image, px, py, pr):
        ''' 
        Variance of intensity values of pixels of a diffraction pattern.

        Parameters
        ----------
        image : array of uint8
            image from which the diffraction pattern has been detected.
        
        px : float64
            x-coordinate of the center of the diffraction pattern.
        
        py : float64
            y-coordinate of the center of the diffraction pattern.
        
        pr : float64
            radius of the diffraction pattern.

        Returns
        -------
        s : float64
            intensity variance

        '''
        # Extract pixels on the circle border
        pxc, pyc = self.get_circle_pixels(px, py, pr)
        pxc = np.array(pxc, dtype=int)
        pyc = np.array(pyc, dtype=int)
        
        # Calculate sum using the filtered values
        s = np.var(image[pxc, pyc])
        return s
    
    
    def convert_coords(self, x, y):
        """
        Convert coordinates between numpy and matplotlib systems.
    
        Parameters
        ----------
        x : int or float
            The x-coordinate in the numpy (column index) format.
        
        y : int or float
            The y-coordinate in the numpy (row index) format.
    
        Returns
        -------
        tuple of (int or float, int or float)
            The converted coordinates in matplotlib format, where:
            - First element corresponds to y (new x in matplotlib).
            - Second element corresponds to x (new y in matplotlib).
        """
        return y, x


    def visualize_results(self, csquare=None, out_file=None, out_dpi=200):
        '''
        Visualize diffractogram and its center after
        center determination + refinement.
    
        Parameters
        ----------
        csquare : int, optional, default is None
            If csquare argument is given,
            only the central square of the diffractogram will be plotted;
            the size of the central square will be equal to csquare argument.
        
        out_file : str, optional, default is None
            If out_file is given,
            save the final plot to image named *out_file*.
        
        out_dpi : int, optional, default is 200
            DPI of the output image file;
            this parameter is relevant only if out_file is given.
    
        Returns
        -------
        None
            The output is the image of the diffractogram
            showing also the central coordinates and refinement ring.
        '''
        
        # (0) Collect center coordinates and radii ----------------------------
        x1, y1 = (self.x1, self.y1)
        r1 = self.center1.r

        x2, y2 = (self.x2, self.y2)
        r2 = self.center2.rr

        # (1) Prepare the image of the diffractogram --------------------------
        image = np.copy(self.to_refine)
        if self.icut is not None:
            im = np.where(image > self.icut, self.icut, image)
        else:
            im = np.copy(image)
            
        # (2) Calculate intensity sum or intensity variance -------------------
        if self.refinement == "var":
            dvar = self.intensity_variance(image, x1, y1, r1)
            rvar = self.intensity_variance(image, x2, y2, r2)
            labeld = f'd-center: [{x1:.1f}, {y1:.1f}]\nint-var: {dvar:.1f}'
            labelr = f'r-center: [{x2:.1f}, {y2:.1f}]\nint-var: {rvar:.1f}'
        else:
            dsum = self.intensity_sum(image, x1, y1, r1)
            rsum = self.intensity_sum(image, x2, y2, r2)
            labeld = f'd-center: [{x1:.1f}, {y1:.1f}]\nint-sum: {dsum:.1f}'
            labelr = f'r-center: [{x2:.1f}, {y2:.1f}]\nint-sum: {rsum:.1f}'
            
        # (3) Calculate xy-limits if we want to see only the central square ---
        if csquare is not None:
            xmin, xmax = (0, image.shape[0])
            ymin, ymax = (0, image.shape[1])
            edge = (xmax - csquare) // 2
            xmin += edge
            xmax -= edge
            ymin += edge
            ymax -= edge
    
        # (4) Final plot: show the diffractogram + refinement results ---------
        fig, ax = plt.subplots(layout='constrained')
        ax.imshow(im, cmap=self.cmap)
    
        ax.set_title(f'Center :: {self.determination}/{self.refinement}')
    
        ax.scatter(x1, y1, 
                   label=labeld, 
                   color='gold', 
                   marker='x', 
                   s=60, lw=1)
        c0 = plt.Circle((x1, y1), r1, 
                        color='gold', 
                        fill=False, 
                        label='detected center', 
                        lw=1)
        ax.add_patch(c0)
    
        ax.scatter(x2, y2,
                   label=labelr, 
                   color='red', 
                   marker='x', lw=1, s=60)
        c1 = plt.Circle((x2, y2), r2, 
                        color='red', 
                        fill=False, 
                        label='refined center', 
                        lw=1)
        ax.add_patch(c1)
    
        ax.legend(loc='upper left', frameon=False, 
                  handler_map={Circle: HandlerCircle()}, 
                  bbox_to_anchor=(.98, 1.01), 
                  bbox_transform=ax.transAxes)
    
        # Remove axes (= tick and ticklabels) around the diffractogram
        ax.axis('off')
        
        # Plot only the central square region if requested
        if csquare is not None:
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
        
        # Save the result to output file if requested
        if out_file is not None:
            fig.savefig(out_file, dpi=out_dpi)
        
        # Show the plot on screen
        plt.show(block=False)


    def ellipse_distortion(self, img, show=True, method=None):
        """
        Corrects elliptical distortion in an image by transforming 4 distorted 
        points (which should ideally lie on a circle) to lie on a true circle.
    
        Parameters
        ----------
        img : np.ndarray
            Input image (2D numpy array) containing an elliptical diffraction 
            pattern.
        
        show : bool, optional
            If True, displays the original and corrected images. 
            Default is True.
        
        method : str, optional
            Point selection method: 'manual' or 'auto'.
            - 'manual': User clicks four points on the image.
            - 'auto': Automatically detects the four brightest Bragg spots 
              (not yet implemented).
    
        Returns
        -------
        corrected : np.ndarray
            The distortion-corrected image.
        """
            
        # Helper function
        def refine_to_local_max(img, pt, window=9):
            """
            Refine a point to the local maximum intensity within a neighborhood.
        
            This function searches for the brightest (maximum intensity) pixel
            within a square window centered around a given point. It is typi-
            cally used to improve the accuracy of an estimated center by snap-
            ping it to the nearest local maximum.
        
            Parameters
            ----------
            img : ndarray
                2D array (grayscale) in which the local maximum is to be found.
                
            pt : array-like of float
                Initial (x, y) point coordinates around which the local maximum 
                is searched.
                
            window : int, optional
                Size of the square window for the local search (default is 9).
                The actual window will be `window x window` pixels, 
                centered at `pt`.
        
            Returns
            -------
            refined_pt : ndarray of float
                Refined point coordinates [x, y] corresponding to the brightest 
                pixel in the local window.
        
            Notes
            -----
            - If the window overlaps the image boundaries, it is clipped 
              appropriately.
            - If the window is completely outside the image, the original point 
              is returned.
            - Assumption : `img` uses (row, column) = (y, x) indexing.
            """
            x, y = int(round(pt[0])), int(round(pt[1]))
            r = window // 2
            x_min, x_max = max(0, x - r), min(img.shape[1], x + r + 1)
            y_min, y_max = max(0, y - r), min(img.shape[0], y + r + 1)
    
            local = img[y_min:y_max, x_min:x_max]
            if local.size == 0:
                return np.array([x, y], dtype=np.float32)
    
            dy, dx = np.unravel_index(np.argmax(local), local.shape)
            refined_x = x_min + dx
            refined_y = y_min + dy
            return np.array([refined_x, refined_y], dtype=np.float32)
        
        # Elliptical correction workflow --------------------------------------
        h, w = img.shape
    
        if method != "manual":
            print("Only user-assisted ('manual') correction available.")
            method = "manual"
    
        if method == 'manual':
            # Manual point selection
            plt.imshow(img, cmap=self.cmap)
            plt.title("Click 4 distorted points (in circular order)")
            raw_pts = np.array(plt.ginput(4, timeout=-1), dtype=np.float32)
            plt.close()
    
            # Refine each point to local max
            distorted_pts = np.array([
                refine_to_local_max(img, pt, window=9) for pt in raw_pts
            ], dtype=np.float32)
    
        # Compute center and average radius
        center = np.mean(distorted_pts, axis=0)
        radii = np.linalg.norm(distorted_pts - center, axis=1)
        avg_radius = np.mean(radii)
    
        # Target points evenly on circle
        angles = np.linspace(0, 2*np.pi, 5)[:-1]
        circle_pts = np.stack([
            center[0] + avg_radius * np.cos(angles),
            center[1] + avg_radius * np.sin(angles)
        ], axis=1).astype(np.float32)
    
        # Perspective transform
        H = cv2.getPerspectiveTransform(distorted_pts, circle_pts)
        
        # Make output square by padding the image to max_dim x max_dim
        max_dim = max(h, w)
        
        # Compute offset to center the original image content in the square 
        x_offset = (max_dim - w) // 2
        y_offset = (max_dim - h) // 2
        
        # Create a translation matrix to shift image content to center
        T = np.array([
            [1, 0, x_offset],
            [0, 1, y_offset],
            [0, 0, 1]
        ], dtype=np.float32)
        
        # Combine the original homography with the translation
        H_total = T @ H  # Matrix multiplication
        
        # Apply the combined transformation
        corrected = cv2.warpPerspective(
            img, H_total, (max_dim, max_dim),
            borderMode=cv2.BORDER_CONSTANT,
            borderValue=0
        )

        # Show results
        if show:
            plt.figure(figsize=(12, 5))
            plt.subplot(1, 2, 1)
            plt.title("Original with refined points")
            plt.imshow(img, cmap=self.cmap)
            plt.scatter(*distorted_pts.T, color='cyan', s=40)
    
            plt.subplot(1, 2, 2)
            plt.title("Corrected image")
            plt.imshow(corrected, cmap=self.cmap)
            plt.show()
    
        return corrected


class CenterDetermination:
    '''
    CenterDetermination object - initial detection of a diffractogram center.

    This class is responsible for identifying the center coordinates of 
    a 2D diffraction pattern image using various detection methods specified 
    by the user. It initializes with the input image and other parameters
    that influence the center detection process.

    Parameters
    ----------
    parent : CenterLocator
        Reference to the parent CenterLocator object, allowing access 
        to shared attributes and methods.
        
    input_image : str, path, or numpy.array
        The input image representing a 2D diffraction pattern, provided as 
        a file path or a NumPy array.
        
    determination : str or None, optional, default=None
        Method used for the initial center determination.
        Options include:
        - 'manual'   : Manual selection of 3 points on a diffraction ring.
        - 'hough'    : Automatic detection using Hough transform.
        - 'intensity': Automatic detection based on the intensity center 
                       of the central region.
        - 'phase'    : Automatic detection using phase correlation to find 
                       the symmetry center.
        - 'ccorr'    : Automatic detection using cross-correlation to find 
                       the symmetry center.
        - 'curvefit' : Automatic detecting using pseudo-Voigt profile fitting 
        
    in_file : str or None, optional, default=None
        Path to the file from which previously saved coordinates will be loaded.
    
    out_file : str or None, optional, default=None
        Path to the file where the determined coordinates will be saved.
        
    sobel : bool, optional, default=False
        If True, apply Sobel filter to the image before processing.
        
    heq : bool or None, optional, default=False
        If True, apply histogram equalization internally (display unchanged).
    
    icut : float or None, optional, default=None
        Cut-off intensity level for processing the image.
        
    cmap : str, optional, default='gray'
        Matplotlib colormap name used for displaying the image.
        
    csquare : int, optional, defaul=50
        Size (in pixels) of the central square used for initial center search.
        
    cintensity : float, optional, default=0.8
        Threshold intensity for finding the intensity center. Pixels with a 
        (relative) intensity lower than cintensity are ignored.
        
    verbose : bool, optional, default=False
        If True, print informational verbose during the detection process.
        
    print_sums : bool, optional, default=False
        If True, print intensity sums after each adjustment (in manual modes).
        
    

    Returns
    -------
    None

    Notes
    -----
    - The class preprocesses the input image and then applies the specified 
      center detection method.
    - The detected center coordinates are stored in instance variables
      `x`, `y`, and `r`, representing the center's x-coordinate, 
      y-coordinate, and radius, respectively.
      
    '''
    
    def __init__(self, parent,
                 input_image,
                 determination = None, 
                 in_file = None,
                 out_file = None,
                 heq = False, 
                 icut = None,
                 cmap = 'gray',
                 csquare=50,
                 cintensity=0.8,
                 verbose = 0,
                 print_sums = False):

        
        ######################################################################
        # PRIVATE FUNCTION: Initialize CenterLocator object.
        # The parameters are described above in class definition.
        ######################################################################
        
        ## (0) Initialize input attributes -----------------------------------
        self.parent = parent
        
        ## (1) Initialize new attributes -------------------------------------
        self.step=0.5 
        
        ## (2) Run functions -------------------------------------------------
        #  (2a) Preprocess data
        self.preprocess(preInit=1)
        
        #  (2b) Center detection methods
        if determination.lower() == "manual":
            self.x, self.y, self.r = self.detection_3points()
            
        elif determination.lower() == "hough":
            self.x, self.y, self.r = self.detection_Hough(
                sobel=self.parent.sobel,
                live_plot=self.parent.live_plot)
            
        elif determination.lower() == "intensity":
            self.x, self.y, self.r = self.detection_intensity(
                self.parent.csquare, 
                self.parent.cintensity,
                sobel=self.parent.sobel)
            
        elif determination.lower() == "phase":
            self.y, self.x, self.r = self.detection_phase(
                sobel=self.parent.sobel,
                masking=self.parent.masking)
         
        elif determination.lower() == "ccorr":
            self.y, self.x, self.r = self.detection_crosscorr(
                sobel=self.parent.sobel,
                disp=False)
        
        elif determination.lower() == "curvefit":
            self.x, self.y, self.r = self.detection_curvefit()
                
        else: 
            print("Selected determination method does not exist.")
            sys.exit()


    def preprocess(self, preInit=0, preHough=0, preManual=0,
                   preVar=0, preSum=0, preInt=0):  
        """
        Preprocess the input image based on the selected automatic detection 
        and refinement methods.
    
        Parameters
        ----------
        preInit : bool, optional, default=0
            Preprocess the original image using initialization methods such as
            histogram equalization (heq) or contrast clipping (icut). This
            is used to enhance the input image prior to any detection or
            correction.
        
        preHough : bool, optional, default=0
            Perform preprocessing required for Hough transform-based automatic 
            detection. This includes edge detection and handling beamstoppers. 
    
        preManual, preVar, preSum, preInt : bool, optional, default=0
            Placeholder flags for future preprocessing needs for manual 
            detection and different refinement methods (variance, sum, 
            intensity-based). Currently unused.
    
        Returns
        -------
        edges : np.ndarray of bool, optional
            Edge map obtained using the Canny edge detector, required for Hough
            transform. Only returned if `preHough` is True.
    
        Notes
        -----
        This function performs different preprocessing steps depending on the
        detection/refinement strategy selected. It is modular and supports
        optional display of intermediate results if `self.parent.verbose` is True.
        """
    
        # Flags
        control_print = 1
        
        # Load original image
        if len(self.parent.image.shape)!=2:
            self.parent.image = np.mean(self.parent.image,axis=2)
            
        image = np.copy(self.parent.image)
        
        ### After initialization: perform an image enhancement if specified
        if preInit == 1:
            # Enhance diffraction pattern to make it more visible
            if self.parent.heq == 1:
                if self.parent.verbose==2:
                    print("Histogram equalized.")
                image = sk.exposure.equalize_adapthist(image)

            self.parent.to_refine = image
            return
        
        ### Hough transform: perform pre-processing necessary for the detection
        if preHough == 1:           
            if self.parent.heq == 0:
                csq = self.central_square(self.parent.to_refine, csquare=80)   
                
                # Beam stopper present in image
                if np.median(csq)<100 and np.median(csq) > 0:
                    if self.parent.verbose==2:
                        print('Beamstopper removed.')
                    max_indices=np.where(self.parent.to_refine>np.median(csq))
        
                    row_idx = max_indices[0]
                    col_idx = max_indices[1]
        
                    self.parent.to_refine[row_idx, col_idx] = 0    
                    
                    max_indices = \
                        np.where(self.parent.to_refine < 0.8*np.median(csq))
                    row_idx = max_indices[0]
                    col_idx = max_indices[1]
        
                    self.parent.to_refine[row_idx, col_idx] = 0   
                    
                    # Detect edges using the Canny edge detector
                    edges = sk.feature.canny(self.parent.to_refine, 
                            sigma=0.2, 
                            low_threshold=2.5*np.median(self.parent.to_refine), 
                            high_threshold=3*np.median(self.parent.to_refine))
                    
                    # Dilate the edges to connect them
                    selem = sk.morphology.disk(5)
                    dilated_edges = sk.morphology.dilation(edges, selem)
                    
                    # Erode the dilated edges
                    # to reduce thickness and smooth the contour
                    connected_edges=sk.morphology.erosion(dilated_edges,selem)

                    if control_print == 1:
                        fig, ax = plt.subplots(nrows=2, ncols=2)
                        ax[0,0].imshow(self.parent.image, origin="upper")
                        ax[0,0].set_title("Original image")
                        ax[0,1].imshow(self.parent.to_refine, origin="upper")
                        ax[0,1].set_title("Hough pre-processed")
                        ax[1,0].imshow(edges, origin="upper")
                        ax[1,0].set_title("Edges")
                        ax[1,1].imshow(connected_edges, origin="upper")
                        ax[1,1].set_title("Connected edges")
                        plt.tight_layout()
                        plt.show(block=False)
                        
                # No beam stopper in image
                else:
                    # Detect edges using the Canny edge detector
                    print('No beamstopper.')
                    edges = sk.feature.canny(self.to_refine, 
                                             sigma=0.2, 
                                             low_threshold=80, 
                                             high_threshold=100)
                    
                    # Dilate the edges to connect them
                    selem = sk.morphology.disk(10)
                    dilated_edges = sk.morphology.dilation(edges, selem)
                    
                    # Erode the dilated edges
                    # to reduce thickness and smooth the contour
                    connected_edges=sk.morphology.erosion(dilated_edges,selem)
                    connected_edges = sk.morphology.remove_small_objects(
                        connected_edges, 
                        min_size=100)
                    
                    if control_print == 1:
                        fig, ax = plt.subplots(nrows=2, ncols=2)
                        ax[0,0].imshow(self.image, origin="upper")
                        ax[0,0].set_title("Original image")
                        ax[0,1].imshow(self.to_refine, origin="upper")
                        ax[0,1].set_title("Hough pre-processed")
                        ax[1,0].imshow(edges, origin="upper",)
                        ax[1,0].set_title("Edges")
                        ax[1,1].imshow(connected_edges,origin="upper")
                        ax[1,1].set_title("Connected edges")
                        plt.tight_layout()
                        plt.show(block=False)
                
            elif self.heq == 1: 
                # Central square extraction
                csq = self.central_square(self.to_refine, csquare=80)   

                # Beam stopper present in image
                if 0.4 <= np.median(csq) <= 0.6:
                    
                    max_indices = \
                        np.where(self.to_refine > 2*np.median(self.to_refine))
    
                    row_idx = max_indices[0]
                    col_idx = max_indices[1]
    
                    self.to_refine[row_idx, col_idx] = 0 
                                                     
                    # Detect edges using the Canny edge detector
                    edges = sk.feature.canny(
                        self.to_refine, 
                        sigma=0.2, 
                        low_threshold=1.5*np.median(self.to_refine), 
                        high_threshold=3*np.median(self.to_refine))

                    
                    if control_print == 1:
                        fig, ax = plt.subplots(nrows=2, ncols=2)
                        ax[0,0].imshow(self.image, origin="upper")
                        ax[0,0].set_title("Original image")
                        ax[0,1].imshow(self.to_refine, origin="upper")
                        ax[0,1].set_title("Hough pre-processed")
                        ax[1,0].imshow(edges, origin="upper")
                        ax[1,0].set_title("Edges")
                      #  ax[1,1].imshow(connected_edges)
                        ax[1,1].set_title("Connected edges")
                        plt.tight_layout()
                        plt.show(block=False)
                    
                # No beam stopper in image
                else:
                    # Detect edges using the Canny edge detector
                    edges = sk.feature.canny(
                        self.to_refine, 
                        sigma=0.2, 
                        low_threshold=2.5*np.median(self.to_refine), 
                        high_threshold=3*np.median(self.to_refine))

                    
                    if control_print == 1:
                        fig, ax = plt.subplots(nrows=2, ncols=2)
                        ax[0,0].imshow(self.image, origin="upper")
                        ax[0,0].set_title("Original image")
                        ax[0,1].imshow(self.to_refine, origin="upper")
                        ax[0,1].set_title("Hough pre-processed")
                        ax[1,0].imshow(edges, origin="upper")
                        ax[1,0].set_title("Edges")
                        ax[1,1].set_title("Connected edges")
                        plt.tight_layout()
                        plt.show(block=False)
            
            return edges
        

    def sobel_filt(self, image, disp=True):   
        """
        The Sobel filter, which is a simple 3×3 kernel convolved with the image 
        which approximates the gradient of intensity in a certain direction. 
        Since the Sobel filter is a simple matrix, it can be applied across 
        datasets without concern for effects caused by changing user defined 
        thresholds, and is much faster and impartial than an iterative method

        Parameters
        ----------
        image : TYPE
            DESCRIPTION.
        disp : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        sobel_magnitude : TYPE
            DESCRIPTION.

        """
        # Apply Sobel filter
        sobel_x = cv2.Sobel(image, 
                            cv2.CV_64F, 
                            1, 0, 
                            ksize=3)  # Horizontal edges
        sobel_y = cv2.Sobel(image, 
                            cv2.CV_64F, 
                            0, 1,
                            ksize=3)  # Vertical edges
        
        # # Take absolute values to remove negatives
        # sobel_x = np.abs(sobel_x)
        # sobel_y = np.abs(sobel_y)

        
        # Compute gradient magnitude
        sobel_magnitude = np.sqrt(sobel_x**2 + sobel_y**2)
        
        # Normalize
        sobel_magnitude = \
            np.uint8(255 * sobel_magnitude / np.max(sobel_magnitude))  
        
        if disp:
            # Display results
            plt.figure(figsize=(12, 4))
            plt.subplot(1, 3, 1), 
            plt.imshow(sobel_x, cmap='gray'), 
            plt.title('Sobel X')
            plt.subplot(1, 3, 2),
            plt.imshow(sobel_y, cmap='gray'),
            plt.title('Sobel Y')
            plt.subplot(1, 3, 3), 
            plt.imshow(sobel_magnitude, cmap='gray'), 
            plt.title('Gradient Magnitude')
            plt.show()
        
        return sobel_magnitude
    
    
    def detection_phase(self, normalize_bg=True, sobel=False, 
                        disp=False, masking=True):
        """
        Detects the center of symmetry in a diffraction image using a 
        combination of weighted intensity averaging and phase cross-correlation.
    
        This method performs the following steps:
            1. Estimates an initial center using a weighted average of pixel 
               intensities.
            2. Crops the image to focus on its central region.
            3. Optionally normalizes the background to mitigate experimental 
               artifacts.
            4. Computes a refined center via inversion symmetry and phase 
               cross-correlation.
            5. Returns the final refined center coordinates.
    
        Parameters
        ----------
        normalize_bg : bool, optional, default=True
            If True, normalize background intensity to reduce artifacts.
            
        sobel : bool, optional, default=False
            If True, applies Sobel filtering (currently unused if not implemented)
            
        disp : bool, optional, default=False
            If True, displays intermediate processing results. 
            
        masking : bool, optional, default=True
            If True, applies masking to exclude unwanted regions from analysis
    
        Returns
        -------
        r_final : float
            Estimated row-coordinate of the image center after refinement.
            
        c_final : float
            Estimated column-coordinate of the image center after refinement.
            
        None : NoneType
            Placeholder for future functionality (e.g., radius estimation).
        """
        # (1) Get image  ------------------------------------------------------
        image = np.copy(self.parent.to_refine)
        
        # if Sobel==True, perform Sobel filtration
        # https://doi.org/10.1016/j.ultramic.2016.12.021
        if sobel:
            image = self.sobel_filt(image,disp)
        
        
        # shift the pixel values (smallest value is zero -> non-negative image)
        image -= image.min()

        # (2) Prepare mask and compute weighted center ------------------------
        if masking:
            mask = self.auto_masking(image) 
            # print("masked")
        else: mask = np.ones_like(image, dtype=bool) # binary mask
        
        weights = image                        # 
        
        rr, cc = np.indices(image.shape)       # matrices of row/col indices

        # Weighted center of mass (brighter pixels contribute more)
        # This gives an initial estimate of the center (r_, c_).
        r_ = int(np.average(rr.flatten(), weights=weights.flatten()))
        c_ = int(np.average(cc.flatten(), weights=weights.flatten()))

        # (3) Crop the image around the estimated center ----------------------
        # Some diffraction patterns are not centered, and so there's a lot of 
        # image area that cannot be used for registration.
        # Radial inversion becomes simple inversion of dimensions
        side_length = \
            floor(min([r_,abs(r_-image.shape[0]),c_,abs(c_-image.shape[1])]))
        rs = slice(r_ - side_length, r_ + side_length)
        cs = slice(c_ - side_length, c_ + side_length)
        image = image[rs, cs]
        mask = mask[rs, cs]
    
        # (4) Normalize background intensity (optional) -----------------------
        # This step removes intensity variations due to experimental artifacts, 
        # such as uneven illumination.
        if normalize_bg:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                image = image.astype(np.float64)  # float64 before division
                image /= gaussian_filter(
                    input=image, sigma=min(image.shape) / 25, truncate=2)
        
        # Replace NaN values with 0s
        image = np.nan_to_num(image, copy=False)


        # (5) Compute inversion symmetry using Fourier methods ----------------
        # an inverted copy of the image (im_i) and mask (mask_i), 
        # flipping both horizontally and vertically ::-1
        im_i = image[::-1, ::-1]
        mask_i = mask[::-1, ::-1]

        # downsampling to speed up processing (for too large images)
        downsampling = 1
        if min(image.shape) > 1024:
            downsampling = 2
    
        # phase cross-correlation finds the best alignment between the image
        # and its inverted copy
        # the shift tells us how much the center is offset from perfect 
        # inversion symmetry.
        shift, *_ = self.calc_fft(phase_cross_correlation)(
            reference_image=image[::downsampling, ::downsampling],
            moving_image=im_i[::downsampling, ::downsampling],
            reference_mask=mask[::downsampling, ::downsampling],
            moving_mask=mask_i[::downsampling, ::downsampling],
        )
        
        # The computed shift is rescaled to match the original image resolution.
        correction = shift * downsampling


        # (6) User information (if required) ----------------------------------
        if (self.parent.verbose or self.parent.final_print):
            self.parent.dText = \
                "Center Determination (PhaseCorrCenter): ({:.3f}, {:.3f})"
        
        # Ensure correction is a tuple or list before extracting elements
        correction_r, correction_c = correction \
            if isinstance(correction, (list, tuple, np.ndarray)) \
                else (correction, correction)
        
        r_final = r_ + correction_r / 2
        c_final = c_ + correction_c / 2
                
        return float(r_final), float(c_final), None
    

    def detection_crosscorr(self, normalize_bg=True, sobel=False, disp=False):
        """
        Detects the center of symmetry in a diffraction image using a 
        combination of weighted intensity averaging and cross-correlation.
    
        This method performs the following steps:
            1. Estimates an initial center using the weighted average of pixel
               intensities.
            2. Crops the image to focus on its central region.
            3. Optionally normalizes the background to reduce experimental
               artifacts.
            4. Computes a refined center by analyzing inversion symmetry 
               using cross-correlation.
            5. Returns the final corrected center coordinates.
    
        Parameters
        ----------
        normalize_bg : bool, optional, default=True
            If True, normalizes the background intensity to reduce illumination 
            artifacts. 
            
        sobel : bool, optional, default=False
            If True, applies Sobel filtering to enhance edges prior to processing. 
            
        disp : bool, optional, default=False
            If True, displays intermediate processing results.
    
        Returns
        -------
        r_final : float
            Estimated row-coordinate of the image center after refinement.
            
        c_final : float
            Estimated column-coordinate of the image center after refinement.
        """
       
        # (1) Get image -------------------------------------------------------
        image = np.copy(self.parent.to_refine)
        
        # Apply Sobel filter if requested
        if sobel:
            image = self.sobel_filt(image, disp=False)
        
        # Shift pixel values to ensure non-negative image
        image -= image.min()
     
    
        # (2) Normalize background intensity (optional) -----------------------
        if normalize_bg:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                image = image.astype(np.float64) 
                image /= gaussian_filter(image, 
                                         sigma=min(image.shape) / 25, 
                                         truncate=2)
        
        # Replace NaN values with 0s
        image = np.nan_to_num(image, copy=False)
    
        # (3) Compute inversion symmetry using cross-correlation --------------
        # Dynamically set the margin to a reasonable fraction of the image size
        margin = int(np.floor(image.shape[0] // 4))
        im_i = image[margin:-margin, margin:-margin]  # Crop slightly
        
        template = im_i[::-1, ::-1]  # Inverted image for symmetry comparison
    
        # Compute cross-correlation using OpenCV template matching
        # it finds areas in an image that best match a given template. 
        # it compares the original image with its inverted version to measure 
        # inversion symmetry
        # im_i (template) slides over image in a pixel-by-pixel manner.
        # at each position, it computes the cross-correlation between 
        # overlapping pixels.
        # result is a 2D matrix, where each value represents the correlation 
        # score at a particular location. Higher values indicate
        # better symmetry at that position.
        
        # Ensure template is smaller than image
        assert template.shape[0] <= image.shape[0] \
            and template.shape[1] <= image.shape[1], \
            "Template is larger than the image"
        
        result = cv2.matchTemplate(image.astype(np.float32), 
                                   template.astype(np.float32), 
                                   method=cv2.TM_CCORR_NORMED)
        
        # Get the location of the best match in `result`
        max_location = np.unravel_index(np.argmax(result), result.shape) 
        correction_r, correction_c = max_location  
        
        if disp: 
            plt.figure()
            plt.imshow(result, cmap="inferno")
            plt.title("Correlation matrix")
            plt.show()
            plt.axis("off")
            plt.tight_layout()
            plt.show()
        
        # Adjust to original image coordinates (account for template cropping)
        corrected_r = correction_r + template.shape[0] // 2
        corrected_c = correction_c + template.shape[1] // 2
        
        # (4) Visualization ---------------------------------------------------
        # Plot with correct alignment
        if disp:
            
            plt.figure()
            if self.parent.icut:
                orig=np.where(self.parent.to_refine>self.parent.icut,
                              self.parent.icut, self.parent.to_refine)
            else: orig = np.copy(self.parent.to_refine)
            
            # Assuming fft_im_i is your Fourier transform image
            center_size = 513
            center_x = orig.shape[1] // 2
            center_y = orig.shape[0] // 2

            # Crop the central part of the spectrum
            orig = orig[
                center_y - center_size // 2:center_y + center_size // 2,
                center_x - center_size // 2:center_x + center_size // 2]
            
            plt.imshow(orig, cmap="gray", origin="lower")
            plt.scatter(corrected_c-center_x//2, corrected_r-center_y//2, 
                        color="white", marker="x", s=100,
                        label = "center")
            my_title =  "Location of the highest correlation coefficient"
            my_title += "in the original image."
            plt.title(my_title)
            plt.axis("off")
            plt.tight_layout()
            
            plt.show()
        
        # (5) User information (if required) ----------------------------------
        if self.parent.verbose or self.parent.final_print:
            self.parent.dText = \
                "Center Determination (CrossCorrCenter): ({:.3f}, {:.3f})"

        return float(corrected_r), float(corrected_c), None

    
    def detection_intensity(self, 
                            csquare, cintensity, plot_results=0, sobel=0):
        '''
        Find center of intensity/mass of an array.
        
        Parameters
        ----------
        arr : 2D-numpy array
            The array, whose intensity center will be determined.
            
        csquare : int, optional, default is 20
            The size/edge of the square in the (geometrical) center.
            The intensity center will be searched only within the central
            square.
            Reasons: To avoid other spots/diffractions and
            to minimize the effect of possible intensity assymetry
            around the center. 
            
        cintensity : float, optional, default is 0.8
            The intensity fraction.
            When searching the intensity center, we will consider only
            pixels with intensity > max.intensity.
            
        Returns
        -------
        xc, yc : float,float
            XY-coordinates of the intensity/mass center of the array.
            Round XY-coordinates if you use them for image/array calculations.
        '''  
        
        # (1) Get image/array size --------------------------------------------
        image = np.copy(self.parent.to_refine)
        
        if sobel:
            image = self.sobel_filt(image, disp=False)
             
        arr = np.copy(image)
        xsize,ysize = arr.shape
        
        # (2) Calculate borders around the central square ---------------------
        xborder = (xsize - csquare) // 2
        yborder = (ysize - csquare) // 2
        
        # (3) Create the central square, --------------------------------------
        # from which the intensity center will be detected
        arr2 = arr[xborder:-xborder,yborder:-yborder].copy()
        
        # (4) In the central square, ------------------------------------------
        # set all values below cintenstity to zero
        arr2 = np.where(arr2>np.max(arr2)*cintensity, arr2, 0)
        
        # (5) Determine the intensity center from image moments ---------------
        # see image moments in...
        # skimage: https://scikit-image.org/docs/dev/api/skimage.measure.html
        # wikipedia: https://en.wikipedia.org/wiki/Image_moment -> Centroid
        # ---
        # (a) Calculate 1st central moments of the image
        M = moments(arr2,1)
        # (b) Calculate the intensity center = centroid according to www-help
        (self.x, self.y) = (M[1,0]/M[0,0], M[0,1]/M[0,0])
        # (c) We have centroid of the central square
        # but we have to recalculate it to the whole image
        (self.x, self.y) = (self.x + xborder, self.y + yborder)
        # (d) Radius of a diffraction ring 
        self.r = self.get_radius(
            self.parent.rtype, 
            image, 
            self.x, 
            self.y,
            disp=False)
        
        # (6) User information (if required) ----------------------------------
        if (self.parent.verbose or self.parent.final_print):
            self.parent.dText =\
                "Center Determination (IntensityCenter): ({:.3f}, {:.3f})"

        # (7) Plot results (if required) --------------------------------------
        if plot_results == 1:
           self.visualize_center(self.x, self.y, self.r) 
                
        # (8) Return the results: ---------------------------------------------
        #  a) XY-coordinates of the center
        #  b) radius of the circle/diff.ring,
        #     from which the center was determined

        return(self.x, self.y, self.r)
    

    def detection_Hough(self, 
                        plot_results=False, sobel=False, live_plot=False):
        '''        
        Perform Hough transform to detect the center of diffraction patterns 
        with optional real-time visualization and a final image showing all 
        detected circles.
    
        Parameters
        ----------
        plot_results : int, binary
            If 1, shows the final detected circle.
        
        sobel : bool, optional, default=False
            If True, applies Sobel filtering to enhance edges prior to processing. 
            
        sobel            
        live_plot : bool, optional, default=False
            If True, animates the circle detection process and shows all 
            detected circles.
    
        Returns
        -------
        self.x : float64
            x-coordinate of the detected center.
            
        self.y : float64
            y-coordinate of the detected center.
            
        self.r : float64
            Radius of the detected circle.
        '''
        
        ## (0) Image preprocessing --------------------------------------------
        im = np.copy(self.parent.to_refine)
                    
        if self.parent.heq == 0:
            if np.sum(im) < 150000:
                im[im > 50] = 0  
            edges = skf.canny(im, 
                              sigma=0.2, 
                              low_threshold=80, 
                              high_threshold=100)
        elif self.parent.heq == 1:
            if np.sum(im) > 40000:
                im[im > 50] = 0  
            edges = skf.canny(im, 
                              sigma=0.2, 
                              low_threshold=0.80, 
                              high_threshold=1)
        if sobel:
            edges = self.sobel_filt(im, disp=False)

        ## (1) Define the radii range for concentric circles ------------------
        min_radius, max_radius, radius_step = 50, 120, 10
        radii = np.arange(min_radius, max_radius + radius_step, radius_step)
    
        # Store detected circles
        detected_circles = []
    
        if live_plot:
            if self.parent.icut is not None:
                image = np.where(im > self.parent.icut, self.parent.icut, im)
            else:
                image = np.copy(im)
                
            fig, ax = plt.subplots()
            plt.axis("off")
            ax.imshow(image, cmap=self.parent.cmap)
            plt.title("Live Hough Circle Detection")
            plt.ion()  # Turn on interactive mode
    
        ## (2) Process each radius one at a time ------------------------------
        for r in radii:
            hough_res = hough_circle(edges, np.array([r]))
            accums, x_peaks, y_peaks, _ = hough_circle_peaks(hough_res, 
                                                             np.array([r]), 
                                                             total_num_peaks=1)
    
            # Store circle if detected
            if len(x_peaks) > 0:
                detected_circles.append((x_peaks[0], y_peaks[0], r))
    
            # Live animation (if enabled)
            if live_plot:
                ax.clear()
                ax.imshow(image, cmap=self.parent.cmap)
                plt.title(f"Detecting Circles: Radius {r}px")
                for x, y, rad in detected_circles:
                    plt.axis("off")
                    circ = plt.Circle((x, y), rad, color='red', 
                                      fill=False, linewidth=1)
                    ax.add_patch(circ)
                plt.pause(0.3)
    
        ## (3) Final detected circle ------------------------------------------
        accums ,self.x,self.y,self.r = hough_circle_peaks(
            hough_circle(edges, radii),
            radii,
            total_num_peaks=5
        )
        
        
        self.x, self.y, self.r = \
            float(self.x[0]), float(self.y[0]), float(self.r[0])
    
        ## (4) User information -----------------------------------------------
        if (self.parent.verbose or self.parent.final_print):
            self.parent.dText = \
                "Center Determination (HoughTransform) : ({:.3f}, {:.3f})".\
                    format(self.x, self.y)
    
        ## (5) Final plot with all detected circles (if enabled) --------------
        if live_plot:
            plt.ioff()  # Turn off interactive mode
            plt.close("all")
        
        
            fig_final, ax_final = plt.subplots()
            ax_final.imshow(image, cmap=self.parent.cmap)
            plt.title("Hough Transform Preliminary Detection")
        
            for x, y, rad in detected_circles:
                circ = plt.Circle((x, y), rad, color='blue', 
                                  fill=False, linewidth=1)
                ax_final.add_patch(circ)
        
            # Highlight the final detected circle (red)
            circ_final = plt.Circle((self.x, self.y), self.r, 
                                    color='red', fill=False, 
                                    linewidth=2, label="Selected Circle")
            plt.scatter(self.x, self.y, color="red", marker='x', s=60, lw=1)
            ax_final.add_patch(circ_final)
        
            # Create a single blue circle just for the legend
            legend_blue = plt.Line2D([0], [0], 
                                     color='blue', 
                                     linewidth=1, 
                                     linestyle='-', 
                                     label="Searching Space")
        
            # Add the legend
            ax_final.legend(handles=[legend_blue, circ_final])
            plt.axis("off")
            plt.tight_layout()

            plt.show()
    
        # Close live plot (to avoid non-closed windows in Jupyter)
        plt.close("all")
        
        ## (6) Final detected circle visualization (if requested) --------------
        if plot_results:
            self.visualize_center(self.x, self.y, self.r)
    
        return self.x, self.y, self.r

   
    def detection_curvefit(self, plot_results=False):
        """
        Detects the center of a diffraction pattern by fitting a 2D pseudo-Voigt 
        function to the central region of the image.
     
        This method:
            - Crops a square region around the image center.
            - Fits a 2D pseudo-Voigt peak to the cropped region.
            - Returns the refined center coordinates in the original image space.
            - Optionally visualizes the result.
            - Falls back to the crop center if fitting fails.
     
        Parameters
        ----------
        plot_results : bool, optional, default=False
            If True, displays the fitted center location on the image. 
     
        Returns
        -------
        x : float
            Refined x-coordinate of the center (column index).
            
        y : float
            Refined y-coordinate of the center (row index).

        """
    
        def pseudo_voigt_2d(xy, amp, x0, y0, sigma_x, sigma_y, eta):
            """
            Computes a 2D pseudo-Voigt function, which is a linear combination 
            of a Gaussian and a Lorentzian function.
        
            Parameters
            ----------
            xy : tuple of numpy.ndarray
                A tuple (x, y) containing coordinate meshgrids.
                
            amp : float
                Amplitude of the peak.
                
            x0, y0 : float
                Coordinates of the center of the peak.
                
            sigma_x, sigma_y : float
                Width (standard deviation) of the peak along the x and y axes.
                
            eta : float
                Mixing parameter (0 = pure Gaussian, 1 = pure Lorentzian).
        
            Returns
            -------
            numpy.ndarray
                Flattened 1D array of function values evaluated at the input 
                coordinates.
            """
            x, y = xy
            r_x = (x - x0) / sigma_x
            r_y = (y - y0) / sigma_y
    
            gaussian = np.exp(-(r_x**2 + r_y**2) / 2)
            lorentzian = 1 / (1 + r_x**2 + r_y**2)
    
            return amp * ((1 - eta) * gaussian + eta * lorentzian)
        
        # Define crop size around the center of the image
        crop_size = 100
        
        # Make a copy of the image to avoid modifying the original
        image = np.copy(self.parent.to_refine)
        
        # Get the image dimensions
        h, w = image.shape
        
        # Estimate the center of the image
        x_center, y_center = w // 2, h // 2
    
        # Calculate the bounds of the cropped square region
        half_crop = crop_size // 2
        x_min, x_max = max(0, x_center - half_crop), min(w, x_center + half_crop)
        y_min, y_max = max(0, y_center - half_crop), min(h, y_center + half_crop)
    
        # Crop the image around the center for curve fitting
        cropped_image = image[y_min:y_max, x_min:x_max]
    
        # Create a meshgrid of x and y coordinates for the cropped region
        x, y = np.meshgrid(np.arange(cropped_image.shape[1]),
                           np.arange(cropped_image.shape[0]))
        
        # Flatten the meshgrid coordinates and image for use with curve fitting
        x_flat, y_flat = x.ravel(), y.ravel()
        image_flat = cropped_image.ravel()

        # Set initial parameter guesses for curve fitting:
        # [amplitude, x0, y0, sigma_x, sigma_y, eta]    
        initial_guess = [
            np.max(cropped_image),
            crop_size // 2,
            crop_size // 2,
            max(1, crop_size / 5),
            max(1, crop_size / 5),
            0.5
        ]
        
        # Define bounds for each parameter during optimization
        bounds = (
            [0, 0, 0, 1, 1, 0],
            [np.inf, crop_size, crop_size, crop_size, crop_size, 1]
        )
    
        try:
            # Fit the pseudo-Voigt function to the cropped image
            popt, _ = curve_fit(
                pseudo_voigt_2d,
                (x_flat, y_flat),
                image_flat,
                p0=initial_guess,
                bounds=bounds,
                max_nfev=2000  # increase if needed
            )
            
            # Extract refined x and y center coordinates from fit
            x_refined, y_refined = popt[1], popt[2]
            used_fallback = False
    
        except RuntimeError:
            # If fitting fails, fall back to geometric center 
            x_refined, y_refined = crop_size // 2, crop_size // 2
            used_fallback = True
            
        # Add the offset of the cropped region to get coordinates in full image
        self.x, self.y = x_refined + x_min, y_refined + y_min
        
        # Print the result if verbose are enabled
        if self.parent.verbose or self.parent.final_print:
            self.parent.dText = (
                f"Center Determination (CurveFitting) :   ({self.x:.3f}, {self.y:.3f})"
                + (" [fallback]" if used_fallback else "")
            )
            
        # Optionally show the result visually
        if plot_results:
            self.visualize_center(self.x, self.y, 100)
    
        return self.x, self.y, None


    def detection_3points(self, plot_results=0):
        '''         
        Semi-automated detection of the diffraction pattern center using manual 
        selection of 3 points.

        This method allows the user to manually select 3 points along a visible 
        diffraction ring in the image to define a circle. Once selected, 
        the center and radius of the circle are automatically computed. After 
        that, the user can manually adjust the center position using
        an interactive interface.

        User Controls (during point selection)
        --------------------------------------
        - **'1'**: Add a point at the current mouse cursor position 
                         (max 3 points).
        - **'2'**: Delete the most recently added point.
        - **'3'**: Delete the point closest to the mouse cursor.
        - **'d'**: When 3 points are selected (DONE), proceed 
        - **Close**    : Terminate the process without detecting a center.
        
        Parameters
        ----------
        plot_results : int, binary, default=0
            Plot the pattern determined by pixels selected by the user.
        
        Returns
        -------
        self.x : float64
            x-coordinate of the detected center
            
        self.y : float64
            y-coordinate of the detected center
            
        self.r : float64
            radius of the detected center
            (if available, othervise returns None)
                                    
        '''
        # (0) Load and prepare image ------------------------------------------
        im = self.parent.to_refine
        
        # Edit contrast with a user-predefined parameter
        if self.parent.icut is not None:
            if self.parent.verbose==2:
                print("Contrast enhanced.")
            im = np.where(im > self.parent.icut, 
                              self.parent.icut, 
                              im)
            
        # (2) Create a figure and display the image ---------------------------
        fig, ax = plt.subplots(figsize=(12, 12)) 
        
        # Allow using arrows to move back and forth between view ports
        plt.rcParams['keymap.back'].append('left')
        plt.rcParams['keymap.forward'].append('right')
 
        plt.title("Select 3 points defining one of diffraction circles", 
                  fontsize = 20)
        ax.imshow(im, cmap = self.parent.cmap, origin="upper")
        ax.axis('off')

        # (3) User information ------------------------------------------------
        instructions = dedent(
            """
            CenterDetermination :: ThreePoints (semi-automated method)
            Select 3 points to define a diffraction circle using keys:
              - '1' : define a point at current cursor position
              - '2' : delete the last point
              - '3' : delete the point closest to the cursor
              - 'd' : done = finished = go to the next step
            Close the figure to terminate. No center will be detected.
            """)

        if (self.parent.verbose==1 or self.parent.verbose==2):
            print(instructions)
       
        # (4) Enable interactive mode -----------------------------------------
        # (figure is updated after every plotting command
        # (so that calling figure.show() is not necessary
        plt.ion()
 
        # (5) Initialization
        # Initialize the list of coordinates
        self.coords = [] 
        
        # Initialize all flags and counters
        calculate_circle_flag = False          # press 'd' event
        termination_flag = False               # close window event
        point_counter = 0                      # number of selected points
        
        # (6) Define the event handler for figure close event -----------------
        def onclose(event):
            nonlocal termination_flag
            termination_flag = True
            if self.parent.verbose==2:
                print('Execution terminated.')
                print("------------------------------------------------------------")

 
        # Connect the event handler to the figure close event
        fig.canvas.mpl_connect('close_event', onclose)
        
        
        # (7) Define the callback function for key press events ---------------
        def onkeypress(event):
            # nonlocal to modify the flag variable in the outer scope
            nonlocal calculate_circle_flag, point_counter, termination_flag
            
            # Store the zoom level
            current_xlim = ax.get_xlim()
            current_ylim = ax.get_ylim()
 
            ## Delete points -- the closest to the cursor
            if event.key == '3':
                point_counter -= 1
                if len(self.coords) > 0:
                    pointer_x, pointer_y = event.xdata, event.ydata
                    distances = [
                        np.sqrt((x - pointer_x)**2 + (y - pointer_y)**2)
                        for x, y in self.coords]
                    closest_index = np.argmin(distances)
                    del self.coords[closest_index]
    
                    ## Redraw the image without the deleted point
                    ax.clear()
                    ax.imshow(im, cmap = self.parent.cmap,
                              origin="upper")
                    for x, y in self.coords:
                        ax.scatter(x, y, 
                                   c='r', marker='x', 
                                   s=self.parent.marker_size)

                    my_plot_title = (
                        "Select 3 points to define "
                        "one of diffraction circles.")
                    plt.title(my_plot_title, fontsize=20)
                    
                    # Retore the previous zoom level
                    ax.set_xlim(current_xlim)
                    ax.set_ylim(current_ylim)
                    ax.axis('off')
                    
                    fig.canvas.draw()
                else:
                    print("No points to delete.")
    
            # Delete recent point (last added) -- independent on the cursor
            if event.key == '2':
                # Check if there are points to delete
                if point_counter > 0:  
                    point_counter -= 1
                    if len(self.coords) > 0:
                        # Delete the last point in the list
                        del self.coords[-1]
    
                        # Redraw the image without the deleted point
                        ax.clear()
                        ax.imshow(im, cmap=self.parent.cmap,
                                  origin="upper")
                        for x, y in self.coords:
                            ax.scatter(x, y,
                                       c='r', marker='x', 
                                       s=self.parent.marker_size)

                        my_plot_title = (
                            "Select 3 points to define "
                            "one of diffraction circles.")
                        plt.title(my_plot_title)
                        
                        # Retore the previous zoom level
                        ax.set_xlim(current_xlim)
                        ax.set_ylim(current_ylim)
                        ax.axis('off')

                        fig.canvas.draw()
                else:
                    print("No points to delete.")
                    
            ## Select points 
            elif event.key == '1':
                # Only allow selecting up to three points
                if point_counter < 3:  
                    # Save the coordinates of the clicked point
                    new_point = (event.xdata, event.ydata)
                    
                    if new_point in self.coords:
                        # Do not allow multiple selection of one point
                        print("The selected point already exists.")
                    else:
                        # Add selected point
                        self.coords.append(new_point)
    
                        # Visualize the selected point on the image
                        ax.scatter(event.xdata, event.ydata, 
                                   c='r', marker='x', 
                                   s=self.parent.marker_size)

                        # Restore the previous zoom level
                        ax.set_xlim(current_xlim)
                        ax.set_ylim(current_ylim)
                        ax.axis('off')

                        fig.canvas.draw()
    
                        point_counter += 1

    
                if len(self.coords) == 3:
                    # Turn off interactive mode
                    plt.ioff()
    

    
            # (8) Calculate circle or terminate -------------------------------
            elif event.key == 'd':
                if len(self.coords) == 3:
                    calculate_circle_flag = True
 
                else:
                    print("Select exactly 3 points to calculate the circle.")
                    fig.canvas.draw()
    
        # Connect the callback function to the key press event
        cid0 = fig.canvas.mpl_connect('key_press_event', onkeypress)

        # Show the plot 
        plt.tight_layout()
        ax.axis('off')

        plt.show(block=False)
      
        # Wait for 'd' key event or close the figure if no points are selected
        while not calculate_circle_flag and not termination_flag:
 
            try:
                plt.waitforbuttonpress(timeout=0.1)
                # Store the zoom level
                current_xlim = ax.get_xlim()
                current_ylim = ax.get_ylim()
 
                # Plot detected diffraction pattern
                if calculate_circle_flag:
 
                    self.calculate_circle(plot_results=0)
                    
                    ax.clear()
                    ax.imshow(im, cmap = self.parent.cmap,
                              origin="upper")
                    # Retore the previous zoom level
                    ax.set_xlim(current_xlim)
                    ax.set_ylim(current_ylim)
                 
                    circle = plt.Circle(
                        (self.x, self.y), self.r, color='r', fill=False)
                    ax.add_artist(circle)
        
                    # Plot center point
                    center, = ax.plot(self.x, self.y, 'rx', markersize=12)
                    plt.title('Manually adjust the position of the center using keys.')
        
                    # Display the image
                    plt.draw()
                    ax.axis('off')

                    plt.show(block = False)

            except KeyboardInterrupt:
                print("Execution manually interrupted by user.")
                break
            except ValueError as e:
                print("ValueError:", e)
                break
           
        # If the termination_flag is True, stop the code
        if termination_flag: 
             print("No points selected. Returned None values.")
             sys.exit()
             return None, None, None
        
        # Disconnect key press events
        fig.canvas.mpl_disconnect(cid0) 
        
        # local variables save
        self.center = center
        
        self.backip = [self.x, self.y, self.r]
        # Manually adjust the calculated center coordinates
        self.x, self.y, self.r = self.adjustment_3points(fig, circle, center)

        # Return the results:
        # a) XY-coordinates of the center
        # b) radius of the circle/diff.ring,
        #    from which the center was determined
        return(self.x, self.y, self.r)
    
    
    def adjustment_3points(self, fig, circle, center, plot_results=0) -> tuple:
        """
        Interactive refinement of the diffraction pattern center after initial
        detection via 3 points.

        This method allows the user to fine-tune the estimated center and radius
        of a diffraction ring (initially determined using 3 manually selected 
        points) through interactive keyboard controls.

        The user interacts directly with the figure window containing the image 
        and the detected circle.

        User Controls (during point selection)
        --------------------------------------
        - Arrow keys (←, →, ↑, ↓): Move the center left, right, up, or down
        - **'+'** : Increase the radius
        - **'-'** : Decrease the radius
        - **'b'** : Increase step size (×5)
        - **'l'** : Decrease step size (/5, with minimum step size 0.5)
        - **'d'**: Done — finalize the center and radius
        - Closing the figure: Cancels refinement and returns the original input 
          center and radius

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            The figure window used for interactive refinement.
            
        circle : matplotlib.patches.Circle
            The circle object initially defined from 3 selected points.
            
        center : tuple
            The initial (x, y) center coordinates estimated from the selected 
            points.
            
        plot_results : bool, optional, default=False
            If True, the results of the adjustment are visualized.

        Returns
        -------
        xy[0] : float
            x-coordinate of the adjusted center of the diffraction pattern.
            
        xy[1] : float
            y-coordinate of the adjusted center of the diffraction pattern.
            
        r : float
            Adjusted radius of the diffraction ring.

        Notes
        -----
        - Interactive adjustments are visualized in real time.
        - Intensity sum at the current center/radius is printed if print_sums
          is True.
        - Arrow key defaults in Matplotlib (e.g., navigation) are temporarily 
          disabled to allow movement control.
        """
                    
        # Remove default left / right arrow key press events
        plt.rcParams['keymap.back'].remove('left')
        plt.rcParams['keymap.forward'].remove('right')
        
        if (self.parent.verbose==1 or self.parent.verbose==2):
            instructions = dedent(
            """
            CenterDetermination :: ThreePoints (interactive adjustment)
            Use these keys:
              - '←' : move left
              - '→' : move right
              - '↑' : move up
              - '↓' : move down
              - '+' : increase circle radius
              - '-' : decrease circle radius
              - 'b' : increase step size
              - 'l' : decrease step size
              - 'd' : refinement done
                  
            DISCLAIMER: For the purpose of the center position adjustment, 
                        the default shortcuts for arrows were removed.
            """)
            print(instructions)
        
        if self.parent.print_sums:
            print("Intensity sums during refinement:")
            
        # Initialize variables and flags
        self.backip = np.array((self.x, self.y))
        xy = np.array((self.x, self.y))
        r = np.copy(self.r)
        termination_flag = False
        
        plt.title("Manually adjust the center position.", fontsize=20)

        plt.ion()
          
        ### Define the event handler for figure close event
        def onclose(event):
            nonlocal termination_flag
            termination_flag = True

        # Connect the event handler to the figure close event
        fig.canvas.mpl_connect('close_event', onclose)
        
        # Define the callback function for key press events
        def onkeypress2(event):
            # Use nonlocal to modify the center position in the outer scope
            nonlocal xy, r, termination_flag
        
            # OTHER KEYS USED IN INTERACTIVE FIGURES
            #   event.key == '1': select a point in self.detection_3points()
            #   event.key == '2': delete the last point in self.detection...
            #   event.key == '3': delete a point in self.detection...
            #   event.key == '+': increase circle radius
            #   event.key == '-': decrease circle radius           
            #   event.key == 'b': increase the step size (big step size)
            #   event.key == 'l': decrease the step size (little step size)
            #   event.key == 'd': proceed in self.detection_3points()
        
            if event.key in ['up', 'down', 'left', 'right', '+', '-']:                    
                if event.key in ['+', '-']:
                    r += 1 if event.key == '+' else -1
                else:
                    # Perform shifts normally
                    if event.key == 'up':
                        xy[1] -= self.step
                    elif event.key == 'down':
                        xy[1] += self.step
                    elif event.key == 'left':
                        xy[0] -= self.step
                    elif event.key == 'right':
                        xy[0] += self.step
                    
                    # Print sum only for arrow keys
                    if self.parent.print_sums:
                        s = self.parent.intensity_sum(self.parent.to_refine, 
                                                      xy[0], xy[1], r)
                        print(f'{s:.2f}')
            
            # Terminate the interactive refinement with 'd' key
            if event.key == 'd':
                termination_flag = True
        
            # Change step size 
            if event.key == 'b':
                self.step = self.step * 5
        
            if event.key == 'l':
                self.step = self.step / 5
                if self.step < 0.5:
                    self.step = 0.5
        
            # Update the plot with the new center position
            circle.set_center((xy[0], xy[1]))  # circle
            circle.set_radius(r)               # radius
            center.set_data([xy[0]], [xy[1]])  # center
        
            plt.title("Manually adjust the center position.", fontsize=20)
        
            # Update the plot
            plt.draw()

                
        # Disconnect the on_key_press1 event handler from the figure
        fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
        
        # Connect the callback function to the key press event
        fig.canvas.mpl_connect('key_press_event', onkeypress2)

        # Enable interaction mode
        plt.ion() 
               
        # Wait for 'd' key press or figure closure
        while not termination_flag:
            try:
                plt.waitforbuttonpress(timeout=0.1)
            except KeyboardInterrupt:
                # If the user manually closes the figure, terminate the loop
                termination_flag = True
                
        # Turn off interactive mode
        plt.ioff()
        
        # Display the final figure with the selected center position and radius
        plt.tight_layout()

        plt.show(block=False)
        
        # If the termination_flag is True, stop the code
        if termination_flag: 
            plt.close()  # Close the figure

        # User information:
        if (self.parent.verbose or self.parent.final_print):
            self.parent.dText = "Center Determination (ThreePoints)    : ({:.3f}, {:.3f})"
        
    
        return xy[0], xy[1], r
        

    def calculate_circle(self, plot_results:int)->tuple[
            float,float,float,tuple[float,float], plt.Circle]:
        """
        Calculate the center and radius of a circle defined by 3 manually 
        selected points.

        Given three user-defined points (typically on a diffraction ring), this 
        method computes the center and radius of the circle that passes through 
        them. Optionally, it visualizes the result including the circle, center, 
        and selected points overlaid on the image.

        Parameters
        ----------
        plot_results : int (binary: 0 or 1)
            If 1, displays the image with the fitted circle and center.
            If 0, skips visualization.

        Returns
        -------
        self.x : float
            x-coordinate of the detected circle center.

        self.y : float
            y-coordinate of the detected circle center.

        self.r : float
            Radius of the detected circle.

        self.center : tuple[float, float]
            The (x, y) coordinates of the circle center.

        self.circle : matplotlib.patches.Circle
            A matplotlib Circle object that can be reused for future plotting.

        Notes
        -----
        - The method uses a geometric approach to calculate the circumcenter
          and circumradius.
        - For best visual results, ensure the selected points are well spaced 
          and belong to the same circular ring.
        """
        # Extract the coordinates of the points        
        x = [self.coords[0][0], self.coords[1][0], self.coords[2][0]]
        y = [self.coords[0][1], self.coords[1][1], self.coords[2][1]]
        
        # Compute the radius and center coordinates of the circle
            # a: the squared length of the side between the second 
            #    and third points (x[1], y[1]) and (x[2], y[2]).
            # b: the squared length of the side between the first 
            #    and third points (x[0], y[0]) and (x[2], y[2]).
            # c: the squared length of the side between the first 
            #    and second points (x[0], y[0]) and (x[1], y[1]).
            # s: the twice the signed area of the triangle formed by 3 points
            
        c = (x[0]-x[1])**2 + (y[0]-y[1])**2
        a = (x[1]-x[2])**2 + (y[1]-y[2])**2
        b = (x[2]-x[0])**2 + (y[2]-y[0])**2
        s = 2*(a*b + b*c + c*a) - (a*a + b*b + c*c) 
        
        # coordinates of the center
        self.x = (a*(b+c-a)*x[0] + b*(c+a-b)*x[1] + c*(a+b-c)*x[2]) / s
        self.y = (a*(b+c-a)*y[0] + b*(c+a-b)*y[1] + c*(a+b-c)*y[2]) / s 
        
        # radius
        ar = a**0.5
        br = b**0.5
        cr = c**0.5 
        self.r = ar*br*cr/((ar+br+cr)*(-ar+br+cr)*(ar-br+cr)*(ar+br-cr))**0.5
        
        # # Print results
        # if self.parent.verbose:
        #     print("CenterEstimator :: manual center detection")
        #     print(f"Center coordinates: {self.x:.2f} {self.y:.2f}")
                    
        if plot_results==1:
            # Create and manage the figure
            fig, ax = plt.subplots()
            manager = plt.get_current_fig_manager()
            manager.window.showMaximized()
            ax.imshow(self.parent.image, cmap = self.parent.cmap,
                      origin="upper")
            
            # Plot center and points
            center, = plt.plot(self.x, self.y, 
                     'rx', 
                     label='Center', 
                     markersize=12)
            plt.scatter(x,y, 
                        marker='x', 
                        color='palevioletred', 
                        label = 'Circle points')
            plt.title('Circle found using 3 manually detected points')
            
            # Circle visualization
            circle = plt.Circle((self.x,self.y), 
                                self.r, 
                                color='palevioletred', 
                                fill=False,
                                label = 'pattern')
            ax.add_artist(circle)
            
            # Set the aspect ratio to equal to have a circular shape
            plt.axis('equal')
            
            plt.legend(loc='lower center', 
                       ncol=2, 
                       bbox_to_anchor=(0.5,-0.1), 
                       mode='expand', 
                       frameon=False)
            plt.axis('off')
            plt.tight_layout()
            plt.show(block=False)

        
        self.center = (self.x, self.y)
        self.circle = plt.Circle((self.x,self.y),self.r)
        

        return self.x, self.y, self.r, self.center, self.circle


    def get_radius(
            self, rtype:int, im:np.ndarray, x:float, y:float, disp:bool=False) -> float:
        """
        Calculate the radius of a circle based on intensity profiles along 
        horizontal and vertical axes.
    
        Parameters
        ----------
        rtype : integer
            Method for radius calculation. Default is 0.
            - rtype=0 : peaks matching method
            - rtype=1 : radial distribution method
            
        im : np.ndarray
            The 2D image array containing the circle. 
            
        x : float
            The x-coordinate of the circle's center.
            
        y : float
            The y-coordinate of the circle's center.
            
        disp : bool, optional
            If True, visualizes the detected intensity profiles and peaks 
            (default is False).
    
        Returns
        -------
        float
            The estimated radius of the circle. Defaults to 100 if no valid 
            radius is detected.
        """
        # Helper function -----------------------------------------------------
        def match_peaks(arr):
            """
            Find the most similar values across the left and right halves of 
            an array, excluding the global maximum.
        
            This function is typically used for analyzing symmetric patterns 
            such as diffraction profiles or intensity histograms. It identifies 
            the most similar pair of peaks on either side of the central maximum, 
            which is excluded from the comparison.
        
            Special case:
            - If the array contains exactly two elements, it assumes they form 
              the best pair and returns their indices directly.
            - If one half is empty after removing the central peak,
              it returns (None, None).
        
            Parameters
            ----------
            arr : np.ndarray or list of float
                1D array of values (intensity profile) where symmetric peak
                positions are to be matched.
        
            Returns
            -------
            best_pair : tuple[int or None, int or None]
                Indices of the two most similar values, one from each half of 
                the array. Returns (None, None) if no valid pair is found.
        
            Notes
            -----
            - The function assumes the most prominent peak (maximum value) lies 
              at or near the center and represents the central feature.
            - The left and right halves are compared pairwise, and the closest 
              match in absolute value difference is returned.
            - This is useful for refining radial symmetry or finding symmetric 
              ring peaks.
            """
            if len(arr) < 2:  # Not enough values to compare
                if self.parent.verbose==2:
                    print("Not enough values to find similar pairs.")
                return None, None
    
            # Special case: if there are exactly 2 peaks, return them
            if len(arr) == 2:
                if self.parent.verbose==2:
                    print("Exactly two peaks detected. Returning them as the best pair.")
                return 0, 1
    
            # Find the index of the highest value
            center_idx = np.argmax(arr)
    
            # Split into left and right halves, excluding the highest value
            left = arr[:center_idx]
            right = arr[center_idx + 1:]
            left_indices = np.arange(center_idx)
            right_indices = np.arange(center_idx + 1, len(arr))
    
            if len(left) == 0 or len(right) == 0:  # Check for empty halves
                if self.parent.verbose==2:
                    print("One of the halves is empty.")
                return None, None
    
            # Initialize variables for tracking the smallest difference
            smallest_diff = np.inf
            best_pair = (None, None)
    
            # Compare each value in the left with every value in the right
            for i, l_val in enumerate(left):
                for j, r_val in enumerate(right):
                    diff = abs(l_val - r_val)
                    if diff < smallest_diff:
                        smallest_diff = diff
                        best_pair = (left_indices[i], right_indices[j])
    
            return best_pair
        
        # Peak matching method ------------------------------------------------
        if rtype == 0:
            self.xpeaks, self.ypeaks = None, None
            self.xyvals, self.yyvals = None, None
        
            x_line = im[int(x), :]
            y_line = im[:, int(y)]
        
            # Define threshold for peak detection
            x_thr = 0.5 * max(x_line)
            y_thr = 0.5 * max(y_line)
        
            # Find peaks with dynamic height thresholds
            self.xpeaks, _ = find_peaks(x_line, 
                                        height=x_thr, 
                                        prominence=1, 
                                        distance=30)
            self.xyvals = x_line[self.xpeaks]
            self.ypeaks, _ = find_peaks(y_line, 
                                        height=y_thr, 
                                        prominence=1, 
                                        distance=30)
            self.yyvals = y_line[self.ypeaks]
        
            # Define half the length of the image
            half_length_x = x_line.shape[0] / 2
            half_length_y = y_line.shape[0] / 2
        
            # Check the additional condition for xpeaks
            if len(self.xpeaks) == 2 and (
                (self.xpeaks[0]<half_length_x and self.xpeaks[1]<half_length_x) or
                (self.xpeaks[0]>half_length_x and self.xpeaks[1]>half_length_x)):
                if self.parent.verbose==2:
                    print("xpeaks condition met: Both peaks are on the same side of the center.")
                self.pairX = None
            else:
                self.pairX = match_peaks(self.xyvals)
        
            # Check the additional condition for ypeaks
            if len(self.ypeaks) == 2 and (
                (self.ypeaks[0]<half_length_y and self.ypeaks[1]<half_length_y) or
                (self.ypeaks[0]>half_length_y and self.ypeaks[1]>half_length_y)):
                if self.parent.verbose==2:
                    print("ypeaks condition met: Both peaks are on the same side of the center.")
                self.pairY = None
            else:
                self.pairY = match_peaks(self.yyvals)
        
            # Determine radius based on available pairs
            if self.pairX is None or None in self.pairX:
                rx_x = None
            else:
                x1 = self.xpeaks[self.pairX[0]]
                x2 = self.xpeaks[self.pairX[1]]
                rx_x = abs(x1 - x2) / 2
        
            if self.pairY is None or None in self.pairY:
                rx_y = None
            else:
                y1 = self.ypeaks[self.pairY[0]]
                y2 = self.ypeaks[self.pairY[1]]
                rx_y = abs(y1 - y2) / 2
        
            if rx_x is not None and rx_y is not None:
                rx = np.mean([rx_x, rx_y])
            elif rx_x is not None:
                rx = rx_x
            elif rx_y is not None:
                rx = rx_y
            else:
                if self.parent.verbose==2:
                    print("No valid pairs detected for radius calculation.")
                return 100  # Default radius or error handling
        
            if disp:
                # Plot xline with peaks
                plt.figure(figsize=(12, 6))
        
                # Plot for xline
                plt.subplot(2, 1, 1)
                plt.plot(x_line, label='xline')
                plt.plot(self.xpeaks, self.xyvals, "ro", label='Peaks')
                plt.title('Peaks in xline')
                plt.xlabel('Index')
                plt.ylabel('Value')
                plt.legend()
        
                # Plot for yline
                plt.subplot(2, 1, 2)
                plt.plot(y_line, label='yline')
                plt.plot(self.ypeaks, self.yyvals, "ro", label='Peaks')  
                plt.title('Peaks in yline')
                plt.xlabel('Index')
                plt.ylabel('Value')
                plt.legend()
        
                plt.tight_layout()
                plt.show()
        
        # Radial distribution method ------------------------------------------
        elif rtype == 1:
            profile = ediff.radial.calc_radial_distribution(
                im, 
                center=(x, y)
                )
            p0, p1 = profile
            idx_max = np.argmax(p1[50:]) + 50
            rx = p0[idx_max]
            max_val = p1[idx_max]
    
            if disp:
                fig, axs = plt.subplots(1, 2, figsize=(8, 3))
                axs[0].imshow(im, cmap=self.parent.cmap)
                axs[0].add_patch(
                    Circle(
                        (y, x), rx, 
                        edgecolor='red', 
                        fill=False, 
                        linewidth=2))
                axs[0].scatter(y, x, color="red", marker="x")
                axs[0].set_title("Initial Estimate")
                axs[0].axis('off')
                axs[1].plot(p0, p1)
                axs[1].set_ylim(0,max_val+20)
                axs[1].axvline(rx, 
                               color='red', 
                               linestyle='--', 
                               label='Initial Radius')
                axs[1].set_title("Radial Profile")
                axs[1].legend()
                plt.tight_layout()
                plt.show()
        
        # Return found radius
        return rx


    def backend_fft(f):
        """
        Decorator to enforce the use of PocketFFTBackend when utilizing the 
        `scipy.fft` module. This ensures that all FFT computations leverage 
        the optimized PocketFFT implementation, improving performance, 
        particularly for multi-threaded execution.
    
        Parameters
        ----------
        f : callable
            The function to be wrapped, ensuring it executes within the 
            PocketFFTBackend context.
    
        Returns
        -------
        callable
            A wrapped function that forces the use of PocketFFTBackend.
        """
        @wraps(f)
        def newf(*args, **kwargs):
            with set_backend(PocketFFTBackend):
                return f(*args, **kwargs)
        return newf
    
    @backend_fft
    def calc_fft(self, fft_function):
        """
        Computes the Fast Fourier Transform (FFT) using the specified function.

        This method ensures that the FFT operation is executed using 
        the PocketFFTBackend, applying optimized settings for improved efficiency.

        Parameters
        ----------
        fft_function : callable
            The FFT function to be applied to the data.

        Returns
        -------
        object
            The result of the FFT computation.
        """
        return fft_function

   
    def auto_masking(self, im, threshold=0.2, show_mask=False):
        """
        Automatically generates a binary mask for an image based on intensity 
        thresholding.
        
        This function determines a lower intensity limit by taking a fraction 
        (threshold) of the median of the highest intensity values in the image. 
        This helps filter out noise while preserving meaningful features, 
        avoiding the influence of hot spots.
    
        Parameters
        ----------
        im : ndarray
            Input grayscale image as a 2D NumPy array.
        
        threshold : float, optional
            A fraction (default is 0.1) of the median of the highest intensity 
            values. Pixels with intensity above this computed limit will be 
            included in the mask.
    
        Returns
        -------
        mask : ndarray
            A boolean mask (same shape as `im`), where `True` indicates pixels 
            above the threshold and `False` represents background or noise.
    
        """
        # Find the median of the highest intensity value of the image to avoid 
        # hot spots (More stable than median scaling)
        lower_limit = np.percentile(im, threshold * 100) 
        
        # Generate a mask
        mask = im > lower_limit
        
        if show_mask:
            plt.figure()
            plt.imshow(mask)
            plt.title("Beam Stopper Masking")
            plt.tight_layout()
            plt.axis("off")
            plt.show()
        
        return mask

        
    def visualize_center(self, x: float, y: float, r: float, rd: str) -> None:
        '''         
        Visualize detected diffraction patterns and mark the center.
        
        Parameters
        ----------
        tit : string
            name of the method used for circle detection
            
        x : float64
            x-coordinate of the detected center
            
        y : float64
            y-coordinate of the detected center
            
        r : float64
            radius of the detected center
        
        Returns
        -------
        None.
                            
        '''
        # Load image
        image = np.copy(self.parent.to_refine)
    
        if rd=="d":
            tit = "Detected center position"
        elif rd=="r":
            tit = "Refined center position"
        else: 
            tit = "Center position"
            
            
        if self.parent.icut is not None:
            im = np.where(image > self.parent.icut, self.parent.icut, image)
        else:
            im = np.copy(image)
            
        # Create a figure and display the image
        fig, ax = plt.subplots()
        
        # Allow using arrows to move back and forth between view ports
        plt.rcParams['keymap.back'].append('left')
        plt.rcParams['keymap.forward'].append('right')
 
        plt.title(tit)
        ax.axis('off')

        # Plot center point
        ax.scatter(x,y,
                label= f'center:  [{x:.1f}, {y:.1f}]',
                marker='x', color="red", s=60)

        plt.legend(loc='upper right')
        
        # Display the image
        ax.imshow(im, cmap = self.parent.cmap, origin="upper")
        plt.axis('off')
        plt.tight_layout()
        plt.show(block=False)



class CenterRefinement:
    '''
    CenterRefinement - Final refinement of the center of a diffraction pattern.

    This class performs refinement of the detected center coordinates in a 2D
    diffraction pattern. It supports multiple refinement methods, both manual
    and automatic, and optionally saves the refined center to a file.


    Parameters
    ----------
    parent : CenterLocator
        Reference to the parent CenterLocator object. Used to access shared
        attributes like initial center estimates and configuration flags.
        
    input_image : str, Path, or numpy.ndarray
        The diffraction pattern image, either as a file path or a 2D NumPy array.
    
        
    refinement : str or None, optional, default=None
        The refinement method to apply. Options include:
        - 'manual': Manual adjustment via user interaction on the selected ring.
        - 'sum'   : Automatic refinement by maximizing the intensity sum 
        - 'var'   : Automatic refinement by minimizing the intensity variance
        If None, no refinement is applied.
          
    out_file : str, optional
        Path to a text file where the refined center coordinates will be saved.
    
    heq : bool, optional, default=False
        Whether to apply histogram equalization to the input image before
        processing.
    
    icut : float, optional, default=None
        Cut-off intensity level for processing the image
        
   cmap : str, optional, default='gray'
        Colormap to be used when displaying the image for manual refinement.
    
    verbose : bool, optional, default=False
        Flag to enable or disable informational verbose during processing. 

    print_sums : bool, optional, default=False
        If True, prints the sum of intensity values for the refined circle 
        after each adjustment, relevant only for manual refinement methods.

    Returns
    -------
    None

    Notes
    -----
    - If an unsupported refinement method is provided, the program exits with
      an error.
    - The coordinates `xx`, `yy`, and radius `rr` can be accessed after
      initialization for further analysis or saving.
    - The refinement strategy relies on initial estimates provided by the parent 
      object.
    - Manual refinement supports click-based center adjustments guided by visual 
      feedback.
    - Sum/variance refinement automatically searches for the most likely center
      by evaluating ring statistics.
    '''
    
    def __init__(self, parent,
                 input_image, 
                 refinement = None,
                 in_file = None,
                 out_file = None,
                 heq = False, 
                 icut = None,
                 cmap = 'gray',
                 verbose = 0,
                 print_sums = False):
        
        ######################################################################
        # PRIVATE FUNCTION: Initialize CenterLocator object.
        # The parameters are described above in class definition.
        ######################################################################
        
        ## (0) Initialize input attributes ------------------------------------
        self.parent = parent
        
        ## (1) Initialize new attributes --------------------------------------
        self.step=0.5
        
        ## (2) Run functions --------------------------------------------------
        if refinement is not None:
            # Flag for later
            self.ret = 1
            par_short = self.parent.center1
        
            if refinement == "manual":
                # Manual refinement method
                if parent.determination == 'manual':
                    self.xx, self.yy, self.rr = \
                        par_short.x, par_short.y, par_short.r
                    par_short.x, par_short.y = \
                        par_short.backip[0], par_short.backip[1]
                    if (self.parent.verbose or self.parent.final_print):
                        self.parent.rText = \
                            "Center Refinement (Interactive)       : ({:.3f}, {:.3f})"
                elif parent.determination == 'intensity':
                    self.yy, self.xx, self.rr = self.ref_interactive(
                        par_short.y, par_short.x, par_short.r)
                else:
                    self.yy, self.xx, self.rr = self.ref_interactive(
                        par_short.x, par_short.y, par_short.r)
        
            elif refinement == "var":
                # Intensity variance refinement
                self.xx, self.yy, self.rr = self.ref_var(
                    par_short.x, par_short.y, par_short.r)
        
            elif refinement == "sum":
                # Intensity sum refinement
                self.xx, self.yy, self.rr = self.ref_sum(
                    par_short.x, par_short.y, par_short.r,
                    live_plot=self.parent.live_plot)
            

            else:
                print("Selected refinement method is not supported.")
                sys.exit()
       
            plt.close("all")          
        else: 
            # Flag for later
            self.ret = 2
            
           
    def ref_interactive(self, px, py, pr):
        """
        Interactive refinement of the diffraction pattern center after initial
        detection via 3 points.
      
        This method allows the user to fine-tune the estimated center and radius
        of a diffraction ring (initially determined using 3 manually selected 
        points) through interactive keyboard controls.
      
        The user interacts directly with the figure window containing the image 
        and the detected circle.
      
        ### Keyboard Controls
        ---------------------
        - Arrow keys (←, →, ↑, ↓): Move the center left, right, up, or down
        - '+' : Increase the radius
        - '-' : Decrease the radius
        - 'b' : Increase step size (×5)
        - 'l' : Decrease step size (/5, with minimum step size 0.5)
        - 'd' : Done — finalize the center and radius
        - Closing the figure: Cancels refinement and returns the original input 
          center and radius
      
        Parameters
        ----------
        fig : matplotlib.figure.Figure
            The figure window used for interactive refinement.
            
        circle : matplotlib.patches.Circle
            The circle object initially defined from 3 selected points.
            
        center : tuple
            The initial (x, y) center coordinates estimated from the selected 
            points.
            
        plot_results : bool, optional, default=False
            If True, the results of the adjustment are visualized.
      
        Returns
        -------
        xy[0] : float
            x-coordinate of the adjusted center of the diffraction pattern.
            
        xy[1] : float
            y-coordinate of the adjusted center of the diffraction pattern.
            
        r : float
            Adjusted radius of the diffraction ring.
      
        Notes
        -----
        - Interactive adjustments are visualized in real time.
        - Intensity sum at the current center/radius is printed if print_sums
          is True.
        - Arrow key defaults in Matplotlib (e.g., navigation) are temporarily 
          disabled to allow movement control.
        """
             
        # (0) Load original image ---------------------------------------------
        im = np.copy(self.parent.to_refine)

        # Edit contrast with a user-predefined parameter
        if self.parent.icut is not None:
            if self.parent.verbose==2:
                print("Contrast enhanced.")
            im = np.where(im > self.parent.icut, 
                              self.parent.icut, 
                              im)
            
        # (1) Initialize variables and flags ----------------------------------
        xy = np.array((px, py))
        r = np.copy(pr)
        termination_flag = False

        # (2) User information ------------------------------------------------
        if (self.parent.verbose==1 or self.parent.verbose==2):
            instructions = dedent("""
            
            Interactive refinement. Use these keys:
                  - '←' : move left
                  - '→' : move right
                  - '↑' : move up
                  - '↓' : move down
                  - '+' : increase circle radius
                  - '-' : decrease circle radius
                  - 'b' : increase step size
                  - 'l' : decrease step size
                  - 'd' : refinement done
                  
            DISCLAIMER: For the purpose of the center shift, the default
            shortcuts for left and right arrows were removed.
            """)
            print(instructions)
        
        if self.parent.print_sums:
            print("Intensity sums during refinement:")
            
        # (3) Create a figure and display the image ---------------------------
        fig, ax = plt.subplots(figsize=(12, 12))
        
        # Allow using arrows to move back and forth between view ports
        plt.rcParams['keymap.back'].append('left')
        plt.rcParams['keymap.forward'].append('right')
        
        circle = plt.Circle(
            (px, py), pr, color='r', fill=False)
        ax.add_artist(circle)

        # Plot center point
        center, = ax.plot(px, py, 'rx', markersize=12)
                    

        plt.title('Manually adjust the center position.', 
                  fontsize=20)

        ax.imshow(im, cmap = self.parent.cmap, origin="upper")
        ax.axis('off')
        
        # (4) Enable interactive mode -----------------------------------------
        plt.ion()
        

        # (5) Display the image -----------------------------------------------
        # fig.set_size_inches(self.fig_width, self.fig_height)
        plt.show(block=False)
        
        # (6) Define the event handler for figure close event -----------------
        def onclose(event):
            nonlocal termination_flag
            termination_flag = True
 
        # Connect the event handler to the figure close event
        fig.canvas.mpl_connect('close_event', onclose)

        # (7) Define the callback function for key press events ---------------
        def onkeypress2(event):
            # Use nonlocal to modify the center position in the outer scope
            nonlocal xy, r, termination_flag
        
            # OTHER KEYS USED IN INTERACTIVE FIGURES
            #   event.key == '1': select a point in self.detection_3points()
            #   event.key == '2': delete the last point in self.detection...
            #   event.key == '3': delete a point in self.detection...
            #   event.key == '+': increase circle radius
            #   event.key == '-': decrease circle radius           
            #   event.key == 'b': increase the step size (big step size)
            #   event.key == 'l': decrease the step size (little step size)
            #   event.key == 'd': proceed in self.detection_3points()
        
            if event.key in ['up', 'down', 'left', 'right', '+', '-']:                    
                if event.key in ['+', '-']:
                    r += 1 if event.key == '+' else -1
                else:
                    # Perform shifts normally
                    if event.key == 'up':
                        xy[1] -= self.step
                    elif event.key == 'down':
                        xy[1] += self.step
                    elif event.key == 'left':
                        xy[0] -= self.step
                    elif event.key == 'right':
                        xy[0] += self.step
                    
                    # Print sum only for arrow keys
                    if self.parent.print_sums:
                        s = self.parent.intensity_sum(self.parent.to_refine, 
                                                      xy[0], xy[1], r)
                        print(f'{s:.2f}')
            
            # Terminate the interactive refinement with 'd' key
            if event.key == 'd':
                termination_flag = True
        
            # Change step size 
            if event.key == 'b':
                self.step = self.step * 5
        
            if event.key == 'l':
                self.step = self.step / 5
                if self.step < 0.5:
                    self.step = 0.5
        
            # Update the plot with the new center position
            circle.set_center((xy[0], xy[1]))  # circle
            circle.set_radius(r)               # radius
            center.set_data([xy[0]], [xy[1]])  # center
        
            plt.title("Manually adjust the center position.", fontsize=20)
        
            # Update the plot
            plt.draw()

            
        # Disconnect the on_key_press1 event handler from the figure
        fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
        
        # Connect the callback function to the key press event
        fig.canvas.mpl_connect('key_press_event', onkeypress2)

        # Enable interaction mode
        plt.ion() 
               
        # Wait for 'd' key press or figure closure
        while not termination_flag:
            try:
                plt.waitforbuttonpress(timeout=0.1)
            except KeyboardInterrupt:
                # If the user manually closes the figure, terminate the loop
                termination_flag = True
         
        # (8) Turn off interactive mode ---------------------------------------
        plt.ioff()
        
        # Display the final figure with the selected center position and radius
        plt.tight_layout()

        plt.show(block=False)
        
        # If the termination_flag is True, stop the code
        if termination_flag: 
            plt.close()  # Close the figure

        # User information:
        if (self.parent.verbose or self.parent.final_print):
            self.parent.rText = "Center Refinement (Interactive)       : ({:.3f}, {:.3f})"

        return xy[0], xy[1], r    
        
    
    def ref_var(self, px, py, pr, plot_results=0):
        """
        Refine the center coordinates (px, py) and radius (pr) of a circular 
        diffraction pattern by minimizing the variance of pixel intensities 
        along the circle's border.

        This method performs iterative local search in the 8-neighbourhood 
        around the current center to find a position and radius that minimize 
        the intensity variance. The refinement stops when convergence criteria 
        are met or after a maximum number of iterations.

        Neighborhood pattern tested (o = current center, x = candidate center):
            x x x   --> (px - dx, py + dy) (px, py + dy) (px + dx, py + dy)
            x o x   --> (px - dx, py)      (px, py)      (px + dx, py)
            x x x   --> (px - dx, py - dy) (px, py - dy) (px + dx, py - dy)

        Parameters
        ----------
        px : float
            Initial x-coordinate of the detected center.
            
        py : float
            Initial y-coordinate of the detected center.
            
        pr : float
            Initial radius of the detected circular pattern.
            
        plot_results : int, optional (default=0)
            Whether to plot the refinement result (1 = yes, 0 = no).

        Returns
        -------
        px : float
            Refined x-coordinate of the center.
            
        py : float
            Refined y-coordinate of the center.
            
        pr : float
            Refined radius of the circular pattern.

        Notes
        -----
        - The function avoids overfitting to local minima by enforcing 
          convergence thresholds and consistency checks.
        - If the refinement worsens the initial variance significantly or 
          suggests an invalid coordinate swap, the original input is preserved.
        - Uses `self.parent.intensity_variance()` for evaluating the criterion.
        """
        
        # Store input for plot
        bckup = [np.copy(px), np.copy(py), np.copy(pr)]
        
        # Load image
        image = np.copy(self.parent.image)

        # Starting values to be modified 
        init_var = self.parent.intensity_variance(image, px, py, pr)
        min_intensity_var = self.parent.intensity_variance(image, px, py, pr)
        best_center = (np.copy(px), np.copy(py))
        best_radius = np.copy(pr)
    
        # Convergence criterion for termination of gradient optimization 
        # (1) small positive value that serves as a threshold to determine 
        #     when the optimization process has converged
        convergence_threshold = 0.1*min_intensity_var
        
        # (2) maximum number of iterations of optimization
        max_iterations = 10
        
        # (3) keep track of the number of consecutive iterations where there 
        #     is no improvement in the objective function beyond 
        #     the convergence threshold
        no_improvement_count = 0
    
    
        # iterative refinement of the center of a circle while keeping
        # the radius constant.
        step = 0.3
        neighbors = [(float(dx), float(dy))
            for dx in np.arange(-1, 1 + step, step)
            for dy in np.arange(-1, 1 + step, step)]
        
        for iteration in range(max_iterations):    
            # Refine center while keeping radius constant
            curr = self.parent.intensity_variance(image, 
                                      best_center[0], 
                                      best_center[1], 
                                      best_radius) 
            # Store intensity sums of the current center's neighborhood
            curr_intensity_var = []
            for dx, dy in neighbors:
                nx, ny = best_center[0] + dx, best_center[1] + dy
                # Check if the point is within the expanded search radius
                curr_intensity_var.append(
                    self.parent.intensity_variance(image, nx, ny, best_radius))
            
            # Find the minimum value coordinates within curr_intensity_var
            cx, _ = np.unravel_index(np.argmin(curr_intensity_var),
                                     [len(curr_intensity_var),1])
                    
            # Check for improvement of criterion -- in each iteration just once,
            # as the algorithm checks the neighbourhood of the best center (in
            # each iteration, the center is updated if possible)
            if min(curr_intensity_var) <= min_intensity_var:                           
                min_intensity_var = max(curr_intensity_var)
                
                # Calculate the new best coordinates of the center
                n = neighbors[cx]
                (nx, ny) = tuple(map(lambda x, y: float(x) + float(y), 
                                     best_center, n))
                best_center = px, py = (np.copy(nx), np.copy(ny))
                
            # Update maximum intensity sum 
            min_intensity_var = self.parent.intensity_variance(image, 
                                                    best_center[0], 
                                                    best_center[1], 
                                                    best_radius) 
            
            # Refine radius if necessary while keeping the center position 
            # constant. It iterates through different radius adjustments
            # to find a radius that maximizes the intensity sum of pixels.
            
            radi_intensity_var = []
            radii = np.arange(-1, 1 + step, step)
            for dr in radii:
                new_radius = best_radius + dr
                radi_intensity_var.append(self.parent.intensity_variance(image, 
                                                             best_center[0], 
                                                             best_center[1], 
                                                             new_radius))
                
            # Find the minimum value coordinates within curr_var
            rx, _ = np.unravel_index(np.argmin(radi_intensity_var),
                                      [len(radi_intensity_var),1])
            
            # Check for improvement of criterion
            if max(radi_intensity_var) < min_intensity_var:
                min_intensity_var = max(radi_intensity_var)
                
                n = radii[rx]
                nr = best_radius+n
                
                best_radius = pr = np.copy(nr)

            
            # Check for convergence and improvement (termination conditions)
            impr = abs(min_intensity_var - curr)
            if impr < convergence_threshold:
                no_improvement_count += 1
                if no_improvement_count == 5:
                    break
        
        # Avoid incorrect/redundant refinement
        ## (1) swapped coordinates
        if ((bckup[0] > bckup[1] and not best_center[0] > best_center[1])
            or  (bckup[0] < bckup[1] and not best_center[0] < best_center[1])):
            best_center = best_center[::-1]
        
        ## (2) worsened final maximum intensity sum than the initial one
        if np.round(init_var,-2) < np.round(min_intensity_var,-2):
            print("Refinement redundant.")
            best_center = np.copy(bckup)
    
        # Print results
        if (self.parent.verbose or self.parent.final_print):
            self.parent.rText = \
                "Center Refinement (IntensityVar)      : ({:.3f}, {:.3f})"
                                  
        return best_center[0], best_center[1], best_radius
    
    
    def ref_sum(self, px, py, pr, plot_results=0, live_plot=False):
        """
        Refine the center coordinates (px, py) and radius (pr) of a circular 
        diffraction pattern by maximizing the summed pixel intensity along 
        a circular ring.

        The method uses an iterative local search strategy based on gradient 
        ascent in the 8-neighbourhood around the current center to find 
        a position and radius that maximize the intensity sum. The center 
        is updated based on the best neighbor until convergence criteria 
        or maximum iteration limits are met.

        Neighborhood pattern tested (o = current center, x = candidate center):
            x x x   --> (px - dx, py + dy) (px, py + dy) (px + dx, py + dy)
            x o x   --> (px - dx, py)      (px, py)      (px + dx, py)
            x x x   --> (px - dx, py - dy) (px, py - dy) (px + dx, py - dy)

        Parameters
        ----------
        px : float
            Initial x-coordinate of the detected center.
            
        py : float
            Initial y-coordinate of the detected center.
            
        pr : float
            Initial radius of the detected circular pattern.
            
        plot_results : int, optional (default=0)
            If set to 1, plots the final result of the refinement.

        Returns
        -------
        px : float
            Refined x-coordinate of the center.
            
        py : float
            Refined y-coordinate of the center.
            
        pr : float
            Refined radius of the circular pattern.

        Notes
        -----
        - Uses `self.parent.intensity_sum()` as the optimization criterion.
        - Alternates between optimizing the center and the radius.
        - Stops after no significant improvement over multiple iterations.
        - If refinement worsens the result, the original values are retained.
        """
        # Store input for plot via self.visualize_refinement()
        bckup = [np.copy(px), np.copy(py), np.copy(pr)]
    
        # Image in which the center is refined
        image = np.copy(self.parent.image)
    
        # Starting values to be modified 
        init_sum = self.parent.intensity_sum(image, px, py, pr)
        max_intensity_sum = self.parent.intensity_sum(image, px, py, pr)
        best_center = (np.copy(px), np.copy(py))
        best_radius = np.copy(pr)
        
        # Convergence criterion for termination of gradient optimization 
        # (1) small positive value that serves as a threshold to determine 
        #     when the optimization process has converged
        convergence_threshold = 0.05*max_intensity_sum
        
        # (2) maximum number of iterations of optimization
        max_iterations = 100
        
        # (3) keep track of the number of consecutive iterations where there 
        #     is no improvement in the objective function beyond 
        #     the convergence threshold
        no_improvement_count = 0
         
        # iterative refinement of the center of a circle while keeping
        # the radius constant.
        step = 0.2
        neighbors = [(float(dx), float(dy))
            for dx in np.arange(-1.0, 1.0 + step, step)
            for dy in np.arange(-1.0, 1.0 + step, step)]
        
        # Live plot initialization 
        if live_plot:
            plt.ion()
            fig, ax = plt.subplots()
            ax.imshow(image, cmap='gray')
            circ_artist = Circle(best_center, best_radius, 
                                 fill=False, color='red', lw=1.5)
            ax.add_patch(circ_artist)
            ax.set_title("Refinement Progress")
            ax.axis("off")
            fig.show()
    
        for iteration in range(max_iterations):    
            # Refine center while keeping radius constant
            curr = self.parent.intensity_sum(image, 
                                      best_center[0], 
                                      best_center[1], 
                                      best_radius)
            
            # Store intensity sums of the current center's neighborhood
            curr_intensity_sum = []
            for dx, dy in neighbors:
                nx, ny = best_center[0] + dx, best_center[1] + dy
                # Check if the point is within the expanded search radius
                curr_intensity_sum.append(self.parent.intensity_sum(image, 
                                                        nx, ny, 
                                                        best_radius))
            
            # Find the maximum value coordinates within curr_sum
            cx, _ = np.unravel_index(np.argmax(curr_intensity_sum),
                                     [len(curr_intensity_sum),1])
                    
            # Check for improvement of criterion - in each iteration just once,
            # as the algorithm checks the neighbourhood of the best center
            # (in each iteration, the center is updated if possible)
            if max(curr_intensity_sum) > max_intensity_sum:                           
                max_intensity_sum = max(curr_intensity_sum)
                
                # Calculate the new best coordinates of the center
                n = neighbors[cx]
                (nx, ny) = tuple(map(lambda x, y: float(x) + float(y), 
                                     best_center, n))
                best_center = px, py = (np.copy(nx), np.copy(ny))
    
            # Update maximum intensity sum 
            max_intensity_sum = self.parent.intensity_sum(image, 
                                                    best_center[0], 
                                                    best_center[1], 
                                                    best_radius)
    
        
            # Refine radius if necessary while keeping the center position 
            # constant. It iterates through different radius adjustments
            # to find a radius that maximizes the intensity sum of pixels.
            
            radi_intensity_sum = []
            radii = np.arange(-1.0, 1.0 + step, step)
            for dr in radii:
                new_radius = best_radius + dr
                radi_intensity_sum.append(self.parent.intensity_sum(image, 
                                                            best_center[0], 
                                                            best_center[1], 
                                                            new_radius))
                
            # Find the maximum value coordinates within curr_sum
            rx, _ = np.unravel_index(np.argmax(radi_intensity_sum),
                                      [len(radi_intensity_sum),1])
    
            # Check for improvement of criterion
            if max(radi_intensity_sum) > max_intensity_sum:
                max_intensity_sum = max(radi_intensity_sum)
                
                n = radii[rx]
                nr = best_radius+n
                
                best_radius = pr = np.copy(nr)
                
            
            # Check for convergence and improvement (termination conditions)
            impr = abs(max_intensity_sum - curr)
            if impr < convergence_threshold:
                no_improvement_count += 1
                if no_improvement_count == 25:
                    break
                
            # --- Live plot update ---
            if live_plot:
                circ_artist.center = best_center
                circ_artist.set_radius(best_radius)
                ax.set_title(
                    f"Iteration {iteration+1} | Sum={max_intensity_sum:.1f}")
                fig.canvas.draw()
                fig.canvas.flush_events()
                plt.pause(0.01)  # small delay for visual update
        
        # Avoid incorrect/redundant refinement
        # ## (1) swapped coordinates
        # if ((bckup[0] > bckup[1] and not best_center[0] > best_center[1])
        #     or  (bckup[0] < bckup[1]
        #     and not best_center[0] < best_center[1])):
        #     best_center = best_center[::-1]
        
        ## (2) worsened final maximum intensity sum than the initial one
        if np.round(init_sum,-2) > np.round(max_intensity_sum,-2):
            print("Refinement redundant.")
            best_center = np.copy(bckup)
    
        # Print results
        if (self.parent.verbose or self.parent.final_print):
            self.parent.rText = \
                "Center Refinement (IntensitySum)      : ({:.3f}, {:.3f})"
        
        if live_plot:
            plt.ioff()
            plt.show()

        return best_center[0], best_center[1], best_radius
    
    
    def diagnose_landscape(self, px, py, pr, metric='std', search_range=2.0, steps=21):
        """
        Visualize the optimization landscape of a center refinement metric
        around a candidate center location and radius.
    
        This function evaluates the selected metric over a grid of positions
        in a 2D neighborhood around the given (px, py) center, using a fixed
        radius `pr`, and visualizes the results as a heatmap with contours.
    
        It can be useful for diagnosing whether the optimization function
        (e.g., standard deviation, median, or sum of pixel intensities) has
        a smooth and well-behaved landscape, which is important for gradient
        optimization methods.
    
        Parameters
        ----------
        px : float
            X-coordinate of the candidate center.
            
        py : float
            Y-coordinate of the candidate center.
            
        pr : float
            Radius of the circular region to evaluate the metric on.
            
        metric : str, optional
            The metric to evaluate. Options are:
            - 'std': standard deviation of pixel values (default)
            - 'median': median of pixel values
            - 'sum': sum of pixel values
            
        search_range : float, optional, default=2.0
            The radius of the square region (in pixels) around the center 
            to explore. The grid will extend ±`search_range` in both x and y 
            directions.
            
        steps : int, optional, default=21
            The number of points to sample along each axis (grid resolution).
            Must be an odd number to ensure the center point is included.
    
        Returns
        -------
        None
            This method produces a matplotlib plot and does not return any value.
        """
        # Helper functions ----------------------------------------------------
        def intensity_sum(self, image, cx, cy, r):
            pixels = self.parent.get_circle_pixels(image, cx, cy, r)
            return np.sum(pixels)
        
        def intensity_std(self, image, cx, cy, r):
            pixels = self.parent.get_circle_pixels(image, cx, cy, r)
            return np.std(pixels)
        
        def intensity_median(self, image, cx, cy, r):
            pixels = self.parent.get_circle_pixels(image, cx, cy, r)
            return np.median(pixels)

        print(f"Diagnosing landscape around ({px:.1f}, {py:.1f}) with r={pr:.1f}...")
        
        # (0) Select the metric function --------------------------------------
        metric_map = {'std': intensity_std, 
                      'median': intensity_median, 
                      'sum': intensity_sum }
        metric_func = metric_map.get(metric, intensity_std)
    
        # (1) Create a grid of x and y coordinates to test --------------------
        x_range = np.linspace(px - search_range, px + search_range, steps)
        y_range = np.linspace(py - search_range, py + search_range, steps)
        score_grid = np.zeros((steps, steps))
    
        # (2) Calculate the metric score at each point on the grid ------------
        for i, y in enumerate(y_range):
            for j, x in enumerate(x_range):
                score_grid[i, j] = metric_func(self.parent.image, x, y, pr)
    
        # (3) Plot the results ------------------------------------------------
        plt.figure(figsize=(4, 4))
        plt.imshow(
            score_grid, 
            origin='lower', 
            extent=[x_range[0], x_range[-1], y_range[0], y_range[-1]]
            )
        plt.colorbar(label=f"Metric Score ('{metric}')")
        plt.contour(x_range, y_range, score_grid, colors='white', alpha=0.5)
        
        # Mark the starting point
        plt.plot(px, py, 'r+', markersize=15, label='Initial Guess')
        plt.title("Optimization Landscape")
        plt.xlabel("X-coordinate")
        plt.ylabel("Y-coordinate")
        plt.legend()
        plt.show()
        

class HandlerCircle(HandlerBase):
    """
    Helper class for creating circular markers in matplotlib legends.

    This class customizes the legend to display circular markers instead of the 
    default. It is intended for internal use within the module and not 
    for general use.

    Methods:
    --------
    create_artists(legend, 
                   orig_handle, 
                   xdescent,
                   ydescent, 
                   width, 
                   height, 
                   fontsize, 
                   trans):
        Creates a circular marker for the legend based on the original handle's 
        properties.

    Parameters for `create_artists`:
        legend : matplotlib.legend.Legend
            The legend instance where the custom marker will be used.
        orig_handle : matplotlib.artist.Artist
            The original handle containing the marker properties 
            (e.g., facecolor, edgecolor).
        xdescent : float
            Horizontal offset adjustment for the marker.
        ydescent : float
            Vertical offset adjustment for the marker.
        width : float
            Width of the legend entry.
        height : float
            Height of the legend entry.
        fontsize : float
            Font size of the legend text.
        trans : matplotlib.transforms.Transform
            Transformation applied to the marker's coordinates.

    Returns:
    --------
    list of matplotlib.patches.Circle
        A list containing a single circular marker artist.
    """

    def create_artists(self, legend, orig_handle, xdescent, ydescent, 
                       width, height, fontsize, trans):
        # Calculate the center and radius of the circle
        x = width / 2
        y = height / 2
        r = min(width, height) / 2

        # Create a circular marker using the properties of the original handle
        marker = Circle((x, y), r,
                        facecolor=orig_handle.get_facecolor(),
                        edgecolor=orig_handle.get_edgecolor(),
                        linewidth=orig_handle.get_linewidth(),
                        transform=trans)

        return [marker]


class PocketFFTBackend:
    """
    High-performance FFT backend leveraging SciPy's PocketFFT with optimized 
    defaults. This backend enhances the efficiency of Fast Fourier Transforms 
    (FFT), particularly for 2D transforms, which may experience up to a 50% 
    speed improvement over standard configurations.

    Attributes
    ----------
    __ua_domain__ : str
        Universal array domain identifier for NumPy's SciPy-based FFT interface.

    Methods
    -------
    __ua_function__(method, args, kwargs)
        Universal array function handler that dispatches FFT operations to the 
        corresponding PocketFFT implementation with optimized threading
        settings.
    """

    __ua_domain__ = "numpy.scipy.fft"

    @staticmethod
    def __ua_function__(method, args, kwargs):
        """
        Dispatches FFT operations to SciPy's PocketFFT with optimized worker 
        settings. This method dynamically selects the appropriate PocketFFT 
        function and applies threading optimizations by adjusting the number 
        of workers.

        Parameters
        ----------
        method : callable
            The FFT function being dispatched.
        args : tuple
            Positional arguments passed to the FFT function.
        kwargs : dict
            Keyword arguments passed to the FFT function, including
            the optional "workers" parameter for parallel execution.

        Returns
        -------
        result : object
            The output of the FFT function if supported; 
            otherwise, returns NotImplemented.
        """
        fn = getattr(_pocketfft, method.__name__, None)
        if fn is None:
            return NotImplemented
        
        workers = kwargs.pop("workers", CPU_COUNT)
        return fn(*args, workers=workers, **kwargs)

        
class IntensityCenter: 
    '''
    Simple center determination for a symmetric diffractogram.
    
    * The center is determined as a center of intensity.
    * This works well for simple, symmetric diffraction patters, which are:
      (i) without beamstopper, (ii) pre-centered, and (iii) powder-like.
    * A real-life example of a simple symmetric diffractogram:
      a good powder electron diffraction pattern from STEMDIFF software.
    * This class is a legacy from previous EDIFF versions;
      it is kept mostly for backward compatibility.
      The functions in this class can be (and should be)
      replaced by a simple call of ediff.center.CenterLocator object.
      
    >>> # Center determination in a simple symmetric diffraction pattern
    >>> # (center = just center_of_intensity, no refinement
    >>>
    >>> # (1) Old way = this (old, legacy) IntensityCenter class:
    >>> xc,yc = ediff.center.IntensityCenter.center_of_intensity(
    >>>     arr, csquare=30, cintensity=0.8)
    >>>
    >>> # (2) New way = newer (and more universal) CenterLocator class:
    >>> xc,yc = ediff.center.CenterLocator(
    >>>     arr, detection_method='intensity', csquare=30, cintensity=0.8)
    '''
    
    @staticmethod
    def center_of_intensity(arr, csquare=20, cintensity=0.8):
        '''
        Find center of intensity/mass of an array.
        
        Parameters
        ----------
        arr : 2D-numpy array
            The array, whose intensity center will be determined.
        csquare : int, optional, default is 20
            The size/edge of the square in the (geometrical) center.
            The intensity center is searched only within the central square.
            Reasons: To avoid other spots/diffractions and
            to minimize the effect of an intensity assymetry around center. 
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
        
        ## Return the final center
        return(xc,yc)
