'''
Module: ediff.center
--------------------
Find the center of a 2D diffraction pattern.

* The center determination may be surprisingly tricky in certain cases.
* Nevertheless, the user just calls the CenterLocator class as shown below.

>>> # Example: How to use CenterLocator and get results?
>>> import ediff as ed
>>>
>>> center = ed.center.CenterLocator(
>>>    input_image='some_diffractogram.png',
>>>    determination='intensity',
>>>    refinement='sum',
>>>    messages=False, final_replot=True)
>>>
>>> print('Determined center coordinates:', center.x1, center.y1)
>>> print('Refined center coordinates   :', center.x2, center.y2)
'''

# CenterDet
# PS 2023-10-06: CentDet update, methods compatibility
# MS 2023-11-26: Improved code formatting and docs + TODO notes for PS
# MS 2024-23-09: Re-desing/draft of CenterLocator => CenterLocator_new
#   # CenterLocator - beter structure, clearer usage, saving/reading to files
#   import ediff as ed
#   CENTER = ed.center.CenterlLocator(args)
#   print(CENTER.x1,CENTER.y1)  # center coords after CenterDetermination
#   print(CENTER.x2,CENTER.y2)  # center coords after CenterDetermination
    

import numpy as np
import skimage as sk
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.legend_handler import HandlerBase

import ediff.io 
import os

from skimage.measure import moments
from skimage.transform import hough_circle, hough_circle_peaks
from scipy.signal import find_peaks
from textwrap import dedent

import sys
import warnings
warnings.filterwarnings("ignore")


class CenterLocator:
    '''
    CenterLocator object for determining and refining the center 
    of a diffraction pattern.

    Parameters
    ----------
    input_image : str, path, or numpy.array
        The input image representing a 2D diffraction pattern, either as a 
        file path or a NumPy array.
    determination : callable, optional, default is None
        Method used for center determination
        = the initial estimate of center coordinates.
    refinement : callable, optional, default is None
        Method used for center refinement
        = refining the calculated center position interactively.
    in_file : str, optional, default is None
        Filename of the text file
        for saving the center coordinates.
    out_file : str, optional, default is None
        Filename of the text file
        for loading the previously saved center coordinates.
    heq : bool, optional, default is False
        Flag to indicate whether to perform histogram equalization 
        on the input image.
        The equalization is done internally,
        the image on the screen remains unchanged.
    icut : float, optional, default, default is None
        Cut-off intensity level for processing the image.
    cmap : str, optional, default is 'gray'
        Colormap to be used for displaying the image. 
    csquare : int, optional, default is 50
        Size of the central square,
        within which we will search the intensity center.
    cintensity : float, optional, default is 0.8
        Threshold intensity
        for finding the intensity center.
        Pixels with a (relative) intensity lower than cintensity are ignored.
    messages : bool, optional, default is False
        Flag to enable or disable informational messages during processing. 
    print_sums : bool, optional, default is False
        If True, prints the sum of intensity values for the refined circle 
        after each adjustment, relevant only for manual methods of center 
        determination and refinement.
    final_replot : bool, optional, default is False
        Flag to indicate whether to replot the final results.
 
    Returns
    -------
    None
        The center coordinates are stored in instance variables
        (x1,y1) for the determined center and
        (x2,y2) for the refined center.
        Look at the initial example at the top of this module
        to see how to use CenterLocator class.
            
    Technical notes
    ---------------
    * The class initializes (and runs) two sub-classes (= processes),
      CenterDetermination and CenterRefinement.
    * The two sub-classes/processes are hidden to common user,
      altough they could be run separately.
    '''
    
    
    def __init__(self,
                 input_image, 
                 determination = None, 
                 refinement = None,
                 in_file = None,
                 out_file = None,
                 heq = False, 
                 icut = None,
                 cmap = 'gray',
                 csquare=50,
                 cintensity=0.8,
                 messages = False,
                 print_sums = False,
                 final_print = True,
                 final_replot=False):
        
        ######################################################################
        # PRIVATE FUNCTION: Initialize CenterLocator object.
        # The parameters are described above in class definition.
        ######################################################################

        ## (0) Initialize input attributes
        self.input_image = input_image
        if refinement is not None:
            self.refinement = refinement.lower()
        else: self.refinement = refinement
        
        if determination is not None:
            self.determination = determination.lower()
        else: self.determination = determination

        self.in_file = in_file
        self.out_file = out_file
        self.heq = heq
        self.icut = icut
        self.cmap = cmap
        self.csquare = csquare
        self.cintensity = cintensity
        self.messages = messages
        self.print_sums = print_sums
        self.final_print = final_print
        self.final_replot = final_replot
        
        ## (1) Initialize new attributes
        self.to_refine = []
        self.dText = []
        self.rText = []
        # Adjust this initial value as desired
        self.marker_size = 100          

        # Allow input images (np.ndarray) and image path
        if isinstance(input_image, np.ndarray):
            self.image = input_image
        else:
            self.image = ediff.io.read_image(self.input_image)
        
        ## (2a) Initialize/run CenterDetermination
        self.center1 = CenterDetermination(self,
                self.input_image,
                self.determination,
                self.heq,
                self.icut,
                self.cmap,
                self.csquare,
                self.cintensity,
                self.messages,
                self.print_sums,
                self.final_replot)
        
        ## (2b) Correct radius
        if self.determination != "manual":
            self.center1.r = self.center1.get_radius(
                self.image, 
                self.center1.x, 
                self.center1.y, 
                disp=False)

        ## (3) Initialize/run CenterRefinement
        self.center2 = CenterRefinement(self,
            self.input_image, 
            self.refinement,
            self.in_file,
            self.out_file,
            self.heq, 
            self.icut,
            self.cmap,
            self.messages,
            self.print_sums,
            self.final_replot)
        
        ## (4a) Collect results
        self.x1 = self.center1.x
        self.y1 = self.center1.y
        
        if self.refinement is not None:
            self.x2 = self.center2.xx
            self.y2 = self.center2.yy
        else: 
            self.x2 = self.center1.x
            self.y2 = self.center1.y
            self.center2.xx = self.center1.x
            self.center2.yy = self.center1.y
            self.center2.rr = self.center1.r
        
        ## (4b) Switching coordinates if necessary
        if (self.determination == "intensity"):
                self.x1, self.y1 = self.convert_coords(self.x1, self.y1)
                self.x2, self.y2 = self.convert_coords(self.x2, self.y2)  
        
        if (self.determination =="hough" and self.refinement == "manual"):
            self.x2, self.y2 = self.convert_coords(self.x2, self.y2)  

        ## (5a) Printing coordinates
        if self.final_print:
            self.dText=str(self.dText)
            self.rText=str(self.rText)
            
            print("----------------------------------------------------------")
            print(self.dText.format(float(self.x1),float(self.y1)))
            
            if self.refinement is not None:
                print(self.rText.format(float(self.x2),float(self.y2)))
                
        ## (5b) Plot results if final_replot
        if final_replot:   
            self.visualize_refinement(
                self.x1, self.y1, self.center1.r,
                (self.x2, self.y2), self.center2.rr)
        
        ## (6) Save results to a .txt file if specified
        if out_file is not None:
            self.save_results()
            
        ## (7) Load results from a .txt file if specified       
        if in_file is not None:
            self.load_results()   
        
    
    def output(self):
        """
        Manage variables that should be send as the output of the center 
        detection. 
        
        If there were set parameters detection_method and 
        correction method during the class initialization, the output will be
        coordinates x, y of the center detected by the detection_method and 
        coordinates x, y of refined center position by the correction method.
        
        If there was not set the correction_method parameter, the function
        outputs x, y coordinates of the detected center and None, None for
        the refined coordinates.

        Returns
        -------
        x : float
            x-coordinate of the center detected via detection_method
        y : float
            y-coordinate of the center detected via detection_method
        xx : float
            x-coordinate of the center detected via refinement_method
        yy : float
            y-coordinate of the center detected via refinement_method
        """
        
        if self.center2.ret == 1:
            # Convert to float
            if type(self.x1) != float:
                self.x1, self.y1, self.x2, self.y2 = \
                    [float(value) for value in (self.x1, self.y1, 
                                                self.x2, self.y2)]
                        
            # Return values of center coordinates
            return (np.round(self.x1,1), np.round(self.y1,1), 
                    np.round(self.x2,1), np.round(self.y2,1))  
        else:
            # Convert to float
            if type(self.x1) != float:
                self.x1, self.y1 = \
                    [float(value) for value in (self.x1,self.y1)]

            # Return values of center coordinates
            return (np.round(self.x1,1), np.round(self.y1,1), None, None)   
    
    
    def save_results(self):
        '''
        Save the current results to a specified file. The results are formatted 
        to four decimal places for clarity.

        
        This method checks if a file path has been provided through 
        the `in_file` attribute. 
        
        - If the specified file exists, the method appends the results 
        (x1, y1, x2, y2) to the end of the file. 
        - If the file does not exist, it creates a new file and writes 
        the results to it.
        
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        Notes
        -----
        - The results are written in the format:
            x1: <value>, y1: <value>
            x2: <value>, y2: <value>
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


    def load_results(self):
        '''
        Load results from a specified text file.
    
        This method reads the coordinates from a text file defined by the 
        `out_file` attribute. It extracts pairs of coordinates (x1, y1) and 
        (x2, y2) from each line in the file. The extracted values are stored in 
        lists, allowing for the retrieval of multiple results.The most recent 
        values of x1, y1, x2, and y2 can be stored as instance variables
        for easy access.
    
        Parameters
        ----------
        None
    
        Returns
        -------
        None
    
        Notes
        -----
        - The results are expected to be in the format:
            x1: <value>, y1: <value>
            x2: <value>, y2: <value>
        - If multiple sets of coordinates are found, they are stored in lists,
        and only the last set can be accessed as instance variables.
        '''        
        
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
    
        else:
            print("Input file does not exist or is not specified.")


    def get_circle_pixels(self, xc, yc, radius, num_points=360):
        '''         
        Get coordinates of pixels defining circle border
    
        Parameters
        ----------
        self.image_path : str
            direct path to a image with diffraction patterns
        xc : float64
            x-coordinate of the detected center
        yc : float64
            y-coordinate of the detected center
        radius : float64
            radius of the detected center
        num_points : float64 
            number of border points. The default is 360
        
        Returns
        -------
        x : array of float64
            x-coordinates of pixels from circle border
        y : array of float64
            y-coordinates of pixels from circle border
            
        '''
        
        # Generate angles from 0 to 2*pi
        theta = np.linspace(0, 2*np.pi, num=num_points)  
                
        # Calculate x,y-coordinates of points on the actual circle border
        x_actual = xc + radius * np.cos(theta)
        y_actual = yc + radius * np.sin(theta)
        
        return x_actual, y_actual
          
    
    def intensity_sum(self, image, px, py, pr):
        ''' 
        Summation of intensity values of pixels of a diffraction pattern.

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
            intensity sum

        '''
        # Extract pixels on the circle border
        pxc, pyc = self.get_circle_pixels(px, py, pr)
        pxc = np.array(pxc, dtype=int)
        pyc = np.array(pyc, dtype=int)
        
        # Calculate sum using the filtered values
        s = np.sum(image[pyc, pxc])/len(pxc)

        return s
    
    
    def intensity_var(self, image, px, py, pr):
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
    
    
    def visualize_refinement(self, px, py, pr, xy, r):
        '''
        Visualize diffraction patterns and center after correction
    
        Parameters
        ----------
        px : float64
            x-coordinate before correction.
        py : float64
            y-coordinate before correction.
        pr : float64
            radius before correction.
        xy : float64
            xy-coordinates after correction.
        r : float64
            radius after correction.
    
        Returns
        -------
        None.
    
        '''
    
        image = np.copy(self.to_refine)
    
        if self.icut is not None:
            im = np.where(image > self.icut, self.icut, image)
        else:
            im = np.copy(image)
    
       # ediff.io.set_plot_parameters(size=(9,9), dpi=100, fontsize=10)
        fig, ax = plt.subplots()
    
        if self.refinement == "var":
            dvar = self.intensity_var(image, px, py, pr)
            rvar = self.intensity_var(image, xy[0], xy[1], r)
            labeld = f'd-center: [{px:.1f}, {py:.1f}]\nint-var: {dvar:.1f}'
            labelr = f'r-center: [{xy[0]:.1f}, {xy[1]:.1f}]\nint-var: {rvar:.1f}'
        else:
            dsum = self.intensity_sum(image, px, py, pr)
            rsum = self.intensity_sum(image, xy[0], xy[1], r)
            labeld = f'd-center: [{px:.1f}, {py:.1f}]\nint-sum: {dsum:.1f}'
            labelr = f'r-center: [{xy[0]:.1f}, {xy[1]:.1f}]\nint-sum: {rsum:.1f}'
    
        # Original Image
        ax.imshow(im, cmap=self.cmap, origin="upper")
        c0 = plt.Circle((px, py), pr,
                        color='r',
                        fill=False,
                        label='detected',
                        linewidth=1)
        ax.add_patch(c0)
        ax.scatter(px, py,
                   label=labeld,
                   color='r',
                   marker='x',
                   s=60, linewidths=1)
    
        # Refined Image
        c1 = plt.Circle(xy, r,
                        color='springgreen',
                        fill=False,
                        label='refined',
                        linewidth=1)
        ax.add_patch(c1)
    
        ax.scatter(xy[0], xy[1],
                   label=labelr,
                   color='springgreen',
                   marker='x',
                   linewidths=1, s=60)
    
        ax.set_title(
            f'Center Location ({self.determination}/{self.refinement})')
          #  fontsize=12)
        
        # Move legend to the top right next to the image
        ax.legend(loc='upper left', frameon=False, #fontsize=10,
                  handler_map={Circle: HandlerCircle()}, ncol=1,
                  bbox_to_anchor=(.98, 1.01), 
                  bbox_transform=ax.transAxes)

    
        ax.axis('on')
    
        #plt.tight_layout()
        plt.show(block=False)

    
    def convert_coords(self, x, y):
        """
        Convert coordinates between numpy and matplotlib systems.
    
        Parameters:
        ----------
        x : int or float
            The x-coordinate in the numpy (column index) format.
        y : int or float
            The y-coordinate in the numpy (row index) format.
    
        Returns:
        -------
        tuple of (int or float, int or float)
            The converted coordinates in matplotlib format, where:
            - First element corresponds to y (new x in matplotlib).
            - Second element corresponds to x (new y in matplotlib).
        """
        return y, x




class CenterDetermination:
    '''
    CenterDetermination object for detecting the center of a diffraction 
    pattern.

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
    determination : str, optional
        The method used for center detection. Options include:
        - 'manual': Manual detection using three points.
        - 'hough': Hough transform method for circle detection.
        - 'intensity': Detection based on intensity thresholds.
    in_file : str, optional
        Filename of the text file from which to load previously saved 
        center coordinates.
    out_file : str, optional
        Filename of the text file to which the detected center coordinates 
        will be saved.
    heq : bool, optional
        Flag to indicate whether to perform histogram equalization on the 
        input image. Default is False.
    icut : float, optional
        Cut-off intensity level for processing the image. Default is None.
    cmap : str, optional
        Colormap to be used for displaying the image. Default is 'gray'.
    csquare : int, optional
        Size of the square for processing. Default is 50.
    cintensity : float, optional
        Threshold intensity for detecting features in the image. Default is
        0.8.
    messages : bool, optional
        Flag to enable or disable informational messages during processing. 
        Default is False.
    print_sums : bool, optional
        If True, prints the sum of intensity values for the detected circle
        after each adjustment, relevant only for manual detection methods.
    final_replot : bool, optional
        Flag to indicate whether to replot the final results. 
        Default is False.

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
                 messages = False,
                 print_sums = False,
                 final_replot=False):
        
        ######################################################################
        # PRIVATE FUNCTION: Initialize CenterLocator object.
        # The parameters are described above in class definition.
        ######################################################################
        
        ## (0) Initialize input attributes
        self.parent = parent
        
        ## (1) Initialize new attributes
        self.step=0.5
        
        ## (2) Run functions
        # (2a) Preprocess data
        self.preprocess(preInit=1)
        
        # (2b) Center detection methods
        if determination.lower() == "manual":
            self.x, self.y, self.r = self.detection_3points()
        elif determination.lower() == "hough":
            self.x, self.y, self.r = self.detection_Hough()
        elif determination.lower() == "intensity":
            self.x, self.y, self.r = self.detection_intensity(
                self.parent.csquare, 
                self.parent.cintensity
                )
                
        else: 
            print("Selected determination method does not exist.")
            sys.exit()


    ## Functions::
    def preprocess(self,preInit=0,preHough=0,preManual=0,
                   preVar=0,preSum=0,preInt=0):  
        """ FOR AUTOMATIC METHODS OPTIMIZATION AND MORE UNIVERSAL SOLUTIONS
        Function for input image preprocessing based on the methods 
        defined in the class initialization - self.detection_method, 
        self.correction_method.

        Parameters
        ----------
        preInit : bool, optional
            Perform preprocessing of the input image (when using icut or heq). 
            This is called automatically every time, if no preprocessing
            specified, the detection and refinement will be performed on 
            original image. The default is 0.
        preHough : bool, optional
            Perform preprocessing for automatic detection via Hough transform. 
            

        Returns
        -------
        manu : NumPy array
            Pre-processed image for the manual detection method
        edges : array of bool
            Detected edges via Canny detector for automatic Hough transform
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
                if self.parent.messages:
                    print("Histogram equalized.")
                image = sk.exposure.equalize_adapthist(image)

                
            # # Edit contrast with a user-predefined parameter
            # if self.parent.icut is not None:
            #     if self.parent.messages:
            #         print("Contrast enhanced.")
            #     image = np.where(image > self.parent.icut, 
            #                      self.parent.icut, 
            #                      image)
                

            self.parent.to_refine = image
            return
        
        ### Hough transform: perform pre-processing necessary for the detection
        if preHough == 1:           
            if self.parent.heq == 0:
                csq = self.central_square(self.parent.to_refine, csquare=80)   
                
                # Beam stopper present in image
                if np.median(csq)<100 and np.median(csq) > 0:
                    if self.parent.messages:
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
                    
                    # Erode the dilated edges to reduce thickness and smooth the contour
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
                    
                    # Erode the dilated edges to reduce thickness and smooth the contour
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
        
    
    def detection_intensity(self, csquare, cintensity, plot_results=0):
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
        
        # (1) Get image/array size
        image = np.copy(self.parent.to_refine)
        
                
        arr = np.copy(image)
        xsize,ysize = arr.shape
        
        # (2) Calculate borders around the central square
        xborder = (xsize - csquare) // 2
        yborder = (ysize - csquare) // 2
        
        # (3) Create the central square,
        # from which the intensity center will be detected
        arr2 = arr[xborder:-xborder,yborder:-yborder].copy()
        
        # (4) In the central square,
        # set all values below cintenstity to zero
        arr2 = np.where(arr2>np.max(arr2)*cintensity, arr2, 0)
        
        # (5) Determine the intensity center from image moments
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
        # (d) Radius of a diffraction ring is hardcoded to 100 here
        self.r = 100
        
        # (6) User information (if required)
        if (self.parent.messages or self.parent.final_print):
            self.parent.dText = "Center Determination (IntensityCenter): ({:.3f}, {:.3f})"


        # (7) Plot results (if required)
        if plot_results == 1:
           self.visualize_center(self.x, self.y, self.r) 
                
        # (8) Return the results:
        # a) XY-coordinates of the center
        # b) radius of the circle/diff.ring,
        #    from which the center was determined
        # ! radius of the circle is just hardcoded here - see TODO note above
        return(self.x, self.y, self.r)
    

    def detection_Hough(self, plot_results=0):
        '''        
        Perform Hough transform to detect center of diffraction patterns.
        This is a method to automatically detect circular diffraction patterns
        
        Parameters
        ----------
        plot_results : int, binary
            Plot the pattern determined by pixels selected by the user.
            Default is 1. To cancel visualization, set plot_results = 0.

        Returns
        -------
        self.x : float64
            x-coordinate of the detected center
        self.y : float64
            y-coordinate of the detected center
        self.r : float64
            radius of the detected center
                                    
        '''
        ## Image preprocessing
        im = np.copy(self.parent.to_refine)
        
        # if the brightness of the image is small enough, pixel values greater
        # than 50 will be set to 0 -- removal of the beam stopper influence
        
        if self.parent.heq == 0:
            if sum(sum(im)) < 150000:
                    max_indices = np.where(im > 50)
        
                    row_idx = max_indices[0]
                    col_idx = max_indices[1]
        
                    im[row_idx, col_idx] = 0    
                
            # Detect edges using the Canny edge detector
            edges = sk.feature.canny(im, 
                                     sigma=0.2, 
                                     low_threshold=80, 
                                     high_threshold=100)
        elif self.heq == 1:
            if sum(sum(im)) > 40000:
                max_indices = np.where(im > 50)

                row_idx = max_indices[0]
                col_idx = max_indices[1]

                im[row_idx, col_idx] = 0    
            
            # Detect edges using the Canny edge detector
            edges = sk.feature.canny(im, 
                                     sigma=0.2, 
                                     low_threshold=0.80, 
                                     high_threshold=1)
        
        
        # Define the radii range for the concentric circles
        # (set empirically based on the available pictures)
        min_radius = 40
        max_radius = 200
        radius_step = 10
        radii = np.arange(min_radius, max_radius + radius_step, radius_step)

        ### Perform the Hough transform to detect circles
        # Circle detection involves converting edge pixels into parameter space, 
        # where each point represents a possible circle center and radius. 
        # The circles are then identified as peaks in the parameter space, 
        # enabling accurate detection of circular shapes in the image.
        hough_res = hough_circle(edges, radii)

        # Extract the circle peaks
        _, self.x, self.y, self.r = hough_circle_peaks(hough_res, 
                                                            radii, 
                                                            total_num_peaks=1)
        
        
        # User information:
        if (self.parent.messages or self.parent.final_print):
            self.parent.dText = "Center Determination (HoughTransform) : ({:.3f}, {:.3f})"
        
        self.x, self.y, self.r = \
            float(self.x[0]), float(self.y[0]), float(self.r[0])

        # (7) Plot results (if required)
        if plot_results == 1:
           self.visualize_center(self.x, self.y, self.r) 
           
        # Return results, convert coordinates to float
        return self.x, self.y, self.r

   
    def detection_3points(self, plot_results=0):
        '''         
        In the input image, select manually 3 points defining a circle using
        a key press event 
            - press '1' to select a point
                
        If the user is not satisfied with the point selection, it can be
        deleted using a key press event:
            - press '2' to delete the most recent
            - press '3' to delete the point closest to the cursor
        
        If the user is satisified with the points selected, the rest 
        of the program will be executed 
            - press 'd' to proceed >DONE<
        
        Coordinates of the center and radius will be calculated automatically
        using method self.calculate_circle()
        
        In addition, the user will be able to manually adjust the original
        center position by using pre-defined keys.
        
        Parameters
        ----------
        plot_results : int, binary
            Plot the pattern determined by pixels selected by the user.
            Default is 1. To cancel visualization, set plot_results = 0.
        
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
        
        # Load image
        im = self.parent.to_refine
        
        # Edit contrast with a user-predefined parameter
        if self.parent.icut is not None:
            if self.parent.messages:
                print("Contrast enhanced.")
            im = np.where(im > self.parent.icut, 
                              self.parent.icut, 
                              im)
            
        # Create a figure and display the image
        fig, ax = plt.subplots(figsize=(12, 12))
        
        # Allow using arrows to move back and forth between view ports
        plt.rcParams['keymap.back'].append('left')
        plt.rcParams['keymap.forward'].append('right')
 
        plt.title("Select 3 points defining one of diffraction circles", 
                  fontsize = 20)
        ax.imshow(im, cmap = self.parent.cmap, origin="upper")
        ax.axis('off')

        # User information:
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

        if self.parent.messages:
            print(instructions)
       
        # Enable interactive mode
        # (figure is updated after every plotting command
        # (so that calling figure.show() is not necessary
        plt.ion()
 
        # Initialize the list of coordinates
        self.coords = [] 
        
        # Initialize all flags and counters
        calculate_circle_flag = False          # press 'd' event
        termination_flag = False               # close window event
        point_counter = 0                      # number of selected points
        
        ### Define the event handler for figure close event
        def onclose(event):
            nonlocal termination_flag
            termination_flag = True
            if self.parent.messages:
                print('Execution terminated.')
                print("------------------------------------------------------------")

 
        # Connect the event handler to the figure close event
        fig.canvas.mpl_connect('close_event', onclose)
        
        
        ### Define the callback function for key press events
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
    

    
            # Calculate circle or terminate
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
        '''
        Adjustment of the center position calculated from 3 points.
        Interactive refinement using keys:

        The user can change the position of the center of the diffraction
        pattern and also the radius of the detected pattern using keys:
            - left / right / top / down arrows : move left / right / top / down
            - '+' : increase radius
            - '-' : decrease radius
            - 'd' : done, termination of the refinement

        If the interactive figure is closed without any modifications,
        the function returns input variables and the proccess terminates.
        
        Parameters
        ----------
        fig : figure.Figure object
            interactive figure in which a diffraction pattern has been
            manually detected.
        circle : patches.Circle object
            circle defined via 3 points manually delected
        center : tuple
            calculated center of the input circle.
        plot_results : boolean
            visualize results. The default is 1 (plot detected center).

        Returns
        -------
        xy : tuple
            x,y-coordinates of the center of the diffraction pattern.
        r : integer
            radius of the diffraction pattern.

        '''            
        # Remove default left / right arrow key press events
        plt.rcParams['keymap.back'].remove('left')
        plt.rcParams['keymap.forward'].remove('right')
        
        if self.parent.messages:
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
        if (self.parent.messages or self.parent.final_print):
            self.parent.dText = "Center Determination (ThreePoints)    : ({:.3f}, {:.3f})"
        
    
        return xy[0], xy[1], r
    
       
    def calculate_circle(self, plot_results:int)->tuple[
            float,float,float,tuple[float,float], plt.Circle]:
        ''' 
        Calculates coordinates of the center and radius of a circle defined via
        3 points determined by the user. Plots the calculated circle, detected 
        points and marks the center.
        
        Parameters
        ----------
        plot_results : int, binary
            Plot the calculated center and circle. To cancel visualization, 
            set plot_results = 0.
        self.coords : array of float64
            Coordinates of 3 manually selected points
        
        Returns
        -------
        self.x : float64
            x-coordinate of the detected center
        self.y : float64
            y-coordinate of the detected center
        self.r : float64
            radius of the detected center
                                    
        '''
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
        # if self.parent.messages:
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
            self, im:np.ndarray, x:float, y:float, disp:bool=False) -> float:
        """
        Calculate the radius of a circle based on intensity profiles along 
        horizontal and vertical axes.
    
        Parameters
        ----------
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
        
        def match_peaks(arr):
            """
            Finds the most similar values across the left and right halves of 
            an array, omitting the highest value. If exactly two peaks are 
            detected, they are returned as the best pair.
            """
            if len(arr) < 2:  # Not enough values to compare
                if self.parent.messages:
                    print("Not enough values to find similar pairs.")
                return None, None
    
            # Special case: if there are exactly 2 peaks, return them
            if len(arr) == 2:
                if self.parent.messages:
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
                if self.parent.messages:
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
            if self.parent.messages:
                print("xpeaks condition met: Both peaks are on the same side of the center.")
            self.pairX = None
        else:
            self.pairX = match_peaks(self.xyvals)
    
        # Check the additional condition for ypeaks
        if len(self.ypeaks) == 2 and (
            (self.ypeaks[0]<half_length_y and self.ypeaks[1]<half_length_y) or
            (self.ypeaks[0]>half_length_y and self.ypeaks[1]>half_length_y)):
            if self.parent.messages:
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
            if self.parent.messages:
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
    
        return rx

        
    def visualize_center(self, x: float, y: float, r: float) -> None:
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
    
        if self.parent.icut is not None:
            im = np.where(image > self.parent.icut, self.parent.icut, image)
        else:
            im = np.copy(image)
            
        # Create a figure and display the image
        fig, ax = plt.subplots()
        
        # Allow using arrows to move back and forth between view ports
        plt.rcParams['keymap.back'].append('left')
        plt.rcParams['keymap.forward'].append('right')
 
        plt.title("Detected center")
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
    CenterRefinement object for refining the center coordinates of 
    a diffraction pattern.

    This class is responsible for refining the center coordinates of a 
    diffraction pattern based on the selected refinement method. It 
    initializes with the input image and other parameters that influence 
    the refinement process.

    Parameters
    ----------
    parent : CenterLocator
        Reference to the parent CenterLocator object, allowing access to 
        shared attributes and methods.
    input_image : str, path, or numpy.array
        The input image representing a 2D diffraction pattern, provided as 
        a file path or a NumPy array.
    refinement : str, optional
        The method used for center refinement. Options include:
        - 'manual': Manual adjustment of the center based on user input.
        - 'var': Variance-based refinement to optimize the center coordinates.
        - 'sum': Sum-based refinement to find the optimal center based on 
                 intensity sums.
    in_file : str, optional
        Filename of the text file from which to load previously saved 
        refined coordinates.
    out_file : str, optional
        Filename of the text file to which the refined center coordinates 
        will be saved.
    heq : bool, optional
        Flag to indicate whether to perform histogram equalization on the
        input image. Default is False.
    icut : float, optional
        Cut-off intensity level for processing the image. Default is None.
    cmap : str, optional
        Colormap to be used for displaying the image. Default is 'gray'.
    messages : bool, optional
        Flag to enable or disable informational messages during processing. 
        Default is False.
    print_sums : bool, optional
        If True, prints the sum of intensity values for the refined circle 
        after each adjustment, relevant only for manual refinement methods.
    final_replot : bool, optional
        Flag to indicate whether to replot the final results after
        refinement. Default is False.

    Returns
    -------
    None

    Notes
    -----
    - The class refines the detected center coordinates based on the 
      specified refinement method.
    - The refined coordinates are stored in instance variables
      `xx`, `yy`, and `rr`, representing the refined center's x-coordinate, 
      y-coordinate, and radius, respectively.
    - If an unsupported refinement method is specified, the class will print 
      an error message and exit.
    '''
    
    def __init__(self, parent,
                 input_image, 
                 refinement = None,
                 in_file = None,
                 out_file = None,
                 heq = False, 
                 icut = None,
                 cmap = 'gray',
                 messages = False,
                 print_sums = False,
                 final_replot=False):
        
        ######################################################################
        # PRIVATE FUNCTION: Initialize CenterLocator object.
        # The parameters are described above in class definition.
        ######################################################################
        
        ## (0) Initialize input attributes
        self.parent = parent
        
        ## (1) Initialize new attributes
        self.step=0.5
        
        ## (2) Run functions
        if refinement is not None:
            self.ret = 1
            par_short = self.parent.center1
        
            if refinement == "manual":
                # Manual refinement method
                if parent.determination == 'manual':
                    self.xx, self.yy, self.rr = \
                        par_short.x, par_short.y, par_short.r
                    par_short.x, par_short.y = \
                        par_short.backip[0], par_short.backip[1]
                    if (self.parent.messages or self.parent.final_print):
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
                    par_short.x, par_short.y, par_short.r)
        
            else:
                print("Selected refinement method is not supported.")
                sys.exit()
                        
        else: 
            self.ret = 2
            
           
    def ref_interactive(self, px, py, pr):
        ''' 
        Manual refinement of the detected diffraction pattern via one of 
        the methods provided in the class CircleDetection.
        
        The user can change the position of the center of the diffraction
        pattern and also the radius of the detected pattern using keys:
            - left / right / top / down arrows : move left / right / top / down
            - '+' : increase radius
            - '-' : decrease radius
            - 'b'/'l' : increase/decrease step size")
            - 'd' : done, termination of the refinement

        If the interactive figure is closed without any modifications,
        the function returns input variables and the proccess terminates.
        
        The results are shown in a figure when the refinement is successful.

        Parameters
        ----------
        px : float64
            x-coordinate of the center
        py : float64
            y-coordinate of the center
        pr : float64
            radius of the circular diffraction pattern

        Returns
        -------
        x : float64
            new x-coordinate of the center
        y : float64
            new y-coordinate of the center
        r : float64
            new radius of the circular diffraction pattern
        '''
        
        # Load original image
        im = np.copy(self.parent.to_refine)

        # Edit contrast with a user-predefined parameter
        if self.parent.icut is not None:
            if self.parent.messages:
                print("Contrast enhanced.")
            im = np.where(im > self.parent.icut, 
                              self.parent.icut, 
                              im)
            
        
        # Initialize variables and flags
        xy = np.array((px, py))
        r = np.copy(pr)
        termination_flag = False

        # User information:
        if self.parent.messages:
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
            
        # Create a figure and display the image
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
        
        # Enable interactive mode
        plt.ion()
        

        # Display the image
        # fig.set_size_inches(self.fig_width, self.fig_height)
        plt.show(block=False)
        
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
        if (self.parent.messages or self.parent.final_print):
            self.parent.rText = "Center Refinement (Interactive)       : ({:.3f}, {:.3f})"

        return xy[0], xy[1], r    
        
    
    def ref_var(self, px, py, pr, plot_results=0):
        '''         
        Adjust center coordinates of a detected circular diffraction pattern.
        The center adjustment is based on variance minimization.
        
        The 8-neighbourhood pixels (x) of the current center (o) 
        will be tested regarding the minimization:
    
        - x x x : (px - dx, py + dy) (px, py + dy) ( px + dx, py + dy)
    
        - x o x : (px - dx, py)      (px, py)      (px + dx, py)
    
        - x x x : (px - dx, py - dy) (px, py - dy) (px + dx, py - dy)
        

        Parameters
        ----------
        self.image : array of uint8
            Input image in which the diffraction pattern is to be found
        px : float64
            x-coordinate of the detected center to be adjusted
        py : float64
            y-coordinate of the detected center to be adjusted
        pr : float64
            radius of the detected center
        plot_results : integer (default = 1)
            Plot Detected center. The default is 1.
        
        Returns
        -------
        px : array of int32
           corrected x-coordinates of pixels from circle border
        py : array of int32
            corrected y-coordinates of pixels from circle border
        pr : array of int32
            radius of the detected center
            
        '''
        
        # Store input for plot
        bckup = [np.copy(px), np.copy(py), np.copy(pr)]
        
        # Load image
        image = np.copy(self.parent.image)

        # Starting values to be modified 
        init_var = self.parent.intensity_var(image, px, py, pr)
        min_intensity_var = self.parent.intensity_var(image, px, py, pr)
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
            curr = self.parent.intensity_var(image, 
                                      best_center[0], 
                                      best_center[1], 
                                      best_radius) 
            # Store intensity sums of the current center's neighborhood
            curr_intensity_var = []
            for dx, dy in neighbors:
                nx, ny = best_center[0] + dx, best_center[1] + dy
                # Check if the point is within the expanded search radius
                curr_intensity_var.append(self.parent.intensity_var(image, 
                                                             nx, ny, 
                                                             best_radius))
            
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
            min_intensity_var = self.parent.intensity_var(image, 
                                                    best_center[0], 
                                                    best_center[1], 
                                                    best_radius) 
            
            # Refine radius if necessary while keeping the center position 
            # constant. It iterates through different radius adjustments to find
            # a radius that maximizes the intensity sum of pixels
            
            radi_intensity_var = []
            radii = np.arange(-1, 1 + step, step)
            for dr in radii:
                new_radius = best_radius + dr
                radi_intensity_var.append(self.parent.intensity_var(image, 
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
        if (self.parent.messages or self.parent.final_print):
            self.parent.rText = "Center Refinement (IntensityVar)      : ({:.3f}, {:.3f})"
                                  
        return best_center[0], best_center[1], best_radius
    
    
    def ref_sum(self, px, py, pr, plot_results=0):
        ''' 
        Adjust center position based on gradient optimization method
        via maximization of intensity sum.
        
        The 8-neighbourhood pixels (x) of the current center (o) 
        will be tested regarding the maximization:
    
        - x x x : (px - dx, py + dy) (px, py + dy) ( px + dx, py + dy)
    
        - x o x : (px - dx, py)      (px, py)      (px + dx, py)
    
        - x x x : (px - dx, py - dy) (px, py - dy) (px + dx, py - dy)
        
    
        Parameters
        ----------
        px : float64
            x-coordinate of the detected center to be adjusted.
        py : float64
            y-coordinate of the detected center to be adjusted.
        pr : float64
            radius of the detected center.
        plot_results : int, optional
            Plot Detected center. 
            The default is 1.
    
        Returns  
        -------
        best_center[0] : float64
            Adjusted x-coordinate of the center.
        best_center[1] : float64
            Adjusted y-coordinate of the center.
        best_radius : float64
            The adjusted radius of the circular diffraction pattern.
        '''
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
                    
            # Check for improvement of criterion -- in each iteration just once,
            # as the algorithm checks the neighbourhood of the best center (in
            # each iteration, the center is updated if possible)
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
            # constant. It iterates through different radius adjustments to find
            # a radius that maximizes the intensity sum of pixels
            
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
                

        
        # Avoid incorrect/redundant refinement
        # ## (1) swapped coordinates
        # if ((bckup[0] > bckup[1] and not best_center[0] > best_center[1])
        #     or  (bckup[0] < bckup[1] and not best_center[0] < best_center[1])):
        #     best_center = best_center[::-1]
        
        ## (2) worsened final maximum intensity sum than the initial one
        if np.round(init_sum,-2) > np.round(max_intensity_sum,-2):
            print("Refinement redundant.")
            best_center = np.copy(bckup)
    
        # Print results
        if (self.parent.messages or self.parent.final_print):
            self.parent.rText = \
                "Center Refinement (IntensitySum)      : ({:.3f}, {:.3f})"
                
        return best_center[0], best_center[1], best_radius
    
    
    
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




    

    

