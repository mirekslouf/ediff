'''
Module ediff.center
-------------------
Find center of 2D diffraction pattern.    
'''

import numpy as np
import skimage as sk
import matplotlib.pyplot as plt

from skimage.io import imread
from skimage.exposure import equalize_adapthist
from skimage.feature import canny
from skimage.measure import moments
from skimage.transform import hough_circle, hough_circle_peaks

import random
import warnings
warnings.filterwarnings("ignore")

# The following extra module and its functions commented out
# Reasons: 
#  1) not necessary, default window sizes are fine and/or adjustable
#  2) we try to use only standard scientific python modules, if possible

# from screeninfo import get_monitors


class CenterDetection:
    '''
    Detection of the center of diffraction patterns.
    
    - Perform Hough transform to detect circles in image
    - Manually selected 3 points defining a circular diffraction pattern 
      to calculate center coordinates
    - Visualize results
    
    SUBCLASS : CenterAdjustment
    - adjust position of the detected center 
        
    Parameters
    ----------
    image_path : string
        direct path to an image with diffraction patterns
    detection_method : string
        Selection of a method for center calculation. String codes are:
            - 'manual' : manual detection via 3 points
            - 'intensity' : detection via maximum intensity
            - 'hough' : automatic detection via Hough transform
    correction_method : string
       Selection of a method for center position correction. Default is None.
       String codes are:
           - 'manual' : manual corection 
           - 'variance' : correction via variance minimization 
           - 'sum' : correction via sum maximization
    heq : boolean
            Allow histogram equalization. The default is 0 (no enhancement)
    icut : boolean
        Allow image enhancement. The default is 0 (no enhancement)

    
    Returns
    -------    
    self.x : float64
        x-coordinate of the detected center
    self.y : float64
        y-coordinate of the detected center
    self.r : float64
        radius of the detected center (if available, othervise returns None)
                    
    '''
    
    def __init__(self, image_path, detection_method, correction_method = None, heq = 0, icut = 0):
        # Initialize attributes
        self.image_path = image_path
        self.image = imread(self.image_path, as_gray = True)
        self.correction_method = correction_method
        self.heq = heq
        self.icut = icut
        
    
        # Get the list of monitors
        # monitors = get_monitors()
        # screen_width = monitors[0].width
        # screen_height = monitors[0].height
        # self.fig_width = screen_width / 2 / 100  # Convert to inches
        # self.fig_height = screen_height / 100  # Convert to inches
        

        # Determine detection method
        if detection_method == 'manual':
            self.x, self.y, self.r = self.detection_3points()
            
        elif detection_method == 'intensity':
            self.x, self.y, self.r = self.detection_intensity()

        elif detection_method == 'hough':
            self.x, self.y, self.r = self.detection_Hough()
            
    
    def detection_3points(self, plot_results=1):
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
        
        Parameters
        ----------
        self.image : array of uint8
            Input image in which the diffraction pattern is to be found
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
            radius of the detected center (if available, othervise returns None)
                                    
        '''
        
        # Load image
        im = np.copy(self.image)
        
        # Create a figure and display the image
        fig, ax = plt.subplots()

        plt.title("Select 3 points defining one of diffraction circles")
        ax.imshow(im)

        
        # User information:
        print("------------------- Manual diffraction pattern detection -----------------")
        print("--------------------------------------------------------------------------")
        print("Select 3 points to define a circular diffraction pattern.")
        print("Use these keys for the selection:")
        print("      - '1' : select a point")
        print("      - '2' : delete the most recent")
        print("      - '3' : delete the point closest to the cursor")
        print("      - 'd' : proceed (selection DONE)")
        print("Close the figure to terminate. No center will be detected.")
        print("--------------------------------------------------------------------------")
        
       
        # Enable interactive mode
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
            print('Execution terminated.')
        
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
                    distances = [np.sqrt((x - pointer_x)**2 + (y - pointer_y)**2) for x, y in self.coords]
                    closest_index = np.argmin(distances)
                    del self.coords[closest_index]
    
                    # Redraw the image without the deleted point
                    ax.clear()
                    ax.imshow(self.image) 
                    for x, y in self.coords:
                        ax.plot(x, y, 'rx')
                    plt.title("Select 3 points defining one of diffraction circles")
                    
                    # Retore the previous zoom level
                    ax.set_xlim(current_xlim)
                    ax.set_ylim(current_ylim)

                    fig.canvas.draw()
                else:
                    print("No points to delete. You must select at lease one point.")
    
            ## Delete recent point (last added) -- independent on the cursor
            if event.key == '2':
                # Check if there are points to delete
                if point_counter > 0:  
                    point_counter -= 1
                    if len(self.coords) > 0:
                        # Delete the last point in the list
                        del self.coords[-1]  
                        print('The most recently selected point deleted.')
                        print('Please select a new one.')
    
                        # Redraw the image without the deleted point
                        ax.clear()
                        ax.imshow(self.image)
                        for x, y in self.coords:
                            ax.plot(x, y, 'rx')
                        plt.title("Select 3 points defining one of diffraction circles")
                        
                        # Retore the previous zoom level
                        ax.set_xlim(current_xlim)
                        ax.set_ylim(current_ylim)

                        fig.canvas.draw()
                else:
                    print("No points to delete. You must select at lease one point.")
                    
            ## Select points 
            elif event.key == '1':
                # Only allow selecting up to three points
                if point_counter < 3:  
                    # Save the coordinates of the clicked point
                    new_point = (event.xdata, event.ydata)
                    
                    if new_point in self.coords:
                        # Do not allow multiple selection of one point
                        print("Warning: The selected point already exists.")
                        print("Please select a new point.")      
                    else:
                        # Add selected point
                        self.coords.append(new_point)
    
                        # Visualize the selected point on the image
                        ax.plot(event.xdata, event.ydata, 'rx')
                        
                        # Retore the previous zoom level
                        ax.set_xlim(current_xlim)
                        ax.set_ylim(current_ylim)

                        fig.canvas.draw()
    
                        point_counter += 1
                else: # 3 points selected
                    print("3 points already selected. Press 'd' to continue.")
    
                if len(self.coords) == 3:
                    # Turn off interactive mode
                    plt.ioff()
    
                    print("3 points selected. Press 'd' to calculate the center position.")
    
            # Calculate circle or terminate
            elif event.key == 'd':
                if len(self.coords) == 3:
                    calculate_circle_flag = True
                else:
                    print("Please select exactly 3 points to calculate the circle.")
                    fig.canvas.draw()
    
        # Connect the callback function to the key press event
        fig.canvas.mpl_connect('key_press_event', onkeypress)
    
        # Show the plot
        # fig.set_size_inches(self.fig_width, self.fig_height)
        plt.tight_layout()
        plt.show(block=False)
    
        # Wait for 'd' key event or close the figure if no points are selected
        while not calculate_circle_flag and not termination_flag:
            try:
                plt.waitforbuttonpress(timeout=0.1)
            except KeyboardInterrupt:
                # If the user manually closes the figure, terminate the loop
                termination_flag = True
            except RuntimeError:
                # If the waitforbuttonpress times out, it means the figure 
                # is closed without any key press events
                termination_flag = True
                print('Detection unsuccessful. Execution terminated by closing the figure.')
                
    
        # If the termination_flag is True, stop the code
        if termination_flag: 
            print("No points selected. Returned None values.")
            return None, None, None
    
        # Calculate the circle if the points are selected
        if calculate_circle_flag:
            self.calculate_circle(plot_results)
            
            # Always adjust the center position after manual detection
            self.ref_interactive(self.x, self.y, self.r)
            
            return self.x, self.y, self.r

        
    def detection_intensity(self, csquare = 20, cintensity=0.5, plot_results = 0):
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
        arr = np.copy(self.image)
        xsize,ysize = arr.shape
        
        # Calculate borders around the central square
        xborder = (xsize - csquare) // 2
        yborder = (ysize - csquare) // 2
        
        # Create central square = cut off the borders
        arr2 = arr[xborder:-xborder,yborder:-yborder].copy()
        
        # In the central square, set all values below cintenstity to zero
        arr2 = np.where(arr2>np.max(arr2)*cintensity, arr2, 0)
        
        # Calculate 1st central moments of the image
        M = moments(arr2,1)
        
        # Calculate the intensity center = centroid according to www-help
        (self.x, self.y) = (M[1,0]/M[0,0], M[0,1]/M[0,0])
        
        # We have centroid of the central square => recalculate to whole image
        (self.x, self.y) = (self.x + xborder, self.y + yborder)
        self.r = 100
        
        # Print results
        print("--------- Diffraction pattern detection via intensity detection ----------")
        print("Central coordinate [ x, y ]: [{:.3f}, {:.3f}]".format(float(self.x), 
                                                                     float(self.y)))
        print("--------------------------------------------------------------------------")

        # Plot result of the Hough transform
        if plot_results == 1:
           self.visualize_center(self.x, self.y, self.r, tit ='Maximum intensity', view_enhancement=0) 
        
        # Return results
        return self.x, self.y, self.r

    
    def detection_Hough(self, canny_plot=0, plot_results=1 ):
        '''        
        Perform Hough transform to detect center of diffraction patterns.
        This is a method to automatically detect circular diffraction patterns
        
        Parameters
        ----------
        self.image : array of uint8
            Input image in which the diffraction pattern is to be found
        canny_plot : integer, binary
            Plot detected edges in input image. Default is 0 (do not plot)
        plot_results : int, binary
            Plot the pattern determined by pixels selected by the user.
            Default is 1. To cancel visualization, set plot_results = 0.

        Returns
        -------
        self.cx : float64
            x-coordinate of the detected center
        self.cy : float64
            y-coordinate of the detected center
        self.crad : float64
            radius of the detected center
                                    
        '''
        
        im = np.copy(self.image)
        # if the brightness of the image is small enough, pixel values greater
        # than 50 will be set to 0 -- removal of the beam stopper influence
        
        if self.heq == 0:
            if sum(sum(im)) < 150000:
                    max_indices = np.where(im > 50)
        
                    row_idx = max_indices[0]
                    col_idx = max_indices[1]
        
                    im[row_idx, col_idx] = 0    
                
            # Detect edges using the Canny edge detector
            edges = canny(im, 
                          sigma=0.2, 
                          low_threshold=80, high_threshold=100)
        elif self.heq == 1:
            if sum(sum(im)) > 150000:
                max_indices = np.where(im > 0.060)

                row_idx = max_indices[0]
                col_idx = max_indices[1]

                im[row_idx, col_idx] = 0    
            
            # Detect edges using the Canny edge detector
            edges = canny(im, 
                          sigma=0.2, 
                          low_threshold=0.80, high_threshold=1)
        
        # Plot detected edges if required
        if canny_plot==1:
            plt.figure()
            plt.imshow(edges)
            plt.title(label='Edges found using Canny detector')
            plt.show(block=False)
        
        
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
        
        # Print results
        print("------------ Diffraction pattern detection via Hough transform -----------")
        print("Central coordinate [ x, y ]: [{:.3f}, {:.3f}]".format(float(self.x), 
                                                                     float(self.y)))
        print("--------------------------------------------------------------------------")

        # Plot result of the Hough transform
        if plot_results==1:
           self.visualize_center(self.x, self.y, self.r, tit ='Hough transform', view_enhancement=1) 
        
        # Return results
        return self.x, self.y, self.r

   
    def calculate_circle(self, plot_results):
        ''' 
        Calculates coordinates of the center and radius of a circle defined via
        3 points determined by the user. Plots the calculated circle, detected 
        points and marks the center.
        
        Parameters
        ----------
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
        
        # Print results
        print(" ")
        print(" ")
        print("------------------- Manual diffraction pattern detection -----------------")
        print("Central coordinate [ x, y ]: [{:.3f}, {:.3f}]".format(float(self.x), 
                                                                     float(self.y)))
        print("--------------------------------------------------------------------------")

        if plot_results==1:
            # Create and manage the figure
            fig, ax = plt.subplots()
            manager = plt.get_current_fig_manager()
            manager.window.showMaximized()
            ax.imshow(self.image)
            
            # Plot center and points
            plt.plot(self.x, self.y, 
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
            ax.add_patch(circle)
            
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
        else:
            plt.close('all')
        
        return self.x, self.y, self.r


    def visualize_center(self, x, y, r, tit, view_enhancement = 0):
        '''         
        Visualize detected diffraction patterns and mark the center.
        
        Parameters
        ----------
        tit : string
            name of the method used for circle detection
        self.image : array of uint8
            Input image in which the diffraction pattern is to be found
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
        image = imread(self.image_path, as_gray=False)

        # Create figure and axes
        fig, ax = plt.subplots()


        # Draw circle
        circle = plt.Circle((np.median(x), np.median(y)), 
                            np.median(r), 
                            color='r', fill=False,
                            label = "Pattern")
        ax.add_patch(circle)

        # Plot center point
        ax.plot(np.median(x), 
                   np.median(y), 
                   color='r', 
                   marker='x', 
                   markersize=12,
                   label = 'Center')

        # Display the image
        ax.imshow(image)
        plt.title(tit)
        plt.legend(loc='lower center', 
                   ncol=2, 
                   bbox_to_anchor=(0.5,-0.1), 
                   mode='expand', 
                   frameon=False)        
        plt.axis('off')
        plt.tight_layout()
        plt.show(block=False)
        
        # Edit contrast with a user-predefined parameter
        if view_enhancement == 1:
            if self.icut != 0:
                print("works")
                # Display the original and equalized images
                plt.figure(figsize=(10, 5))
                plt.subplot(1, 2, 1)
                plt.imshow(image, 
                            cmap='gray')
                plt.title('Original Image')
    
                plt.subplot(1, 2, 2)
                plt.imshow(self.image, 
                            cmap='gray')
                plt.title('Enhanced Image, icut = {:.1f}'.format(self.icut))            
                plt.show(block=False)
            
            elif self.heq == 1:
                # Display the original and equalized images
                plt.figure(figsize=(10, 5))
                plt.subplot(1, 2, 1)
                plt.imshow(image, 
                            cmap='gray')
                plt.title('Original Image')
    
                plt.subplot(1, 2, 2)
                plt.imshow(self.image, 
                            cmap='gray')
                plt.title('Equalized Histogram')            
                plt.show(block=False)
                
         

    def central_square(self, arr, csquare, xcenter=None, ycenter=None):
        ''' 
        Return central square from an array
        
        Parameters
        ----------
        arr : 2D-numpy array
            The original array from which the central_square will be extracted
        csquare : int, optional, default is 20
            The size/edge of the square in the (geometrical) center.
            The intensity center will be searched only within the central square.
            Reasons: To avoid other spots/diffractions and
            to minimize the effect of possible intensity assymetry around center. 
        xcenter : float64
            x-coordinate of array center. Deafault is None
        ycenter : float64
            y-coordinate of array center. Deafault is None

        Returns
        -------
        arr2 : 2D-numpy array
            central square extracted from input array
        '''
        
        xsize, ysize = arr.shape
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
 
    

class CenterRefinement(CenterDetection):
    ''' SUBCLASS of CircleDetection
    ----------
    
    Automatic adjustment of the center position of diffraction patterns,
    which was found using methods defined in the class CircleDetection
    
    
    Parameters
    ----------
    image_path : string
        direct path to an image with diffraction patterns
    detection_method : string
        Selection of a method for center calculation. String codes are:
            - 'manual' : manual detection via 3 points
            - 'intensity' : detection via maximum intensity
            - 'hough' : automatic detection via Hough transform
    correction_method : string
       Selection of a method for center position correction. Default is None.
       String codes are:
           - 'manual' : manual corection 
           - 'variance' : correction via variance minimization 
           - 'sum' : correction via sum maximization
    heq : boolean
        Allow histogram equalization. The default is 0 (no enhancement)
    icut : boolean
        Allow image enhancement. The default is 0 (no enhancement)


    Returns
    -------
    self.xx : int32
        adjusted x-coordinate of the detected center
    self.yy : int32
        adjusted y-coordinate of the detected center
    self.rr : int32
        radius of the detected center
    '''
    
    def __init__(self, image_path, detection_method, heq = 0, icut = 0, correction_method = None):
        # Call the constructor of the base class to initialize its methods
        super().__init__(image_path, detection_method)
        self.heq = heq
        self.icut = icut

        # Enhance diffraction pattern to make it more visible
        if self.heq == 1:
           # print('a')
            self.image = equalize_adapthist(self.image)
            # plt.figure()
            # plt.imshow(self.image)
            # plt.show(block=False)
            
        # Edit contrast with a user-predefined parameter
        if icut != 0:
            self.image = np.where(self.image > self.icut, self.icut,  self.image)

            
            
            
        if correction_method is not None:
            self.ret = 1
            if correction_method == 'manual':
                self.xx, self.yy, self.rr = self.ref_interactive(self.x, self.y, self.r)
            elif correction_method == 'variance':
                self.xx, self.yy, self.rr = self.ref_minimize_var(self.x, self.y, self.r)
            elif correction_method == 'sum':
                self.xx, self.yy, self.rr = self.ref_maximize_sum(self.x, self.y, self.r)
        else:
            self.ret = 2
        
        


    def output(self):
        if self.ret == 1:
            return (self.x, self.y, self.xx, self.yy)  
        else:
            return (self.x, self.y, None, None)
        
    
    def ref_interactive(self, px, py, pr):
        ''' 
        Manual refinement of the detected diffraction pattern via one of 
        the methods provided in the class CircleDetection.
        
        The user can change the position of the center of the diffraction
        pattern and also the radius of the detected pattern using keys:
            
            - 'm' : while pressed down, MOVE MODE is ON 
            - left / right / top / down arrows : move left / right / top / down
            - '+' : increase radius
            - '-' : decrease radius
            - 'd' : done, termination of the refinement
        To navigate between views, use keys '<' / '>' to move back / forth

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
.

        Returns
        -------
        xy : tuple (float64, float64)
            new xy-coordinates of the center
        r : float64
            new radius of the circular diffraction pattern
        '''
        
        plt.close('all')
        # Load original image
        im = np.copy(self.image)
        
        # Initialize variables and flags
        xy = np.array((px, py))
        r = np.copy(pr)
        termination_flag = False
        global move_mode
        move_mode = False

   
        print(" ")
        print("--------------------------------------------------------------------------")
        print("Interactive refinement. Use these keys:")
        print("      - 'm' : press down to have the MOVE MODE ON. ")
        print("              do not press down to have the MOVE MODE OFF. ")
        print("      - 'left arrow' : move left")
        print("      - 'right arrow' : move right")
        print("      - 'top arrow' : move up")
        print("      - 'bottom arrow' : move down")
        print("      - '+' : increase circle radius")
        print("      - '-' : decrease circle radius")
        print("      - 'd' : refinement done")
        print("To navigate between views, use keys '<' / '>' to move back / forth. ")
        print("--------------------------------------------------------------------------")
        
        # Create a figure and display the image
        fig, ax = plt.subplots()
        plt.title("Press keys to adjust the center position")
        ax.imshow(im)
        
        # Enable interactive mode
        plt.ion()
        
        # Plot detected diffraction pattern
        circle = plt.Circle((px, py), pr, color='r', fill=False)
        ax.add_patch(circle)

        # Plot center point
        center, = ax.plot(px, py, 'rx', markersize=12)
        plt.title('Manually adjust the position of the center using keys.')

        # Display the image
        # fig.set_size_inches(self.fig_width, self.fig_height)
        plt.show(block=False)
        
        # Define the event handler for figure close event
        def onclose(event):
            nonlocal termination_flag
            termination_flag = True
            print('Execution terminated by user.')
            
        # Connect the event handler to the close event
        fig.canvas.mpl_connect('close_event', onclose)
    
        def on_key_release(event):
            global move_mode
            if event.key == 'm':
                move_mode = False
                #print('Move mode OFF')
            else:   
                try:
                    plt.rcParams['keymap.back'].remove('left')
                    plt.rcParams['keymap.forward'].remove('right')
                    plt.rcParams['keymap.back'].append(',')
                    plt.rcParams['keymap.forward'].append('.')
                except ValueError:
                    # Do nothing if 'left' or 'right' is not in the list
                    pass  

        # Connect the event handler to the key release event
        fig.canvas.mpl_connect('key_release_event', on_key_release)
        
        # Define the callback function for key press events
        def onkeypress(event):
            # Use nonlocal to modify the center position in the outer scope
            nonlocal xy, r, termination_flag
            global move_mode 

            # OTHER KEYS USED IN INTERACTIVE FIGURES
            #   event.key == '1': select a point in self.detection_3points()
            #   event.key == '2': delete the most recent point in self.detection_3points()
            #   event.key == '3': delete a point in self.detection_3points()
            #   event.key == 'd': proceed in self.detection_3points()
            
            if event.key == 'm':
                #print('Move mode ON')
                move_mode = not move_mode  # Toggle move_mode
        
            if move_mode:  # Only perform actions when move_mode is ON

                if event.key in ['up', 'down', 'left', 'right', '+', '-']:
                    if event.key in ['+', '-']:
                        r += 1 if event.key == '+' else -1
                     #   print('Radius', 'increased.' if event.key == '+' else 'decreased.')
                    else:
                        # Perform shifts normally
                        if event.key == 'up':
                            xy[1] -= 1
                           # print('Moved up')
                        elif event.key == 'down':
                            xy[1] += 1
                           # print('Moved down')
                        elif event.key == 'left':
                            xy[0] -= 1
                           # print('Moved left')
                        elif event.key == 'right':
                            xy[0] += 1
                           # print('Moved right')

            # Terminate the interactive refinement with 'd' key
            if event.key == 'd':
                    termination_flag = True
                    print("--------------------------------------------------------------------------")
                    print("Refinement done.")
                    print("--------------------------------------------------------------------------")

            # Update the plot with the new center position
            circle.set_center((xy[0], xy[1]))  # circle
            circle.set_radius(r)               # radius
            center.set_data([xy[0]], [xy[1]])  # center

            plt.title('Manually adjust the position of the center using keys.')
         
            # Update the plot
            plt.draw() 
        
        # Connect the callback function to the key press event
        fig.canvas.mpl_connect('key_press_event', onkeypress)

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
        # fig.set_size_inches(self.fig_width, self.fig_height)
        plt.tight_layout()
        plt.show(block=False)

        plt.close('all')

        ## Display original and refined images in one figure
        # self.visualize_refinement(px, py, pr, xy, r, view_enhancement=0)
              
        # Print results
        print(" ")
        print("--------------- Manual correction of radius and coordinates --------------")
        print("Central coordinate [ x, y ]: [{:.3f}, {:.3f}]".format(float(xy[0]), 
                                                                     float(xy[1])))
        print("--------------------------------------------------------------------------")
        
        
        return xy[0], xy[1], r

        
    def ref_minimize_var(self, px, py, pr, plot_results = 1):
        '''         
        Adjust center coordinates of a detected circular diffraction pattern.
        The center adjustment is based on variance minimization.

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
            Plot detected center of the diffraction pattern. The default is 1.
        
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
        
        # Load original image
        im = np.copy(self.image)
        
        # if the brightness of the image is small enough, pixel values greater
        # than 50 will be set to 0 -- removal of the beam stopper
        
        # CASE 1: the user want to enhance the image properties
        if self.heq != 0 or self.icut !=0:
            if sum(sum(im)) < 150000:
                    max_indices = np.where(im > 50)
                    row_idx = max_indices[0]
                    col_idx = max_indices[1]
        
                    im[row_idx, col_idx] = 0    
                
            # Detect edges using the Canny edge detector
            edges = canny(im, 
                          sigma=0.2, 
                          low_threshold=80, high_threshold=100)
        
        # CASE 2: the user does not want to enhance the image properties    
        else:
            if sum(sum(im)) > 150000:
                max_indices = np.where(im > 0.060)
                row_idx = max_indices[0]
                col_idx = max_indices[1]

                im[row_idx, col_idx] = 0    
            
            # Detect edges using the Canny edge detector
            edges = canny(im, 
                          sigma=0.9, 
                          low_threshold=0.80, high_threshold=1)
        
        
        # terminate condition of while cycle
        stopper = 0
        
        # initialization / calculation of variables
        px_coords, py_coords = self.get_circle_pixels(px, py, pr)
        px_coords = np.array(px_coords, dtype=int)
        py_coords = np.array(py_coords, dtype=int)
        
        # Filter pixel values based on the condition for pixel values
        filtered_values = edges[px_coords, py_coords][edges[px_coords, py_coords] >= 20]
        
        # Calculate variance using the filtered values
        pi_var = np.var(filtered_values)   # currenly the lowest variance
        
        # Coordinates of the input center
        min_position = (px, py)            # currently the circle center with pi_var 
        
        while stopper < 100:
            var_curr = []       # current variance of pixel values 
                                # (of the actual circle position)
            coords_curr = []    # current coordinates of the best center position 
                                # lowest variance (of the actual circle position)

        # Adjust center position based on minimazation of variance
        # the 8-neighbourhood pixels (x) of the current center (o) will be 
        # tested regarding the minimization:
        #   x x x     -->  (px - dx, py + dy) (px, py + dy) ( px + dx, py + dy)
        #   x o x     -->  (px - dx, py)      (px, py)       (px + dx, py)
        #   x x x     -->  (px - dx, py - dy) (px, py - dy)  (px + dx, py - dy)
            
            for dx in [-1, 0, 1]:                 # increment of x-coordinates
                for dy in [-1, 0, 1]:             # increment of y-coordinates
                    if dx == 0 and dy == 0:       # current center coordinates
                        continue
                    
                    # Select new center-candidate
                    nx, ny = px + dx, py + dy
                    
                    # Extract pixels on the circle border
                    pxc, pyc = self.get_circle_pixels(nx, ny, pr)
                    pxc = np.array(pxc, dtype=int)
                    pyc = np.array(pyc, dtype=int)
                    
                    # Filter pixel values based on the condition
                    fil_val = edges[pxc, pyc][edges[pxc, pyc] >= 20]
                    
                    # Calculate variance using the filtered values
                    v = np.var(fil_val)
                    
                    # Compare variance of the new circle with the current best 
                    if v < pi_var: 
                        var_curr.append(v)
                        coords_curr.append((nx, ny))
                        
                        # Find the lowest variance
                        w = np.argmin(var_curr)
                        
                        # Set the new lowest variance
                        pi_var = var_curr[w]
                        
                        # Set the new position of the center
                        min_position = coords_curr[w]
                        
                    else:   stopper += 1
            
            # Center position adjustment
            if (px, py) == min_position:
                px, py = px+random.choice([0, -1, 1]), py+random.choice([0, -1, 1])
            else: px, py = min_position
        
        plt.close('all')
        
        if plot_results == 1:
            self.visualize_refinement(bckup[0],  bckup[1], bckup[2], (px,py), pr, view_enhancement=0)

        # Print results
        print(" ")
        print("----- Correction of radius and coordinates via variance minimization -----")
        print("Central coordinate [ x, y ]: [{:.3f}, {:.3f}]".format(float(px), 
                                                                     float(py)))
        print("--------------------------------------------------------------------------")
        

        return px, py, pr
    
    
    def ref_maximize_sum(self, px, py, pr, plot_results=1):
        ''' 
        Adjust center position based on maximization of intensity sum.
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
            Plot detected center of the diffraction pattern. The default is 1.
    
        Returns
        -------
        px : float64
            Adjusted x-coordinate of the center.
        py : float64
            Adjusted y-coordinate of the center.
        pr : float64
            The adjusted radius of the circular diffraction pattern.
        '''
    
        # Store input for plot
        bckup = [np.copy(px), np.copy(py), np.copy(pr)]
    
        # Load original image
        # Load original image
        image = np.copy(self.image)
        image = canny(image, sigma=0.1, low_threshold=50, high_threshold=80)

        # initialization
        stopper = 0             # termination condition
        no_improv = 0           # no improvement counter
        max_intensity_sum = self.intensity_sum(image, px, py, pr)
        best_center = (px, py)
        best_radius = pr
    
        # Keep track of previous center and radius
        last_center = (px, py)
        last_radius = pr
    
        # Number of random starting points
        num_starting_points = 10
    
        # Termination condition
        while stopper < 100:
            neighbors = [(dx, dy) for dx in [-1, 0, 1] 
                                  for dy in [-1, 0, 1] 
                                  if dx != 0 or dy != 0]
    
            # Use multiple starting points
            for i in range(num_starting_points):
                px_start = px + random.choice([*range(-1,1)])
                py_start = py + random.choice([*range(-1,1)])
                # print((px_start, py_start))
    
                for dx, dy in neighbors:
                    nx, ny = px_start + dx, py_start + dy
    
                    # Calculate the intensity sum for the new candidate center
                    curr_intensity_sum = self.intensity_sum(image, nx, ny, pr)
    
                    if curr_intensity_sum > max_intensity_sum:
                        max_intensity_sum = curr_intensity_sum
                        best_center = (nx, ny)
    
                # Adjust the radius
                for dr in [-1, -0.5, 0, 0.5, 1]:
                    new_radius = pr + dr
    
                    # Calculate the intensity sum for the current center 
                    # and new candidate radius
                    curr_intensity_sum = self.intensity_sum(image, 
                                                            px_start, 
                                                            py_start, 
                                                            new_radius)
    
                    if curr_intensity_sum > max_intensity_sum:
                        max_intensity_sum = curr_intensity_sum
                        best_center = (px_start, py_start)
                        best_radius = new_radius
    
            # Termination condition:
            # check if the current best center and radius are the same
            # as the previous ones
            if best_center == last_center and best_radius == last_radius:
                no_improv += 1
                #  print(no_improv)
                px = px + random.choice([0, -1, 1])
                py = py + random.choice([0, -1, 1])
                if no_improv == 20:
                    break
    
            last_center = best_center
            last_radius = best_radius
            px, py = best_center
            pr = best_radius
    
            stopper += 1
    
        plt.close('all')
        
        if plot_results == 1:
            self.visualize_refinement(bckup[0],  bckup[1], bckup[2], 
                                      (px,py), pr, view_enhancement=0)

        # Print results
        print(" ")
        print("--- Correction of radius and coordinates via intensity sum maximization --")
        print("Central coordinate [ x, y ]: [{:.3f}, {:.3f}]".format(float(px), 
                                                                     float(py)))
        print("--------------------------------------------------------------------------")
                
        return px, py, pr
    
    
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
        
        # Calculate x,y-coordinates of points
        x = xc + radius * np.cos(theta)  
        y = yc + radius * np.sin(theta)
        
        return x, y
          
    
    def intensity_sum(self, image, px, py, pr):
        ''' 
        Summation of intensity values of pixels definifing a diffraction pattern.

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
        s = np.sum(image[pxc, pyc])
        return s
    
    
    def visualize_refinement(self, px, py, pr, xy, r, view_enhancement = 0):
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
        view_enhancement : boolean
            Plot original image and enhanced for visual comparison. The default
            is 0 (do not plot)

        Returns
        -------
        None.

        '''
        # Load original image
        image = imread(self.image_path, as_gray = False)
        
        # Display original and refined images in one figure
        fig, ax = plt.subplots(nrows=1, ncols=2)

        
        ax[0].imshow(image)
        c0 = plt.Circle((px, py), pr, 
                        color='r', 
                        fill=False,
                        label = 'pattern')
        ax[0].add_patch(c0)
        ax[0].scatter(px,py, 
                      label='center', 
                      color='r', 
                      marker='x', 
                      s=100)
        ax[0].set_title('Detected center of the diffraction pattern.')
        ax[0].legend(loc='lower center', 
                   ncol=2, 
                   bbox_to_anchor=(0.5,-0.1), 
                   mode='expand', 
                   frameon=False) 
        ax[0].axis('off')
        
        ax[1].imshow(image)
        c1 = plt.Circle(xy, r,
                        color='r', 
                        fill=False,
                        label = 'pattern')
        ax[1].add_patch(c1)
        ax[1].scatter(xy[0],xy[1], 
                      label='center', 
                      color='r', 
                      marker='x',
                      s=100)
        ax[1].set_title('Corrected center position.')
        ax[1].legend(loc='lower center', 
                   ncol=2, 
                   bbox_to_anchor=(0.5,-0.1), 
                   mode='expand', 
                   frameon=False)  
        ax[1].axis('off')
        plt.tight_layout()
        plt.show(block=False)
        
        # Edit contrast with a user-predefined parameter
        if view_enhancement == 1:
            if self.icut != 0:
                print("works")
                # Display the original and equalized images
                plt.figure(figsize=(10, 5))
                plt.subplot(1, 2, 1)
                plt.imshow(image, 
                            cmap='gray')
                plt.title('Original Image')
    
                plt.subplot(1, 2, 2)
                plt.imshow(self.image, 
                            cmap='gray')
                plt.title('Enhanced Image, icut = {:.1f}'.format(self.icut))            
                plt.show(block=False)
            
            elif self.heq == 1:
                # Display the original and equalized images
                plt.figure(figsize=(10, 5))
                plt.subplot(1, 2, 1)
                plt.imshow(image, 
                            cmap='gray')
                plt.title('Original Image')
    
                plt.subplot(1, 2, 2)
                plt.imshow(self.image, 
                            cmap='gray')
                plt.title('Equalized Histogram')            
                plt.show(block=False)

        
#%% OLD SIMPLE FUNCTIONS

# Old interface, kept just for backward compatibility
# To be removed in one of the next versions

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
