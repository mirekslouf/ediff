'''
Module ediff.center
-------------------
Find center of 2D diffraction pattern. 

CenterDet, aktualizace 06-10-23 CentDet update, methods compatibility
'''

import numpy as np
import skimage as sk
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.legend_handler import HandlerBase

import ediff.io

from skimage.measure import moments
from skimage.transform import hough_circle, hough_circle_peaks


import warnings
import sys
warnings.filterwarnings("ignore")


class CenterEstimator:
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
    heq : boolean, optional, default is False
        Allow histogram equalization.
        The equalization is virtual, the original image is unchanged.
    icut : integer, optional, default is None
        Intensity cut to enhance image contrast.
        Example: icut=300, all pixels intensities >300 are set to 300.
    cmap : str, optional, default is 'gray'
        Name of matplotlib colormap.
        Any valid colormap ('viridis','magma') can be used.
    messages : bool, optional, default is False
        Print to terminal additional messages about program run.
    
    Returns
    -------    
    self.x : float64
        x-coordinate of the detected center
    self.y : float64
        y-coordinate of the detected center
    self.r : float64
        radius of the detected center (if available, othervise returns None)
        
    
    tady jsem udelala zmenu
                    
    '''
    
    def __init__(self, input_image, 
                       detection_method = None, 
                       correction_method = None,
                       heq = False, 
                       icut = None,
                       cmap = 'gray',
                       messages = False):
        
        # Initialize attributes
        self.input_image = input_image
        
        # Allow input images (np.ndarray) and image path
        if isinstance(input_image, np.ndarray):
            self.image = input_image
        else:
            self.image = ediff.io.read_image(self.input_image)
            
        self.correction_method = correction_method
        self.heq = heq
        self.icut = icut
        self.to_refine = []
        self.cmap = cmap
        self.messages = messages
        
        # Initialize the marker size and zoom level
        self.marker_size = 100  # Adjust this initial value as desired
        
        # Run functions
        self.preprocess_images(preInit = 1)

        # Determine detection method
        if detection_method == 'manual':
            self.x, self.y, self.r = self.detection_3points()
            
        elif detection_method == 'intensity':
            self.x, self.y, self.r = self.detection_intensity()

        elif detection_method == 'hough':
            self.x, self.y, self.r = self.detection_Hough()  
        else:
            print("Incorrect method for detection selected")
            sys.exit()
            
                
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
        im = self.to_refine
    
        # Create a figure and display the image
        fig, ax = plt.subplots(figsize=(12, 12))
        
        # Allow using arrows to move back and forth between view ports
        plt.rcParams['keymap.back'].append('left')
        plt.rcParams['keymap.forward'].append('right')
 
        plt.title("Select 3 points defining one of diffraction circles", 
                  fontsize = 20)
        ax.imshow(im, cmap = self.cmap)
        ax.axis('off')

 
        # User information:
        if self.messages:
            print("--- Manual diffraction pattern detection ---")
            print()
            print("Select 3 points to define a circle in the diffractogram:")
            print("Use these keys for the selection:")
            print("  - '1' : define a point at current cursor position")
            print("  - '2' : delete the last point")
            print("  - '3' : delete the point closest to the cursor")
            print("  - 'd' : done = finished = go to the next step")
            print()
            print("Close the figure to terminate. No center will be detected.")
       
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
            if self.messages:
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
                    distances = [
                        np.sqrt((x - pointer_x)**2 + (y - pointer_y)**2)
                        for x, y in self.coords]
                    closest_index = np.argmin(distances)
                    del self.coords[closest_index]
    
                    # Redraw the image without the deleted point
                    ax.clear()
                    ax.imshow(self.to_refine, cmap = self.cmap)
                    for x, y in self.coords:

                        ax.scatter(x, y, 
                                   c='r', marker='x', 
                                   s=self.marker_size)

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
    
            ## Delete recent point (last added) -- independent on the cursor
            if event.key == '2':
                # Check if there are points to delete
                if point_counter > 0:  
                    point_counter -= 1
                    if len(self.coords) > 0:
                        # Delete the last point in the list
                        del self.coords[-1]
                        if self.messages:
                            print('The most recently selected point deleted.')
                            print('Please select a new one.')
    
                        # Redraw the image without the deleted point
                        ax.clear()
                        ax.imshow(self.to_refine, cmap = self.cmap)
                        for x, y in self.coords:
                            ax.scatter(x, y,
                                       c='r', marker='x', 
                                       s=self.marker_size)

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
                                   s=self.marker_size)

                        # Retore the previous zoom level
                        ax.set_xlim(current_xlim)
                        ax.set_ylim(current_ylim)
                        ax.axis('off')

                        fig.canvas.draw()
    
                        point_counter += 1
                else: # 3 points selected
                    print("3 points already selected.")
    
                if len(self.coords) == 3:
                    # Turn off interactive mode
                    plt.ioff()
    
                    if self.messages:
                        print("3 points selected.")
                        print("Press 'd' to calculate the center position.")
    
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
                    ax.imshow(self.to_refine, cmap = self.cmap)
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
        
        
        # Manually adjust the calculated center coordinates
        self.x, self.y, self.r = self.adjustment_3points(fig, circle, center)

 
        return self.x, self.y, self.r

    
    def adjustment_3points(self, fig, circle, center, plot_results=0):
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
        
        if self.messages:
            print("--- Interactive center refinement ---")
            print("Use these keys:")
            print("  - 'leftArrow'  : move left")
            print("  - 'RightArrow' : move right")
            print("  - 'UpArrow'    : move up")
            print("  - 'DownArrow'  : move down")
            print("  - '+' : increase circle radius")
            print("  - '-' : decrease circle radius")
            print("  - 'd' : refinement done")
            print("Technical note:")
            print("  the default shortcuts for LeftArrow and RightArrow")
            print("  temporarily overriden, but corresponding icons do work.")
        
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
            if self.messages:
                print('Execution terminated.')
 
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
            #   event.key == 'd': proceed in self.detection_3points()
            
            if event.key in ['up', 'down', 'left', 'right', '+', '-']:
                if event.key in ['+', '-']:
                    r += 1 if event.key == '+' else -1
                else:
                    # Perform shifts normally
                    if event.key == 'up':
                        xy[1] -= 0.5
                       # print('Moved up')
                    elif event.key == 'down':
                        xy[1] += 0.5
                       # print('Moved down')
                    elif event.key == 'left':
                        xy[0] -= 0.5
                       # print('Moved left')
                    elif event.key == 'right':
                        xy[0] += 0.5
                       # print('Moved right')

            # Terminate the interactive refinement with 'd' key
            if event.key == 'd':
                termination_flag = True
                if self.messages:
                    print("--- Refinement done. ---")

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

        # Print results
        if self.messages:
            print("CenterEstimator :: manual detection + adjustment")
            print(f"Center coordinates: {xy[0]:.2f} {xy[1]:.2f}")
        
    
        return xy[0], xy[1], r
        
        
    def detection_intensity(self, csquare=20, cintensity=0.5, plot_results=0):
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
        image = np.copy(self.to_refine)
        arr = np.copy(image)
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
        
        # User information:
        if self.messages:
            print("--------- Diffraction pattern detection via intensity detection ----------")
            print("Central coordinate [ x, y ]: [{:.3f}, {:.3f}]".format(float(self.x), 
                                                                         float(self.y)))
            print("--------------------------------------------------------------------------")

        # Plot result of the Hough transform
        if plot_results == 1:
           self.visualize_center(self.x, self.y, self.r) 
        
        
        # Return results
        return self.x, self.y, self.r

    
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
        im = np.copy(self.to_refine)
        
        # if the brightness of the image is small enough, pixel values greater
        # than 50 will be set to 0 -- removal of the beam stopper influence
        
        if self.heq == 0:
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
        if self.messages:
            print("------------ Diffraction pattern detection via Hough transform -----------")
            print("Central coordinate [ x, y ]: [{:.3f}, {:.3f}]".format(float(self.x), 
                                                                         float(self.y)))
            print("--------------------------------------------------------------------------")
        
        
        self.x, self.y, self.r = self.x[0], self.y[0], self.r[0]
        # Return results, convert coordinates to float
        return self.x, self.y, self.r

   
    def calculate_circle(self, plot_results):
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
        
        # Print results
        if self.messages:
            print("CenterEstimator :: manual center detection")
            print(f"Center coordinates: {self.x:.2f} {self.y:.2f}")
                    
        if plot_results==1:
            # Create and manage the figure
            fig, ax = plt.subplots()
            manager = plt.get_current_fig_manager()
            manager.window.showMaximized()
            ax.imshow(self.image, cmap = self.cmap)
            
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
        else:
            pass
            #plt.close('all')
        
        self.center = (self.x, self.y)
        self.circle = plt.Circle((self.x,self.y),self.r)
        

        return self.x, self.y, self.r, self.center, self.circle


    def visualize_center(self, x, y, r):
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
        im = self.to_refine
    
        # Create a figure and display the image
        fig, ax = plt.subplots(figsize=(12, 12))
        
        # Allow using arrows to move back and forth between view ports
        plt.rcParams['keymap.back'].append('left')
        plt.rcParams['keymap.forward'].append('right')
 
        plt.title("Detected center", 
                  fontsize = 20)
        ax.axis('off')

        # Plot center point
        ax.scatter(x,y,
                label= f'center:  [{x:.1f}, {y:.1f}]',
                marker='rx', s=100)

        plt.legend(loc='upper right', fontsize=20)
        
        # Display the image
        ax.imshow(im, cmap = self.cmap)
        plt.axis('off')
        plt.tight_layout()
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
 
    

class CenterLocator(CenterEstimator):
    ''' SUBCLASS of CircleEstimator
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
    
    def __init__(self, image_path,
                 detection_method, correction_method=None, 
                 heq=False, icut=None, cmap='gray',
                 messages=False, final_replot=False):
        
        # (1) Call the constructor of the base class to initialize its methods
        super().__init__(image_path,
            detection_method, correction_method,
            heq=heq, icut=icut, cmap=cmap, messages=messages)
        

        # (2) Define additional parameter
        self.final_replot = final_replot
        self.messages = messages
        self.image_path = image_path


        # (3) Run correction method and get refined parameters
        if correction_method is not None:
            self.ret = 1
            if correction_method == 'manual':
                if detection_method == 'manual':
                    self.xx, self.yy, self.rr = self.x, self.y, self.r
                    self.x, self.y, self.r = self.backip[0], self.backip[1], self.r
                    if self.final_replot:
                        self.visualize_refinement(
                            self.backip[0], self.backip[1], self.r, 
                            (self.xx, self.yy), self.r)
                else:
                    self.xx, self.yy, self.rr = \
                        self.ref_interactive(self.x, self.y, self.r)
            elif correction_method == 'variance':
                self.xx, self.yy, self.rr = \
                    self.ref_minimize_var(self.x, self.y, self.r)
            elif correction_method == 'sum':
                self.xx, self.yy, self.rr = \
                    self.ref_maximize_sum(self.x, self.y, self.r)
            else:
                print("Incorrect method for correction selected")
                sys.exit()

        else:
            self.ret = 2


    def preprocess_images(self, preInit=0, preHough=0, preManual=0, preVar = 0, 
                          preSum = 0, preInt=0):
        """
        FOR AUTOMATIC METHODS OPTIMIZATION AND MORE UNIVERSAL SOLUTIONS
        >>> to be finished later, as it is not the most important thing now
            
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
        image = np.copy(self.image)
        
        ### After initialization: perform an image enhancement if specified
        if preInit == 1:
            # Enhance diffraction pattern to make it more visible
            if self.heq == 1:
                if self.messages:
                    print("Histogram equalized.")
                image = sk.exposure.equalize_adapthist(image)
                # plt.figure()
                # plt.imshow(image)
                # plt.show(block=False)
                
            # Edit contrast with a user-predefined parameter
            if self.icut is not None:
                if self.messages:
                    print("Contrast enhanced.")
                image = np.where(image > self.icut, self.icut, image)
                

            self.to_refine = image
            return
        
        #######################################################################
        # Hough transform: perform pre-processing necessary for the detection
        if preHough == 1:
           # self.imHough = np.copy(self.to_refine)
           
            if self.heq == 0:
                csq = self.central_square(self.to_refine, csquare=80)   
                
                # Beam stopper present in image
                if np.median(csq)<100 and np.median(csq) > 0:
                    if self.messages:
                        print('Beamstopper removed.')
                    max_indices = np.where(self.to_refine > np.median(csq))
        
                    row_idx = max_indices[0]
                    col_idx = max_indices[1]
        
                    self.to_refine[row_idx, col_idx] = 0    
                    
                    max_indices = np.where(self.to_refine < 0.8*np.median(csq))
                    row_idx = max_indices[0]
                    col_idx = max_indices[1]
        
                    self.to_refine[row_idx, col_idx] = 0   
                    
                    # Detect edges using the Canny edge detector
                    edges = sk.feature.canny(self.to_refine, 
                                              sigma=0.2, 
                                              low_threshold=2.5*np.median(self.to_refine), 
                                              high_threshold=3*np.median(self.to_refine))
                    
                    # Dilate the edges to connect them
                    selem = sk.morphology.disk(5)
                    dilated_edges = sk.morphology.dilation(edges, selem)
                    
                    # Erode the dilated edges to reduce thickness and smooth the contour
                    connected_edges = sk.morphology.erosion(dilated_edges, selem)

                    if control_print == 1:
                        fig, ax = plt.subplots(nrows=2, ncols=2)
                        ax[0,0].imshow(self.image)
                        ax[0,0].set_title("Original image")
                        ax[0,1].imshow(self.to_refine)
                        ax[0,1].set_title("Hough pre-processed")
                        ax[1,0].imshow(edges)
                        ax[1,0].set_title("Edges")
                        ax[1,1].imshow(connected_edges)
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
                    
                    # # Erode the dilated edges to reduce thickness and smooth the contour
                    connected_edges = sk.morphology.erosion(dilated_edges, selem)
                    connected_edges = sk.morphology.remove_small_objects(connected_edges, 
                                                                         min_size=100)
                    
                    if control_print == 1:
                        fig, ax = plt.subplots(nrows=2, ncols=2)
                        ax[0,0].imshow(self.image)
                        ax[0,0].set_title("Original image")
                        ax[0,1].imshow(self.to_refine)
                        ax[0,1].set_title("Hough pre-processed")
                        ax[1,0].imshow(edges)
                        ax[1,0].set_title("Edges")
                        ax[1,1].imshow(connected_edges)
                        ax[1,1].set_title("Connected edges")
                        plt.tight_layout()
                        plt.show(block=False)
                
            elif self.heq == 1: 
                # Central square extraction
                csq = self.central_square(self.to_refine, csquare=80)   

                # Beam stopper present in image
                if 0.4 <= np.median(csq) <= 0.6:
                    
                    max_indices = np.where(self.to_refine > 2*np.median(self.to_refine))
    
                    row_idx = max_indices[0]
                    col_idx = max_indices[1]
    
                    self.to_refine[row_idx, col_idx] = 0 
                                                     
                    # Detect edges using the Canny edge detector
                    edges = sk.feature.canny(self.to_refine, 
                                              sigma=0.2, 
                                              low_threshold=1.5*np.median(self.to_refine), 
                                              high_threshold=3*np.median(self.to_refine))

                    
                    # Erode the dilated edges to reduce thickness and smooth the contour
                    #  connected_edges = sk.morphology.erosion(dilated_edges, selem)
                    
                    if control_print == 1:
                        fig, ax = plt.subplots(nrows=2, ncols=2)
                        ax[0,0].imshow(self.image)
                        ax[0,0].set_title("Original image")
                        ax[0,1].imshow(self.to_refine)
                        ax[0,1].set_title("Hough pre-processed")
                        ax[1,0].imshow(edges)
                        ax[1,0].set_title("Edges")
                      #  ax[1,1].imshow(connected_edges)
                        ax[1,1].set_title("Connected edges")
                        plt.tight_layout()
                        plt.show(block=False)
                    
                # No beam stopper in image
                else:
                    # Detect edges using the Canny edge detector
                    edges = sk.feature.canny(self.to_refine, 
                                              sigma=0.2, 
                                              low_threshold=2.5*np.median(self.to_refine), 
                                              high_threshold=3*np.median(self.to_refine))

                   # connected_edges = sk.morphology.erosion(dilated_edges, selem)
                    
                    if control_print == 1:
                        fig, ax = plt.subplots(nrows=2, ncols=2)
                        ax[0,0].imshow(self.image)
                        ax[0,0].set_title("Original image")
                        ax[0,1].imshow(self.to_refine)
                        ax[0,1].set_title("Hough pre-processed")
                        ax[1,0].imshow(edges)
                        ax[1,0].set_title("Edges")
                     #   ax[1,1].imshow(connected_edges)
                        ax[1,1].set_title("Connected edges")
                        plt.tight_layout()
                        plt.show(block=False)
            
            return edges
        
                    
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
            x-coordinate of the center detected via correction_method
        yy : float
           y-coordinate of the center detected via detection_method
           
        """
        
        if self.ret == 1:
            # Convert to float
            if type(self.x) != float:
                self.x, self.y, self.xx, self.yy, self.r, self.rr = \
                    [float(value) for value in (self.x, self.y, 
                                                self.xx, self.yy, 
                                                self.r, self.rr)]
            
            # Round radius
            self.r, self.rr = np.round(self.r), np.round(self.rr)
            
            # Return values of center coordinates
            return (np.round(self.x,1), np.round(self.y,1), 
                    np.round(self.xx,1), np.round(self.yy,1))  
        else:
            # Convert to float
            if type(self.x) != float:
                self.x, self.y, self.r = \
                    [float(value) for value in (self.x, self.y, self.r)]

            # Return values of center coordinates
            return (np.round(self.x,1), np.round(self.y,1), None, None)
        
    
    def ref_interactive(self, px, py, pr):
        ''' 
        Manual refinement of the detected diffraction pattern via one of 
        the methods provided in the class CircleDetection.
        
        The user can change the position of the center of the diffraction
        pattern and also the radius of the detected pattern using keys:
            - left / right / top / down arrows : move left / right / top / down
            - '+' : increase radius
            - '-' : decrease radius
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
        im = np.copy(self.to_refine)
        
        # Initialize variables and flags
        xy = np.array((px, py))
        r = np.copy(pr)
        termination_flag = False

        # User information:
        if self.messages:
            print(" ")
            print("--------------------------------------------------------------------------")
            print("Interactive refinement. Use these keys:")
            print("      - 'left arrow' : move left")
            print("      - 'right arrow' : move right")
            print("      - 'top arrow' : move up")
            print("      - 'bottom arrow' : move down")
            print("      - '+' : increase circle radius")
            print("      - '-' : decrease circle radius")
            print("      - 'd' : refinement done")
            print("DISCLAIMER: for the purpose of the center shift, the default shortcuts ")
            print("for left and right arrows were removed.")
            print("--------------------------------------------------------------------------")
        
        # Create a figure and display the image
        fig, ax = plt.subplots(figsize=(12, 12))
        
        # Allow using arrows to move back and forth between view ports
        plt.rcParams['keymap.back'].append('left')
        plt.rcParams['keymap.forward'].append('right')
        circle = plt.Circle(
            (px, py), pr, color='r', fill=False)
        ax.add_artist(circle)

        # Plot center point
        center, = ax.plot(self.x, self.y, 'rx', markersize=12)
                    

        plt.title('Manually adjust the position of the center using keys.', 
                  fontsize=20)

        ax.imshow(im, cmap = self.cmap)
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
            if self.messages:
                print('Execution terminated.')
 
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
            #   event.key == 'd': proceed in self.detection_3points()
            
            if event.key in ['up', 'down', 'left', 'right', '+', '-']:
                if event.key in ['+', '-']:
                    r += 1 if event.key == '+' else -1
                else:
                    # Perform shifts normally
                    if event.key == 'up':
                        xy[1] -= 0.5
                       # print('Moved up')
                    elif event.key == 'down':
                        xy[1] += 0.5
                       # print('Moved down')
                    elif event.key == 'left':
                        xy[0] -= 0.5
                       # print('Moved left')
                    elif event.key == 'right':
                        xy[0] += 0.5
                       # print('Moved right')

            # Terminate the interactive refinement with 'd' key
            if event.key == 'd':
                termination_flag = True
                if self.messages:
                    print("--- Refinement done. ---")

            # Update the plot with the new center position
            circle.set_center((xy[0], xy[1]))  # circle
            circle.set_radius(r)               # radius
            center.set_data([xy[0]], [xy[1]])  # center

            plt.title('Manually adjust the position of the center using keys.', 
                      fontsize=20)
         
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

        # Print results
        if self.messages:
            print("CenterEstimator :: manual detection + adjustment")
            print(f"Center coordinates: {xy[0]:.2f} {xy[1]:.2f}")
        
        # Plot results
        if self.final_replot:
            self.visualize_refinement(px, py, pr, xy, r)
                
        return xy[0], xy[1], r

        
    def ref_minimize_var(self, px, py, pr, plot_results = 0):
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
        image = np.copy(self.image)

        # Starting values to be modified 
        init_var = self.intensity_var(image, px, py, pr)
        min_intensity_var = self.intensity_var(image, px, py, pr)
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
            curr = self.intensity_var(image, 
                                      best_center[0], 
                                      best_center[1], 
                                      best_radius) 
            # Store intensity sums of the current center's neighborhood
            curr_intensity_var = []
            for dx, dy in neighbors:
                nx, ny = best_center[0] + dx, best_center[1] + dy
                # Check if the point is within the expanded search radius
                curr_intensity_var.append(self.intensity_var(image, 
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
            min_intensity_var = self.intensity_var(image, 
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
                radi_intensity_var.append(self.intensity_var(image, 
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

        # Refinement visualization
        if self.final_replot:
            self.visualize_refinement(
                   bckup[0], bckup[1], bckup[2], (best_center), best_radius)
    
        # Print results
        if self.messages:
            print("CenterLocator: manual detection + adjustment:")
            print(f"Estimated center {px:.2f} {py:.2f}")
                
        return best_center[0], best_center[1], best_radius
    
    
    def ref_maximize_sum(self, px, py, pr, plot_results=1):
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
        image = np.copy(self.image)

        # Starting values to be modified 
        init_sum = self.intensity_sum(image, px, py, pr)
        max_intensity_sum = self.intensity_sum(image, px, py, pr)
        best_center = (np.copy(px), np.copy(py))
        best_radius = np.copy(pr)
        
        # Convergence criterion for termination of gradient optimization 
        # (1) small positive value that serves as a threshold to determine 
        #     when the optimization process has converged
        convergence_threshold = 0.05*max_intensity_sum
        
        # (2) maximum number of iterations of optimization
        max_iterations = 50
        
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
            curr = self.intensity_sum(image, 
                                      best_center[0], 
                                      best_center[1], 
                                      best_radius)
            
            # Store intensity sums of the current center's neighborhood
            curr_intensity_sum = []
            for dx, dy in neighbors:
                nx, ny = best_center[0] + dx, best_center[1] + dy
                # Check if the point is within the expanded search radius
                curr_intensity_sum.append(self.intensity_sum(image, 
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
            max_intensity_sum = self.intensity_sum(image, 
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
                radi_intensity_sum.append(self.intensity_sum(image, 
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
        ## (1) swapped coordinates
        # if ((bckup[0] > bckup[1] and not best_center[0] > best_center[1])
        #     or  (bckup[0] < bckup[1] and not best_center[0] < best_center[1])):
        #     best_center = best_center[::-1]
        
        ## (2) worsened final maximum intensity sum than the initial one
        if np.round(init_sum,-2) > np.round(max_intensity_sum,-2):
            print("Refinement redundant.")
            best_center = np.copy(bckup)

        # Refinement visualization
        if plot_results == 1:
            if self.final_replot:
                self.visualize_refinement(
                    bckup[0], bckup[1], bckup[2], (best_center), best_radius)
    
        # Print results
        if self.messages:
            print("CenterLocator: manual detection + adjustment:")
            print(f"Estimated center {px:.2f} {py:.2f}")
                
        return best_center[0], best_center[1], best_radius
    
    
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

        fig, ax = plt.subplots(figsize=(12, 12))

        # Original Image
        ax.imshow(image, cmap=self.cmap)
        c0 = plt.Circle((px, py), pr,
                        color='r',
                        fill=False,
                        label='detected',
                        linewidth=1)
        ax.add_patch(c0)
        ax.scatter(px, py,
                   label= f'd-center: [{px:.1f}, {py:.1f}]',
                   color='r',
                   marker='x',
                   s=100, linewidths=1)

        # Refined Image
        c1 = plt.Circle(xy, r,
                        color='springgreen',
                        fill=False,
                        label='refined',
                        linewidth=1)
        ax.add_patch(c1)
        ax.scatter(xy[0], xy[1],
                   label= f'r-center:  [{xy[0]:.1f}, {xy[1]:.1f}]',
                   color='springgreen',
                   marker='x',
                   linewidths=1, s=100)

        ax.set_title('Center Detection and Correction', fontsize=20)
        ax.legend(loc='upper left', frameon=True, fontsize=20, 
                  handler_map={Circle: HandlerCircle()}, ncol=2,
                  bbox_to_anchor=(0.085, 1.00), bbox_transform=ax.transAxes)  
        ax.axis('off')

        plt.tight_layout()
        plt.show(block=False)

class HandlerCircle(HandlerBase):
    '''
    Help class for visualization. Not to be used outside of the module.
    '''
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        # Create a circular marker
        x = width / 2
        y = height / 2
        r = min(width, height) / 2
        marker = Circle((x, y), r, facecolor=orig_handle.get_facecolor(),
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
            
        Note
        ----
        Testovaci zmena.
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
        # IMPORTANT: switch xc,yc => arrays(rows,cols) x images(cols,rows
        (xc,yc) = (yc,xc)
        # Return the final center
        return(xc,yc)