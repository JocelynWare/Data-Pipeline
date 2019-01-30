'''
Created on 29 Jan 2019

@author: Jacob Parnell

'''

## documented changes in reduction script ##
## from basic veloce reduction script to MQ single-fibre adaptation ##
# as of 29 Jan 2019
# will run quick extraction method, will not run optimal extraction

# (1) basic_reduction_script.py:
#       Commented and tidied up for presentation
# (2) line 346, 351, 607 and 613 in extraction.py:
#       The path for variable 'fibparms' under the 'if simu:' condition now set to our directory
# (3) line 256 in process_scripts.py:
#       Added passed parameter 'simu=True', need to change to False for real data
# (4) line 737/738 in calibration.py:
#       Commented out 'TOTALEXP' and uncommented 'exptime', must check if real data have these headers
# (5) line 157 in basic_reduction_script:
#       Copied step four but code can be run for optimal extraction via change of parameter 'ext_method'
# (6) line 155 in basic_reduction_script:
#       Changed parameter 'fibs' to 'single' for our single fibre veloce adaptation
# (7) lines 68-79 in calibration.py:
#       Defined function altered to include int() as error raises in Python 3, not in Python 2
# (8) ***Need to delete created files before re-running code, otherwise 'stellar_list' imports new files too***
#       Alternatively, could rewrite the code to save in a new path
# (9) lines 52/53, 217/218, 274/275, 359/360, 621/622, 974/975, 1141/1142 and 1148/1149 in extraction.py:
#       Added lines of code to skip over 'order_01' - this is as issue arose with veloce required a 'fix' which
#       is due to the order not falling on the chip entirely, and may not be an issue with MQLTS, hence skipped
# (10) lines 554/555, 627/628 and 693/694 in order_tracing.py:
#       Added lines of code to skip over 'order_01' - this is as issue arose with veloce required a 'fix' which
#       is due to the order not falling on the chip entirely, and may not be an issue with MQLTS, hence skipped
# (11) plot_test.py:
#       Added a new .py file that contains code to convert the 2-dim extracted images to 1-dim spectra. Refined it
#       so it will iterate over all elements in a list of all the extracted stellar 2-dim images, and also iterate
#       over all orders in each of those files.
# (12) lines 112 and 121 in basic_reduction_script.py:
#       Added in those lines of code so that we can run plot_test.py without having to run the entire
#       basic_reduction_script. Outputs the flux in a dictionary as a .npy file.


## brief explanation of reduction process - step by step ##
## basic_reduction_script.py ##
## quick method only ##

## PREP
# lines 9-14:
#   Import important python functions that are useful in this reduction
# lines 16-21:
#   Import/load other sub-routines (as .py files) which can be called to perform different tasks. These files must be
#   in the same directory as the reduction script

## STEP (0) - GET FILE INFO
# lines 32-34:
#   Imports either the simulated data, OR real data, from the directory that you have designated as 'path'
# lines 37 and 38:
#   Generates a dummy image of the first file from a list that you had just generated. The purpose is so that for our
#   tests on simulated data, we are able to generate a master bias and master dark of the same shape

## STEP (1) - BAD PIXEL MASK
# Do not have code here for our reduction, but similar code can be adapted from full veloce code

## STEP (2) - CALIBRATION
# line 50:
#   Set the gain to be a 4x1 list of ones ([1,1,1,1]) or unit gain for simplistic purposes with our simulations
## (i) - BIAS
# lines 72 and 73:
#   Setting the master bias to be an array of zeros, with the same size as the dummy image of the simulation generated
#   earlier. The same is done for the readout noise mask, but is an array of ones - can multiply by any integer, and
#   here we have set the array to be all values of 3
## (ii) - DAKRS
# lines 84 and 85:
#   Similar process as above, but setting the master dark to an array of zeros, of the same size as the dummy image
#   created earlier. We then delete the value assigned to the dummy image
## (iii) - WHITES
#  line 89:
#   Here we create the bias/dark subtracted flat-field or master white frame, and provide its errors also. It reads
#   into the 'process_whites' function which is located in 'process_scripts.py'. This routine processes all whites from
#   a given list of files. It corrects the orientation of the image and crops the overscan regions, as well as subtracts
#   both the master bias frame [in ADU], and the master dark frame [in e-] from every image, before combining them
#   to create a master white/flat frame. This function does not do cosmic-ray removal or background subtraction.
#   NOTE: the input image has units of ADU, but the output image has units of electrons.

## STEP (3) - ORDER TRACING
# line 96:
#   Roughly finds the orders, and traces them (for a flat field Echelle spectrum). It reads the 'find_stripes' function
#   located in 'order_tracing.py'. Starting in the central column, the algorithm identifies peaks and traces each
#   stripe to the edge of the detector by following the brightest pixels along each order. It then fits a polynomial
#   to each stripe. To improve algorithm stability, the image is first smoothed with a Gaussian filter. It not only
#   eliminates noise, but also ensures that the cross-section profile of the flat becomes peaked in the middle, which
#   helps to identify the center of each stripe. Choose gauss_filter accordingly. To avoid false positives, only peaks
#   above a certain (relative) intensity threshold are used.
# lines 100-102:
#   Assigns a physical diffraction order to order-fit polynomials and bad-region masks (set to be a dummy function
#   for the time being). Reads into the two defined functions located in 'order_tracing.py'. Essentially generates a
#   dictionary with the orders and corresponding order-fit polynomials, and a second dictionary with the same orders
#   but a boolean response for bad-region masks. Then saves the polynomial dictionary to a .npy file.
# line 105:
#   Extracts the stripes of user-defined width from the science image centred on the polynomial fits defined in step
#   (1). It reads the 'extract_stripes' function defined in 'order_tracing.py'. Extracts the stripes from the original
#   2D spectrum (set as the master white) to a sparse array, containing only relevant pixels. This function marks all
#   relevant pixels for extraction. Using the provided dictionary P_id it iterates over all stripes in the image
#   and saves a sparse matrix for each stripe.
# lines 108-110:
#   Generates three dictionaries containing orders and their corresponding pixel numbers (dispersion direction),
#   extracted flux, and uncertainty of this extracted flux (incl. photon and read-out noise). Reads the
#   'extract_spectrum_from_indices' function from 'extraction.py'. This routine is simply a wrapper code for the
#   different extraction methods. We are using the quick extraction method which provides a quick-look reduction of
#   an echelle spectrum, by simply adding up the flux in a pixel column perpendicular to the dispersion direction.

## STEP (4) - PROCESS SCIENCE IMAGES
# line 153-155:
#   Inputs the 'stellar_list' of simulated (or real) data, and dictionary of polynomials (along with other parameters)
#   into the function 'process_science_images' located in 'process_scripts.py'. Key parameters in these lines include:
#   slit_height=7, scalable=True, ext_method='quick', fibs='single', timit=True.
#   Currently does the bias/dark subtraction, flat-fielding (removing pixel-to-pixel sensitivity variations), extraction
#   of stripes, and extracts the 1-dim spectra. Does not do cosmic ray removal, background extraction/estimation,
#   get relative intensities of different fibres, wavelength solution, and barycentric solution.