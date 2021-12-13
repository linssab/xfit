import numpy as np
import matplotlib.pyplot as plt
from .SimpleMath import savgol_filter

def strip(an_array,cycles,width):
    """
    Strips the peaks contained in the input array following an
    iteractive peak clipping method.

    Reference: 
    [1] van Grieken, R.E.; Markowicz, A.A. Handbook of X-ray spectrometry,
    Marcel Dekker, Inc., 2002.
    [2] Clayton CA, Hines JW, Elkins PD. Anal Chem 59:2506, 1987.

    ----------------------------------------------------------------------

    INPUT:
        an_array; np.array
        cycles; int
        width; int
    OUTPUT:
        an_array; np.array
    """

    ##################################################################
    # W IS THE WIDTH OF THE FILTER. THE WINDOW WILL BE (1*W)+1       #
    # W VALUE MUST BE LARGER THAN 2 AND ODD, SINCE 3 IS THE MINIMUM  #
    # SATISFACTORY POLYNOMIAL DEGREE TO SMOOTHEN THE DATA            #
    ##################################################################
    
    size = an_array.shape[0]
    for k in range(cycles):
        if k >= cycles-8:
            width = int(width/1.4142135623730950) #square root of 2
        for l in range(0, size):
            if l-width < 0: low = 0
            else: low = l-width
            if l+width >= size: high = size-1
            else: high = l+width
            m = (an_array[low] + an_array[high])*0.5
            if m < 1: m = 1
            if an_array[l] > m: an_array[l] = m
    return an_array

def peakstrip(an_array, cycles, width, savgol_window, order):
    """
    Calculates the continuum/background contribution of the input spectrum following
    the SNIPBG method.

    Reference:
    [1] van Grieken, R.E.; Markowicz, A.A. Handbook of X-ray spectrometry,
    Marcel Dekker, Inc., 2002.

    --------------------------------------------------------------------------------

    INPUT:
        an_array; np.array
        cycles; int
        width; int
        savgol_window: int, > order
        order: int, odd
    OUTPUT:
        snip_bg; np.array
    """

    TEST = False

    ############################
    # initialize snip_bg array #
    ############################

    snip_bg = np.zeros(an_array.shape[0])
    
    ########################
    # square root the data #
    ########################

    sqr_data = an_array**0.5

    ###############################################
    # apply savgol filter to the square root data #
    ###############################################

    smooth_sqr = savgol_filter(sqr_data,savgol_window,order)
    
    for i in range(smooth_sqr.shape[0]): 
        if smooth_sqr[i] < 0: smooth_sqr[i] = 0
    
    ################## NUMBA PRESENTS SERIOUS ISSUES DEPENDING ON VERSION ################
    # THIS SNIPPET WORKS WITH NUMBA 0.42.1 (win10 only) and 0.45.1 (win7 AND win10)      #
    # numba functions do not return an output, the output must be passed as an argument. #
    # in this case, smooth_sqr is CHANGED.                                               #
    ######################################################################################
     
    strip(smooth_sqr,cycles,width)

    ##################################################################################
    # change snip_bg to the transformed smooth_sqr after applying the strip function #
    ##################################################################################

    snip_bg = smooth_sqr[:] 

    ##################
    # transform back #
    ##################

    snip_bg **= 2
       
    if TEST == True:
        plt.semilogy(snip_bg,label="background")
        plt.semilogy(an_array,label="spectrum")
        plt.legend()
        plt.show()

    return snip_bg

