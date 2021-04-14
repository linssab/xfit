import numpy as np
import os, sys
from .Constants import *

def reta(A,B,size):
    return np.array([(A*i) + B for i in range(size)])

def savgol_filter(y, window_size, order, deriv=0, rate=1):
    from math import factorial

    """ This function was taken from https://gist.github.com/krvajal
    and compared to scipy.signal package, a reliable widely known package
    for signal analysis.

    References
    [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688 """

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")

    order_range = range(order+1)
    half_window = (window_size -1) // 2
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def gaus(x, E_peak, gain, Noise, Fano, *A):
    """ Model function to be fitted.

    ---------------------------------------------

    INPUT:
        x; <1D-array> - contains the calibrated energy axis
        E_peak; <1D-array> - contains the energy where peaks where found
        gain; <float> - energy gain per channel in electronvolts
        Noise; <float> - Noise factor
        Fano; <float> - Fano factor
        A; <1D-array> - Gaussian amplitude array """

    s = np.sqrt(((Noise/2.3548200450309493)**2)+3.85*Fano*E_peak) #np.sqrt works for arrays
    return gain*np.sum(
            A/(s*2.5066282746310002)*np.exp(-np.square(x[:,None]-E_peak)/(2*(s**2))),1)

def linregress(x, y, sigmay=None, full_output=False):
    """
    Linear fit to a straight line following P.R. Bevington:
    "Data Reduction and Error Analysis for the Physical Sciences"
    Parameters
    ----------
    x, y : array_like
        two sets of measurements.  Both arrays should have the same length.
    sigmay : The uncertainty on the y values
    Returns
    -------
    slope : float
        slope of the regression line
    intercept : float
        intercept of the regression line
    r_value : float
        correlation coefficient
    if full_output is true, an additional dictionary is returned with the keys
    sigma_slope: uncertainty on the slope
    sigma_intercept: uncertainty on the intercept
    stderr: float
        square root of the variance
    """

    x = np.asarray(x, dtype=np.float).flatten()
    y = np.asarray(y, dtype=np.float).flatten()
    N = y.size
    if sigmay is None:
        sigmay = np.ones((N,), dtype=y.dtype)
    else:
        sigmay = np.asarray(sigmay, dtype=np.float).flatten()
    w = 1.0 / (sigmay * sigmay + (sigmay == 0))

    n = S = w.sum()
    Sx = (w * x).sum()
    Sy = (w * y).sum()
    Sxx = (w * x * x).sum()
    Sxy = ((w * x * y)).sum()
    Syy = ((w * y * y)).sum()
    # SSxx is identical to delta in Bevington book
    delta = SSxx = (S * Sxx - Sx * Sx)

    tmpValue = Sxx * Sy - Sx * Sxy
    intercept = tmpValue / delta
    SSxy = (S * Sxy - Sx * Sy)
    slope = SSxy / delta
    sigma_slope = np.sqrt(S /delta)
    sigma_intercept = np.sqrt(Sxx / delta)

    SSyy = (n * Syy - Sy * Sy)
    r_value = SSxy / np.sqrt(SSxx * SSyy)
    if r_value > 1.0:
        r_value = 1.0
    if r_value < -1.0:
        r_value = -1.0

    if not full_output:
        return slope, intercept, r_value

    ddict = {}
    # calculate the variance
    if N < 3:
        variance = 0.0
    else:
        variance = ((y - intercept - slope * x) ** 2).sum() / (N - 2)
    ddict["variance"] = variance
    ddict["stderr"] = np.sqrt(variance)
    ddict["slope"] = slope
    ddict["intercept"] = intercept
    ddict["r_value"] = r_value
    ddict["sigma_intercept"] = np.sqrt(Sxx / SSxx)
    ddict["sigma_slope"] = np.sqrt(S / SSxx)
    return slope, intercept, r_value, ddict

def tophat(y,v,w):
    """ Applies a tophat filter to a 1D-array following the recipe
    given reference below:
    [1] Van Grieken, R.E.; Markowicz, A.A. Handbook of X-ray spectrometry, 
    Marcel Dekker, Inc., 2002.
    
    --------------------------------------------------------------------

    INPUT:
        y; 1D-array
        v; side window width (int)
        w; central window width (int)
    OUTPUT:
        new_y; filtered pectrum (1D-array)
        delta; variance (1D-array)
        peaks; identified peaks indexes (1D-list)
    """
    
    def hk_(k,v,w):

        if -v-(w/2) <= k < -w/2: return -1/(2*v)
        elif -w/2 <= k <= w/2: return 1/w
        elif w/2 < k <= v+(w/2): return -1/(2*v)
        else: return 0

    new_y = np.zeros([y.shape[0]])
    delta = np.zeros([y.shape[0]])
    peaks = []
    for i in range(y.shape[0]-(int(v+w/2))):
        y_star = 0
        y_delta = 0
        for k in range(int(-v-w/2),int(v+w/2)):
            hk = hk_(k,v,w)
            y_star += hk*y[i+k]
            y_delta += abs((hk**2)*y[i+k])
        new_y[i] = y_star
        delta[i] = y_delta**.5
    for j in range(1,y.shape[0]-1):
        if new_y[j-1] <= new_y[j] and new_y[j] > new_y[j+1]:
            peaks.append(j)
    return new_y, delta, peaks

def findpeaks(spectrum,w=9,r=1):
    v = int(w/2)+1
    y,var,indexes = tophat(spectrum,v,w)
    r*=(np.max(y)/(np.max(np.sqrt(var)*r)))/130
    selected_peaks = []
    for i in indexes:
        if y[i]>r*(var[i]**0.5) and i>120:
            selected_peaks.append(i)
    return selected_peaks

def FN_fit_gaus(spec,spec_bg,e_axis,gain):
    """ Finds a series of relevant peaks between 2100 and 30000 eV and fits them
    with a gaussian fuction to obtain the global Fano and Noise factors.

    ----------------------------------------------------------------------------
    INPUT:
        spec; <1D-array> - global spectrum to be used in the peak sampling and fit
        spec_bg; <1D-array> - continuum of the above spectrum
        e_axis; <1D-array> - calibrated energy axis in KeV
        gain; <float> - calibration gain in KeV
    OUTPUT:
        Fano; <float> - Fano factor of the input spectrum
        Noise; <float> - Noise factor of the input spectrum """

    from scipy.optimize import curve_fit
    import copy

    ####################################

    Fano = 0.114
    Noise = 80
    gain = gain*1000
    energyaxis = e_axis*1000
    zero = energyaxis[0]
    y_array, y_cont = spec, spec_bg

    ####################################

    #########################
    # perform deconvolution #
    #########################

    w = PEAK_TOLERANCE
    v = int(w/2)+1
    r = 2 
    peaks = []
    y_conv, y_delta, pre_peaks = tophat(y_array,v,w)
    for i in pre_peaks:
        if y_conv[i-1]<=y_conv[i]>=y_conv[i+1]: 
            if y_conv[i]>(r*y_delta[i]) and \
                    y_conv[i]>y_cont[i]*CONTINUUM_SUPPRESSION: 
                peaks.append(i)

    #########################

    ######################################
    # filters peaks and exclude outliers #
    ######################################

    filtered_peaks, E_peaks, y_peaks = [],[],[]
    for i in peaks:
        if 2100<=energyaxis[i]<=30000:
            E_peaks.append(energyaxis[i])
            filtered_peaks.append(i)
            y_peaks.append(y_array[i])

    peaks = np.asarray(filtered_peaks)
    E_peaks = np.asarray(E_peaks) 
    y_peaks = np.asarray(y_peaks)

    ######################################

    guess = y_peaks*np.sqrt(
            ((Noise/2.3548)**2)+(3.85*Noise*E_peaks))*np.sqrt(2*np.pi)/gain
    guess = np.insert(guess,0,[Noise,Fano],axis=0)
    uncertainty = np.sqrt(y_array).clip(1)
    
    try:
        popt_gaus, pcov_gaus = curve_fit(lambda x,Noise,Fano,*A: gaus(
                energyaxis,E_peaks,gain,Noise,Fano,
                    *A) + y_cont,
                energyaxis,
                y_array,
                p0=guess,
                sigma=uncertainty,
                maxfev=FIT_CYCLES)
    except: 
        raise ValueError("Failed to fit fano and noise. Continuing with default")
        return Fano, Noise

    #print(popt_gaus.shape)
    #print("NOISE",popt_gaus[0])
    #print("FANO",popt_gaus[1])

    #y_gaus = np.asarray(
    #        [gaus(energyaxis,E_peaks[i],gain,*popt_gaus[[0,1,2+i]]) \
    #                for i in range(E_peaks.size)])+ y_cont
    #for i in y_gaus:
    #    plt.semilogy(energyaxis,i,label="Fit")
    #plt.semilogy(energyaxis,y_array,label="Data")
    #plt.semilogy(energyaxis,y_cont,label="Continuum")
    #plt.legend()
    #plt.show()
    print(popt_gaus[1], popt_gaus[0])

    return popt_gaus[1], popt_gaus[0]
