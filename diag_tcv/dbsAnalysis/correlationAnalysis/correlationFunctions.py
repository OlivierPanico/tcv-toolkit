# %% Full coherence analysis function 
import numpy as np 
import matplotlib.pyplot as plt
from dataAnalysis.utils.utils import get_closest_ind, normalize_array_1d, my_linearRegression
from dataAnalysis.utils.array_splitting import custom_split_1d
from dataAnalysis.utils.plot_utils import plot_1d, my_text, my_legend
from dataAnalysis.spectral.spectralAnalysis import custom_csd, custom_time_coherence
import scipy
import scipy.signal as signal
from scipy.optimize import curve_fit
import math as mp

'''
Input: zref, zhop, dt
Output: corr, spectral_coh
'''


def normalize_complex_1d(a):
    a_r = normalize_array_1d(a.real)
    a_i = normalize_array_1d(a.imag)
    return a_r + 1j*a_i

# Step 0: normalization of the signals
def normalize_signals(zref, zhop):
    zref_norm = normalize_complex_1d(zref)
    zhop_norm = normalize_complex_1d(zhop)
    return zref_norm, zhop_norm

# Step 0bis: choose mode
def choose_mode(zref_norm, zhop_norm, mode='full'):
    '''
    mode='full' => compute coherence on the full complex signal
    mode='amp' => compute coherence on the amplitude signal
    mode='real' => compute coherence on the real signal
    mode='imag' => compute coherence on the imag signal
    mode='phase' => compute coherence on the phase signal
    '''
    if mode == 'amp':
        zref_norm = np.sqrt(zref_norm.real**2 + zref_norm.imag**2)
        zhop_norm = np.sqrt(zhop_norm.real**2 + zhop_norm.imag**2)
    elif mode=='real':
        zref_norm = zref_norm.real
        zhop_norm = zhop_norm.real
    elif mode=='imag':
        zref_norm = zref_norm.imag
        zhop_norm = zhop_norm.imag
    elif mode=='phase':
        zref_norm = np.angle(zref_norm)
        zhop_norm = np.angle(zhop_norm)
    else:
        pass
    return zref_norm, zhop_norm

# Step 1: compute psd and csd
def compute_psd_csd(zref_norm, zhop_norm, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True):
    fref, psd_ref = custom_csd(zref_norm, zref_norm, dt=dt, nperseg=nperseg, noverlap=noverlap, window=window,norm=True, remove_mean=remove_mean)
    fhop, psd_hop = custom_csd(zhop_norm, zhop_norm,dt=dt, nperseg=nperseg, noverlap=noverlap, window=window, norm=True, remove_mean=remove_mean)
    fcsd, csd = custom_csd(zref_norm, zhop_norm, dt=dt, nperseg=nperseg, noverlap=noverlap, window=window, norm=True, remove_mean=remove_mean)
    return fref, psd_ref, fhop, psd_hop, fcsd, csd

# Step 2: normalization of csd by psd to obtain spectral coherence
def compute_spectral_coherence(psd_ref, psd_hop, csd):
    spectral_coh = abs(csd)**2/(psd_ref*psd_hop)
    return spectral_coh

# Step 3: ifftshift + inverse Fourier transform to obtain the correlation function
def compute_correlation_function(csd, dt, nperseg = 1024, mode='full'):
    csd_ifftshift = np.fft.ifftshift(csd)
    corr_from_csd = np.fft.ifft(csd_ifftshift)/dt       
    tcorr_spec = scipy.signal.correlation_lags(nperseg, nperseg, mode='same')*dt
    return corr_from_csd, tcorr_spec

# Step 4: fftshift on the correlation function to replace peak at the center
def shift_correlation_function(corr_from_csd):
    corr = np.fft.fftshift(corr_from_csd)
    return corr


# second possibility: use the scipy.signal.correlate function
def scipy_correlation_function(zref_norm, zhop_norm, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True):
    
    zref_split = custom_split_1d(zref_norm, nperseg, noverlap=noverlap, zero_padding=False)
    zhop_split = custom_split_1d(zhop_norm, nperseg, noverlap=noverlap, zero_padding=False)
    nbseg = len(zref_split[:,0])
    corr = np.zeros((nperseg), dtype=zref_norm.dtype)
    
    for i in range(nbseg):
        if remove_mean:
            zref_seg = zref_split[i,:]-np.mean(zref_split[i,:])
            zhop_seg = zhop_split[i,:]-np.mean(zhop_split[i,:])
        else:
            zref_seg = zref_split[i,:]
            zhop_seg = zhop_split[i,:]
        
        corr_seg=signal.correlate(zref_seg, zhop_seg, mode='same') / np.sqrt(signal.correlate(zref_seg,zref_seg, mode='same')[int(nperseg/2)] * signal.correlate(zhop_seg,zhop_seg, mode='same')[int(nperseg/2)])
        corr+=corr_seg/nbseg
    
    delay_arr = signal.correlation_lags(nperseg,nperseg, mode='same')*dt
               
    return delay_arr, corr


# Full coherence analysis function
def full_coherence_analysis(zref, zhop, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True, plot=False,ax=None, mode='full', verbose=False):

    # prepare signals
    zref_norm, zhop_norm = normalize_signals(zref, zhop)
    zref_norm, zhop_norm = choose_mode(zref_norm, zhop_norm, mode=mode)

    # fourier space
    fref, psd_ref, fhop, psd_hop, fcsd, csd = compute_psd_csd(zref_norm, zhop_norm, dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean)
    spectral_coh = compute_spectral_coherence(psd_ref, psd_hop, csd)
    
    # from fourier to real space
    corr_from_csd, tcorr_spec = compute_correlation_function(csd, dt, nperseg=nperseg)
    corr = shift_correlation_function(corr_from_csd)
    
    # for tests (correlation directly in real space using scipy.signal.correlate)
    tcorr_scipy, corr_scipy = scipy_correlation_function(zhop_norm,zref_norm, dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean)
    
    # Estimate the maximum of the correlation functions
    a_lorentz, amp_err, max_spectral_coh_raw = estimate_max_spectral_coh(fcsd, spectral_coh, plot=plot)
    max_corr_delay, max_corr = get_max_corr_delay(tcorr_spec, corr)
    max_corr_scipy_delay, max_corr_scipy = get_max_corr_delay(tcorr_scipy, corr_scipy)
    
    if plot:
        # Plotting the results
        if ax is None:
            fig, ax = plot_1d([], [], grid=True)
            # fig2, ax2 = plot_1d([], [], grid=True)
        # ax.plot(tcorr_spec*1e6,(corr), color='darkorchid', label='correlation function', markersize=2)
        ax.plot(tcorr_spec*1e6,np.sqrt(corr.real**2+corr.imag**2), label='amplitude', marker='')
        ax.set_xlabel(r'delay $[\mu_s]$')
        ax.set_ylabel('correlation')
        # plt.legend()
        # plt.title('Spectral correlation function')
        ax.set_xlim(-5, 5)
        ax.set_ylim(-0.1, 1)
        
        # ax2.plot(fcsd/1e6, spectral_coh, label='spectral coherence', marker='')
        # ax2.set_xlim(-2,2)
        # ax2.set_ylim(-0.1, 1.1)
        # ax2.set_xlabel(r'frequency $[MHz]$')
        # ax2.set_ylabel('coherence')

    
    dictFullCohAnalysis = dict()
    dictFullCohAnalysis['tcorr_spec'] = tcorr_spec
    dictFullCohAnalysis['corr'] = corr
    dictFullCohAnalysis['fcsd'] = fcsd
    dictFullCohAnalysis['spectral_coh'] = spectral_coh
    dictFullCohAnalysis['tcorr_scipy'] = tcorr_scipy
    dictFullCohAnalysis['corr_scipy'] = corr_scipy
    dictFullCohAnalysis['max_fit_spectral_coh'] = a_lorentz
    dictFullCohAnalysis['err_max_fit_spectral_coh'] = amp_err
    dictFullCohAnalysis['max_raw_spectral_coh'] = max_spectral_coh_raw
    dictFullCohAnalysis['max_corr_delay'] = max_corr_delay
    dictFullCohAnalysis['max_corr'] = max_corr
    dictFullCohAnalysis['max_corr_scipy_delay'] = max_corr_scipy_delay
    dictFullCohAnalysis['max_corr_scipy'] = max_corr_scipy
    
    return dictFullCohAnalysis


# step 5: get maximum of correlation & associated delay 
def get_max_corr_delay(tcorr, corr):
    max_corr = np.max(abs(corr))
    max_corr_ind = get_closest_ind(abs(corr), max_corr)
    max_corr_delay = tcorr[max_corr_ind]
    return max_corr_delay, max_corr

# Define a Lorentzian function
def lorentzian(x, a, x0, FWHM):
    return a / ((x-x0)**2 / (FWHM/2)**2 + 1)

def custom_lorentz_fit_wrapper(xdata, ydata, p0=None, verbose=False, **kwargs):
    
    # remove any NaNs from the data:
    mask = np.isnan(ydata) | np.isnan(xdata)
    xdata = xdata[~mask]
    ydata = ydata[~mask]

    # initial guess for fit parameters
    bounds = kwargs.pop('bounds', None)
    
    if p0 is None:
        # Initial guess for the Gaussian fit parameters
        p0 = [np.max(ydata), np.mean(xdata), np.std(xdata)]
        
    if not bounds:
        # we just need to make sure that amplitude and FWHM are positive (otherwise, the subsequent Taylor fit might fail due to poor initial guess)
        bounds = ([0, -np.inf, 0], [np.inf, np.inf, np.inf])
    
    popt, pcov = curve_fit(lorentzian, xdata, ydata, p0=p0, bounds=bounds, **kwargs)
    
    return popt, pcov

def estimate_max_spectral_coh(fcsd, spectral_coh, plot=False):

    xdata=fcsd/1e6
    ydata=abs(spectral_coh)



    #We remove the zero frequency
    ydata_for_fit = ydata.copy()
    ydata_for_fit[len(xdata)//2]= mp.nan
    
    # cleaning large frequencies
    ydata[xdata<-1]=0
    ydata[xdata>1]=0
    ydata_for_fit[xdata<-1]=0
    ydata_for_fit[xdata>1]=0
    
    # bounds and initial guess => that is the tricky part
    bounds = ([0, -2, 0], [max(ydata), 2 , 5]) #bounds for the fit parameters: amplitude, center, FWHM
    p0 = [max(ydata), 0, 2] #initial guess for the fit parameters: amplitude, center, FWHM


    # fit
    popt, pcov = custom_lorentz_fit_wrapper(xdata, ydata_for_fit,p0=p0,bounds=bounds, verbose=False)

    # error is estimated as the max of the standard deviation and the amplitude of the spectral coherence
    param_errors = np.sqrt(np.diag(pcov))
    a_lorentz_err, x0_err, FWHM_err = param_errors
    amp_err = max(a_lorentz_err, np.std(ydata))

    a_lorentz = popt[0]
    x0 = popt[1]
    FWHM = popt[2]
    max_spectral_coh_raw = np.max(abs(ydata))
    
    if abs(max_spectral_coh_raw-a_lorentz)/max_spectral_coh_raw > 0.3:
        print('Warning: the maximum of the spectral coherence is not well estimated by the Lorentzian fit: raw max = {:.2f}, est. max = {:.2f} +/- {:.2f}'.format(max_spectral_coh_raw, a_lorentz, amp_err))
    
    
    if plot:
        ylorentz = lorentzian(xdata, a_lorentz, x0, FWHM)
        
        fig, ax = plot_1d([], [], grid=True)
        ax.plot(xdata, ydata, color='red', label='spectral coherence')
        ax.plot(xdata, ylorentz, color='black', label='lorentzian fit')
        ax.set_ylabel('Coherence')
        my_legend(ax, loc='upper right')
        ax.set_xlim(-2, 2)
        ax.set_xlabel('Frequency [MHz]')
        ax.set_ylim(-0.1, 1.1)
        my_text(ax, 0.2, 0.9, 'raw max = {:.2f}'.format(np.max(abs(ydata))), color='black', fontsize=12)
        my_text(ax, 0.2, 0.75, 'est. max = {:.2f} +/- {:.2f}'.format(a_lorentz, amp_err), color='black', fontsize=12)
        
    return a_lorentz, amp_err, max_spectral_coh_raw
    
    
