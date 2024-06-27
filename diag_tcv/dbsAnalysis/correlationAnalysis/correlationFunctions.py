# %% Full coherence analysis function 
import numpy as np 
import matplotlib.pyplot as plt
from dataAnalysis.utils.utils import get_closest_ind, normalize_array_1d
from dataAnalysis.utils.array_splitting import custom_split_1d
from dataAnalysis.utils.plot_utils import plot_1d
from dataAnalysis.spectral.spectralAnalysis import custom_csd, custom_time_coherence
import scipy
import scipy.signal as signal

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

# Step 1: compute psd and csd
def compute_psd_csd(zref_norm, zhop_norm, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True):
    fref, psd_ref = custom_csd(zref_norm, zref_norm, dt=dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean)
    fhop, psd_hop = custom_csd(zhop_norm, zhop_norm,dt=dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean)
    fcsd, csd = custom_csd(zref_norm, zhop_norm, dt=dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean)
    return fref, psd_ref, fhop, psd_hop, fcsd, csd

# Step 2: normalization of csd by psd to obtain spectral coherence
def compute_spectral_coherence(psd_ref, psd_hop, csd):
    spectral_coh = abs(csd)**2/(psd_ref*psd_hop)
    return spectral_coh

# Step 3: ifftshift + inverse Fourier transform to obtain the correlation function
def compute_correlation_function(csd, dt):
    csd_ifftshift = np.fft.ifftshift(csd)
    corr_from_csd = np.fft.ifft(csd_ifftshift)/dt/2
    return corr_from_csd

# Step 4: fftshift on the correlation function to replace peak at the center
def shift_correlation_function(corr_from_csd):
    corr = np.fft.fftshift(corr_from_csd)
    return corr


# second possibility: use the scipy.signal.correlate function
def scipy_correlation_function(zref_norm, zhop_norm, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True):
    '''
    usually does not work as well as the custom function
    '''
    
    zref_split = custom_split_1d(zref_norm, nperseg, noverlap=noverlap, zero_padding=False)
    zhop_split = custom_split_1d(zhop_norm, nperseg, noverlap=noverlap, zero_padding=False)
    nbseg = len(zref_split[:,0])
    corr = np.zeros((nperseg), dtype=zref_norm.dtype)
    
    for i in range(nbseg):
        zref_seg = zref_split[i,:]-np.mean(zref_split[i,:])
        zhop_seg = zhop_split[i,:]-np.mean(zhop_split[i,:])
        
        
        corr_seg=signal.correlate(zref_seg, zhop_seg, mode='same') / np.sqrt(signal.correlate(zref_seg,zref_seg, mode='same')[int(nperseg/2)] * signal.correlate(zhop_seg,zhop_seg, mode='same')[int(nperseg/2)])
        corr+=corr_seg/nbseg
    
    delay_arr = signal.correlation_lags(nperseg,nperseg, mode='same')*dt
               
    return delay_arr, corr


# Full coherence analysis function
def full_coherence_analysis(zref, zhop, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True, plot=False, mode='full', verbose=False):
    '''
    mode='full' => compute coherence on the full complex signal
    mode='amp' => compute coherence on the amplitude signal
    mode='real' => compute coherence on the real signal
    mode='imag' => compute coherence on the imag signal
    mode='phase' => compute coherence on the phase signal
    '''
    
    zref_norm, zhop_norm = normalize_signals(zref, zhop)
    
        
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

    fref, psd_ref, fhop, psd_hop, fcsd, csd = compute_psd_csd(zref_norm, zhop_norm, dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean)
    spectral_coh = compute_spectral_coherence(psd_ref, psd_hop, csd)
    corr_from_csd = compute_correlation_function(csd, dt)
    corr = shift_correlation_function(corr_from_csd)
    tcorr_spec = scipy.signal.correlation_lags(nperseg, nperseg, mode='same')*dt
    
    # for tests
    tcorr_scipy, corr = scipy_correlation_function(zref_norm, zhop_norm, dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean)
    
    if plot:
        # Plotting the results
        
        fig, ax = plot_1d([], [], grid=True)
        ax.plot(tcorr_spec*1e6,(corr), color='darkorchid', label='correlation function', markersize=2)
        ax.plot(tcorr_spec*1e6,np.sqrt(corr.real**2+corr.imag**2), color='indigo', label='amplitude', markersize=2)
        ax.set_xlabel(r'delay $[\mu_s]$')
        ax.set_ylabel('correlation')
        plt.legend()
        plt.title('Spectral correlation function')
        # plt.xlim(-20, 20)
    
    if verbose:
        print('Full coherence analysis done')
    
    return tcorr_spec, corr, fcsd, spectral_coh


