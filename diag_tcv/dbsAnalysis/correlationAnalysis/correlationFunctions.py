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
Goal: compute the spatial correlation function 
Input: zref, zhop, dt
Output: corr as a function of time, max corr 


There are few different methods to compute the correlation function 

Method 1 (Panico phd thesis)
1) normalize signals zref and zhop
2) take amplitude of zref and zhop 
3) Perform psd and csd (take care with normalizations)
4) ifft of the csd (take care with normalizations) => correlation vs time
5) fit the correlation function to estimate the maximum of correlation
6) Do this for all frequencies in a plateau will give the spatial correlation function 

Remarks: 
- taking the amplitude of zref and zhop leads to a double slope in the spatial correlation function that is interpreted as avalanche signal
- taking amplitude signal removes time delays in the time correlation function (maximum is centered on zero time delay)



Method 2
1) normalize signals zref and zhop
2) take the full complex signal of zref and zhop
3) Perform psd and csd (take care with normalizations)
4) ifft of the amplitude of the csd => correlation vs time
5) take maximum of correlation 
6) Do this for all frequencies in a plateau will give the spatial correlation function 

Remarks: 
- no double slopes 
- performs better in case of noise: maximum of correlation is lowered but the spatial behaviour is still coherent 
- in some case we can infer a time delay 

Alternative to method 2: 
1), 2), 3) are similar 
4) we perform a taylor fit of the csd similar to when we estimate the doppler shift for the velocity
5) ifft of the fit 
6) maximum of correlation 
7) do this for all frequencies 

Remarks: 
- Increase slightly the overall correlation as compared to method 2
- But does not reach very low values 
- still no double slopes
- cannot fit the csd that originates from amplitude signal because the signal is not adapted to the fit functions

Method 3: 
Perform directly the correlation in real space using scipy.signal.correlate

remarks: 
- gives similar results as method 2
- can be used as a cross-check

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
        zref_norm = np.exp(1J*np.angle(zref_norm))
        zhop_norm = np.exp(1J*np.angle(zhop_norm))
    else:
        pass
    return zref_norm, zhop_norm


### ======================================================== ###
### method 1: 
### ======================================================== ###

# Step 1: compute psd and csd
def compute_psd_csd(zref_norm, zhop_norm, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True):
    fref, psd_ref = custom_csd(zref_norm, zref_norm, dt=dt, nperseg=nperseg, noverlap=noverlap, window=window,norm=True, remove_mean=remove_mean)
    fhop, psd_hop = custom_csd(zhop_norm, zhop_norm,dt=dt, nperseg=nperseg, noverlap=noverlap, window=window, norm=True, remove_mean=remove_mean)
    fcsd, csd = custom_csd(zref_norm, zhop_norm, dt=dt, nperseg=nperseg, noverlap=noverlap, window=window, norm=True, remove_mean=remove_mean)
    return fref, psd_ref, fhop, psd_hop, fcsd, csd

# Step 2: normalization of csd by psd to obtain spectral coherence
def compute_spectral_coherence(psd_ref, psd_hop, csd):
    spectral_coh = abs(csd)**2/(psd_ref*psd_hop)
    # spectral_coh = abs(csd)/np.sqrt(psd_ref*psd_hop)
    return (spectral_coh)

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
    a_lorentz, amp_err, max_spectral_coh_raw = estimate_max_spectral_coh(fcsd, spectral_coh, mode=mode, plot=plot)
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
    dictFullCohAnalysis['csd'] = csd
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



### ======================================================== ###
### method 2: 
### ======================================================== ###
'''
We use DBS fit methods that's why there is a specobj structure
'''
from DBS.processing.sigprocessing import get_noise_level, get_center_of_gravity
from DBS.processing.fit_utils import gaussian, lorentzian

from DBS.processing.sigprocessing import nantrapz, plot_semi_log10, get_noise_level, get_frequency_range, get_spikes_in_PSD, get_even_odd_spectrum, get_center_of_gravity
import scipy.signal as signal
from matlabtools import Struct

ODD_EVEN_THRESHOLD = 0.1

def preprocessing_csd(zref, zhop, dt, noise_freq_range=[0.9*4.9e6, 4.9e6], nperseg=1024, noverlap=512, fmin=50e3, fmax=10e6, csdmode='amp'):
    """Compute PSD and related spectral quantites of a given complex signal z. Returns a `specobj` structure containing the results."""

    fref, Pref = custom_csd(zref, zhop, dt=dt, nperseg=nperseg, noverlap=noverlap, window=None, norm=False, remove_mean=False)
    fhop, Phop = custom_csd(zref, zhop, dt=dt, nperseg=nperseg, noverlap=noverlap, window=None, norm=False, remove_mean=False)
    fcsd, Pcsd = custom_csd(zref, zhop, dt=dt, nperseg=nperseg, noverlap=noverlap, window=None, norm=True, remove_mean=False)
    Pcsd_full = Pcsd.copy()
    
    if csdmode=='amp':
        Pcsd = np.sqrt(Pcsd.real**2 + Pcsd.imag**2) # take the magnitude of the csd
    elif csdmode=='real':  
        Pcsd = Pcsd.real # take only the real part of the csd
    elif csdmode=='imag':
        Pcsd = Pcsd.imag # take only the imaginary part of the csd
    elif csdmode=='phase':
        # Pcsd = np.angle(Pcsd) # take only the phase of the csd
        Pcsd = np.exp(1J*np.angle(Pcsd)) # take only the phase of the csd

    
    P_noise = get_noise_level(fcsd, Pcsd, *noise_freq_range).squeeze()

    # remove frequencies below/above some threshold:
    include_mask = get_frequency_range(fcsd, fmin, fmax) # convert to int to avoid problems when exporting/importing 
    
    # remove spikes in the csd:
    is_spike_mask = get_spikes_in_PSD(fcsd, Pcsd, threshold=3)
    include_mask = include_mask & ~is_spike_mask

    # symmetric/asymmetric components of the spectrum
    P_even, P_odd, fsym = get_even_odd_spectrum(fcsd,Pcsd, assume_centered=False)
    P_odd_positive = ( P_odd + np.abs(P_odd) ) / 2
    # _include_mask = get_frequency_range(fsym, fmin, fmax) # somehow this leads to wrong CoG estimation, so not use it for now

    # center of gravity of the asymmetric spectrum:
    f_cog = get_center_of_gravity(fsym, P_odd_positive)
    # f_cog_PH = get_center_of_gravity(diag_data.fsym, diag_data.P_dop)

    # compute the power of the odd and even components:
    power_odd =  nantrapz(P_odd_positive, x=fsym)
    power_even =  nantrapz(P_even, x=fsym)

    # print('noise level [V^2/Hz]: {:.3e} (PH: {:.3e})'.format(P_noise, diag_data.Bdop))
    # print('center of gravity [MHz]: {:.3f} (PH: {:.3f})'.format(f_cog/1e6, f_cog_PH/1e6))
    # print('odd/even power [V^2] {:.3e} / {:.3e} = {:.3e} (PH: {:.3e})'.format(power_odd, power_even, power_odd/power_even, diag_data.Roe))
    
    specobjcsd = Struct()
    specobjcsd.validated = 0
    specobjcsd.dt = dt
    specobjcsd.f = fcsd
    specobjcsd.Pcsd = Pcsd #chosen
    specobjcsd.Pcsd_full = Pcsd_full #full complex pcsd
    specobjcsd.include_mask = include_mask
    specobjcsd.P_noise = P_noise
    specobjcsd.fsym = fsym
    specobjcsd.P_even = P_even
    specobjcsd.P_odd = P_odd
    specobjcsd.P_odd_positive = P_odd_positive
    specobjcsd.f_cog = f_cog
    specobjcsd.power_odd = power_odd
    specobjcsd.power_even = power_even
    specobjcsd.odd_even_threshold = ODD_EVEN_THRESHOLD

    return specobjcsd

def perform_specobj_fits_csd(specobjcsd, include_mask=None, p0=None, reinitialize=False, verbose=False):
    
    from DBS.processing import perform_fits
    s = specobjcsd  # shortcut
    
    if include_mask is None:
        include_mask = s.include_mask # get_frequency_range(s.f, fmin=50e3, fmax=8e6)
    
    include_mask = include_mask.astype(bool)
    
    if not hasattr(s, 'fsym'):
        reinitialize = True
    
    xscale = 1e5
    yscale = np.nanmax(s.Pcsd[include_mask]) 
    xdata = s.f[include_mask] / xscale
    ydata = s.Pcsd[include_mask] / yscale
    noise = s.P_noise / yscale
    dt = specobjcsd.dt

    s.include_mask = include_mask
    s.xscale = xscale
    s.yscale = yscale
        
    try:
        odd_even_ratio = s.power_odd / s.power_even
        # print(f'odd/even power [V^2] {s.power_odd:.3e} / {s.power_even:.3e} = {odd_even_ratio:.3e}')
        
        fit_results_full = Struct()
        
        # Gaussian fit of the odd part (only if the odd/even ratio is above a certain threshold):
        threshold = s.odd_even_threshold if hasattr(s, 'odd_even_threshold') else ODD_EVEN_THRESHOLD
        if odd_even_ratio > threshold:
        
            xdata_odd = s.fsym / xscale
            ydata_odd = s.P_odd_positive / yscale
            from DBS.processing.fit_utils import gauss_lorentz_fit_wrapper
            curve_func, curve_params = gauss_lorentz_fit_wrapper(xdata_odd, ydata_odd, curve_type='gaussian')
            res_odd = Struct({'func': curve_func, 'params': curve_params})
            fit_results_full.update({'odd': res_odd})
            
            if p0 is None:
                p0 = np.hstack([np.nanmax(ydata),  res_odd['params'][1:]])
                
                if verbose:
                    print(f'init guess odd:\n',p0)# keep the center and FWHM of the Gaussian fit of the odd part as initial guess for the Gaussian fit of the whole spectrum
                
                    print(f'results odd:\n', fit_results_full.odd.params)
        else:
            if p0 is None:
                
                cog = get_center_of_gravity(xdata, ydata)
                p0 = [np.nanmax(ydata), cog, 0.5e6 / xscale] # initial guess for the Gaussian fit parameters
                if verbose:
                    print(f'init guess raw:\n',p0)

        try:
            fit_results_full.update(perform_fits(xdata, ydata, noise, dt * xscale, p0=p0, verbose=verbose))
            
            if verbose:
                print(f'results Gauss:\n', fit_results_full.gaussian.params)
                print('results Lorentz\n', fit_results_full.lorentzian.params)
                
            # create also a stripped down version of the fits (just the fit parameters, not the functions), that will be later saved to file:
            s.fit_params = Struct()
            for curve_type, fit_result in fit_results_full.items():
                s.fit_params[curve_type] = fit_result['params']
        except RuntimeError:
            if verbose:
                print(f'Fitting failed for specobj with header:\n{s.header}')
    
    except Exception as e:
        traceback.print_exc()
        warnings.warn(f'Fitting failed for specobj with header:\n{s.header}')
    

def show_spec_csd(specobjcsd, ax=None, axis_labels=True, 
              legend=True,
              show_curves=['rawPSD', 'fittedPSD', 'oddPSD', 'fit_odd', 'fit_gaussian', 'fit_lorentzian', 'fit_taylor', 'annotations'],
              verbose=False,
              indicate_estimate=False,
              **plot_kwargs):
    """ Plots the spectral data contained in specobj.
    :param specobj: object containing the spectral data
    :param ax: matplotlib (or other) axis object to plot the data
    :param axis_labels: whether to show axis labels
    :param legend: whether to show the legend
    :param show: list of strings indicating which curves to show
    :return: plot_dict: dictionary containing the plot objects
    """
    colors = {
            'gaussian': 'blue',
            'lorentzian': 'darkgreen',
            'taylor': 'black',
            'odd': 'magenta',
            'even': 'lightgreen'
        }
    import matplotlib.pyplot as plt
    from DBS.processing.fit_utils import fitfuncs
    plot_dict = {}    
    s = specobjcsd # shortcut
    
    if ax is None:
        fig, ax = plt.subplots()
    
    include_mask = s.include_mask
    xscale = s.xscale
    yscale = s.yscale
    xdata = s.f[include_mask] / xscale
    ydata = s.Pcsd[include_mask] / yscale
    
    for show_c in show_curves:
        
        if show_c == 'rawPSD':
            # determine the bounds of the plot:
            fmin, fmax = np.min(xdata), np.max(xdata)
            ind_f = (s.f / xscale > fmin) & (s.f / xscale < fmax)
            _xraw, _yraw = np.copy(s.f) / 1e6, np.copy(s.P)
            _xraw[~ind_f] = np.nan
            _yraw[~ind_f] = np.nan
            lrawPSD = plot_semi_log10(ax, _xraw, _yraw, label='raw PSD', alpha=0.7, **plot_kwargs)
            plot_dict['lrawPSD'] = lrawPSD
            
        elif show_c == 'fittedPSD':
            lfittedPSD = plot_semi_log10(ax, xdata * xscale / 1e6, ydata * yscale, color='red', label='fitted part', **plot_kwargs)
            plot_dict['lfittedPSD'] = lfittedPSD
        
        elif show_c == 'oddPSD' and hasattr(s, 'fsym'):
            loddPSD = plot_semi_log10(ax, s.fsym / 1e6, s.P_odd_positive + s.P_noise, color=colors['odd'], label='odd part', ls='--', **plot_kwargs)
            plot_dict['loddPSD'] = loddPSD
        
        elif show_c == 'evenPSD' and hasattr(s, 'fsym'):
            levenPSD = plot_semi_log10(ax, s.fsym / 1e6, s.P_even, color=colors['even'], label='even part', ls='--', **plot_kwargs)
            plot_dict['levenPSD'] = levenPSD

        elif show_c in ['fit_odd', 'fit_gaussian', 'fit_lorentzian', 'fit_taylor']:
            curve_type = show_c.split('_')[1]
            if f'lfit_{curve_type}' in plot_dict.keys(): # skip if already plotted
                if verbose:
                    print(f'skipping {show_c} as it is already plotted')
            else:
                if hasattr(s, 'fit_params') and curve_type in s.fit_params.keys():                    
                    
                    curve_func   = fitfuncs[curve_type]
                    curve_params = s.fit_params[curve_type]
                    xfit = np.linspace(np.min(xdata), np.max(xdata), 1000)
                    
                    if curve_type=='gaussian':
                        print('Gaussian fit params:', curve_params)
                    
                    if curve_type == 'taylor':
                        yfit = curve_func(xfit, *curve_params, dt=s.dt * xscale)
                    else:
                        yfit = curve_func(xfit, *curve_params)
                        
                    x = xfit * xscale / 1e6
                    y = yfit * yscale + s.P_noise
                    l = plot_semi_log10(ax, x, y, label=f'{curve_type} fit',color=colors[curve_type], **plot_kwargs)
                    plot_dict[f'lfit_{curve_type}'] = l
                
                else:
                    if verbose:
                        print(f'No fit results found for {curve_type}, skipping')
            
            
    

    if not hasattr(s, 'fit_params'):
        if verbose:
            print('No fit results found, skipping fit plots')
    else:
        
        if 'annotations' in show_curves:
    
            fDop_fits = np.array([fitpar[1] for fitpar in s.fit_params.values()]) * xscale
            
            # add line indicating the Center of Gravity of the odd part of the spectrum:
            cog_vline = ax.axvline(s.f_cog/1e6, color=colors['odd'], ls='-.', label='odd CoG')
            
            # indicate noise level:
            noise_hline = ax.axhline(10*np.log10(s.P_noise), color='black', ls='--', alpha=0.5, lw=0.5, label='noise level')
            
            plot_dict['cog_vline'] = cog_vline
            plot_dict['noise_hline'] = noise_hline
            

    
    ax.set_ylim(bottom=10*np.log10(0.5*s.P_noise))
    
    if axis_labels:
        ax.set_xlabel('Frequency [MHz]')
        ax.set_ylabel('PSD [dB]')
    
    if legend:
        ax.legend()
        
    ax.grid(True, which='both', ls='--', alpha=0.5)
    ax.axvline(0, color='black', lw=0.5)
    
    return plot_dict


def get_correlation_from_fit(specobjcsd, fit_type='gaussian'):
    from scipy.fft import ifft, fftshift, ifftshift
    dt = specobjcsd.dt
    f = specobjcsd.f
    Pcsd = specobjcsd.Pcsd
    Pcsd_full = specobjcsd.Pcsd_full
    P_noise = specobjcsd.P_noise
    xscale = specobjcsd.xscale
    yscale = specobjcsd.yscale

    xdata = f / xscale
    ydata = Pcsd / yscale

    xfit = np.linspace(np.min(xdata), np.max(xdata), len(xdata))
    
    if fit_type=='gaussian':
        from DBS.processing.fit_utils import gaussian
        fit_func = gaussian
    elif fit_type=='taylor':
        from DBS.processing.fit_utils import taylor
        fit_func = taylor
    elif fit_type=='lorentzian':
        from DBS.processing.fit_utils import lorentzian
        fit_func = lorentzian
    else:
        raise ValueError(f'Unknown fit type: {fit_type}')
    
    fit_params = specobjcsd.fit_params[fit_type]
    
    if fit_type == 'taylor':
        data_fit = fit_func(xfit, *fit_params, dt=dt*xscale, nfft=4096)
    else:
        data_fit = fit_func(xfit, *fit_params)

    x = xfit * xscale / 1e6
    
    y = data_fit * yscale + P_noise

    #ifft_from_fit = np.fft.ifftshift(y)
    #corr_from_fit = np.fft.ifft(ifft_from_fit)/dt
    #tcorr_spec = signal.correlation_lags(len(xfit), len(xfit), mode='same')*dt
    #corr_from_fit = np.fft.fftshift(corr_from_fit)

    ifft_from_csd_full = np.fft.ifftshift(Pcsd_full)
    corr_from_csd_full = np.fft.ifft(ifft_from_csd_full)/dt
    tcorr = signal.correlation_lags(len(f), len(f), mode='same')*dt
    corr_from_csd_full = np.fft.fftshift(corr_from_csd_full)

    return tcorr, corr_from_csd_full


def correlation_wrapper(zref, zhop, dt, nperseg=1024, noverlap=512, noise_freq_range=[0.9*4.9e6, 4.9e6], fmin=50e3, fmax=10e6, fit_type='taylor', csdmode = 'full', verbose=False):
    
    specobjcsd = preprocessing_csd(zref, zhop, dt, noise_freq_range=noise_freq_range, nperseg=nperseg, noverlap=noverlap, fmin=fmin, fmax=fmax, csdmode=csdmode)
    perform_specobj_fits_csd(specobjcsd, include_mask=None, p0=None, reinitialize=False, verbose=verbose)
    tcorr, corr_from_csd_full = get_correlation_from_fit(specobjcsd, fit_type=fit_type)

    return specobjcsd, tcorr, corr_from_csd_full


### ======================================================== ###
### method 3: use the scipy.signal.correlate function
### ======================================================== ###
def scipy_correlation_function(zref_norm, zhop_norm, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True, method='auto'):
    '''
    method: 'auto', 'direct', 'fft'
    '''
    
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





### ======================================================== ###
### Utils for fitting in method 1
### ======================================================== ###

# Define a Lorentzian function
def lorentzian(x, a, x0, FWHM):
    return a / ((x-x0)**2 / (FWHM/2)**2 + 1)


def gaussian(x, a, x0, sigma, sigma_is_FWHM=True):
    if sigma_is_FWHM:
        sigma = sigma / (2 * np.sqrt(2 * np.log(2)))
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def custom_lorentz_fit_wrapper(xdata, ydata, curve_type=lorentzian, p0=None, verbose=False, **kwargs):
    
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
    
    popt, pcov = curve_fit(curve_type, xdata, ydata, p0=p0, bounds=bounds, **kwargs)
    
    return popt, pcov

def estimate_max_spectral_coh(fcsd, spectral_coh, mode='amp', plot=False, **kwargs):

    default_kwargs = {'color':'red', 'label':'spectral coherence'}
    default_kwargs.update(kwargs)
    
    xdata=fcsd/1e6
    ydata=abs(spectral_coh)



    #We remove the zero frequency
    ydata_for_fit = ydata.copy()
    ydata_for_fit[len(xdata)//2]= mp.nan
    
    # cleaning large frequencies
    if mode == 'amp':
        ydata[xdata<-1]=0
        ydata[xdata>1]=0
        ydata_for_fit[xdata<-1]=0
        ydata_for_fit[xdata>1]=0
        
    # bounds and initial guess => that is the tricky part
    if mode != 'amp':
        bounds = ([0, -5, 0], [max(ydata), 5 , 10]) #bounds for the fit parameters: amplitude, center, FWHM
    else:
        bounds = ([0, -2, 0], [max(ydata), 2 , 5]) #bounds for the fit parameters: amplitude, center, FWHM
    p0 = [max(ydata), 0, 2] #initial guess for the fit parameters: amplitude, center, FWHM


    # fit
    if mode != 'amp':
        popt, pcov = custom_lorentz_fit_wrapper(xdata, ydata_for_fit, curve_type=gaussian, p0=p0,bounds=bounds, verbose=False)
        # error is estimated as the max of the standard deviation and the amplitude of the spectral coherence
        amp = popt[0]
        x0 = popt[1]
        sigma = popt[2]
        yfit = gaussian(xdata, *popt)
        param_errors = np.sqrt(np.diag(pcov))
        a_gauss_err, x0_err, sigma_err = param_errors
        residuals = ydata_for_fit - yfit
        # amp_err = max(a_gauss_err, np.std(residuals))
        amp_err = max(a_gauss_err, np.std(residuals), abs(amp-np.max(abs(ydata))))
        
    else:
        popt, pcov = custom_lorentz_fit_wrapper(xdata, ydata_for_fit, curve_type=lorentzian, p0=p0,bounds=bounds, verbose=False)
         # error is estimated as the max of the standard deviation and the amplitude of the spectral coherence
        param_errors = np.sqrt(np.diag(pcov))
        a_lorentz_err, x0_err, FWHM_err = param_errors
        amp = popt[0]
        x0 = popt[1]
        FWHM = popt[2]
        yfit = lorentzian(xdata, *popt)
        residuals = ydata_for_fit - yfit
        
        amp_err = max(a_lorentz_err, np.std(residuals), abs(amp-np.max(abs(ydata))))
        
    max_spectral_coh_raw = np.max(abs(ydata))
    
    if abs(max_spectral_coh_raw-amp)/max_spectral_coh_raw > 0.3:
        print('Warning: the maximum of the spectral coherence is not well estimated by the fit: raw max = {:.2f}, est. max = {:.2f} +/- {:.2f}'.format(max_spectral_coh_raw, amp, amp_err))
    
    
    if plot:
        
        
        fig, ax = plot_1d([], [], grid=True)
        ax.plot(xdata, ydata, **default_kwargs)
        if mode != 'amp':
            ax.plot(xdata, yfit, color='black', label='gaussian fit')
        else:
            ax.plot(xdata, yfit, color='black', label='lorentzian fit')
        ax.set_ylabel('Coherence')
        my_legend(ax, loc='upper right')
        ax.set_xlim(-2, 2)
        ax.set_xlabel('Frequency [MHz]')
        ax.set_ylim(-0.1, 1.1)
        my_text(ax, 0.2, 0.9, 'raw max = {:.2f}'.format(np.max(abs(ydata))), color='black', fontsize=12)
        my_text(ax, 0.2, 0.75, 'est. max = {:.2f} +/- {:.2f}'.format(amp, amp_err), color='black', fontsize=12)
    
        return amp, amp_err, max_spectral_coh_raw
    
    return amp, amp_err, max_spectral_coh_raw
    
    
