#%%


from diag_tcv.dbsAnalysis.correlationAnalysis.correlationAnalysis import CorrelationAnalysis, plot_correlation_slopes

import diag_tcv.dbsAnalysis.correlationAnalysis.correlationFunctions as cF
import numpy as np 
import matplotlib.pyplot as plt
import scipy
from dataAnalysis.utils.plot_utils import plot_1d, my_text, my_legend

shotobj=CorrelationAnalysis(81069, numDemod=True)
shotobj.get_normalized_data_isweep(5)

#%% Test coherence analysis with noise
dt = shotobj.processedData['sweep5']['dt']
zref = shotobj.processedData['sweep5']['z_list_ref'][3,:]
plot_1d(zref)
tcorr_no_noise, corr_no_noise, fcsd_no_noise, spectral_coh_no_noise, tcorr_scipy_no_noise, corr_scipy_no_noise= cF.full_coherence_analysis(zref, zref, dt=dt, mode='amp', plot=False)

noise_amp_real=0
noise_amp_imag=1
noise = noise_amp_real*np.random.normal(-1e-2, 1e-2, len(zref)) + noise_amp_imag*1j*np.random.normal(-1e-2, 1e-2, len(zref))

zref_noise = zref+noise
# tcorr_spec2, corr2, fcsd2, spectral_coh2, tcorr_scipy2, corr_scipy2= cF.full_coherence_analysis(zref, zref_noise, dt=dt, mode='amp', plot=False)
 # prepare signals
zref, zref_noise = cF.normalize_signals(zref, zref_noise)
zref, zref_noise = cF.choose_mode(zref, zref_noise, mode='amp')

# fourier space
fref, psd_ref, fhop, psd_hop, fcsd, csd = cF.compute_psd_csd(zref, zref_noise, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True)

# from fourier to real space
corr_from_csd, tcorr_spec = cF.compute_correlation_function(csd, dt, nperseg=1024)
corr = cF.shift_correlation_function(corr_from_csd)

import math as mp 
ind_zero_freq = len(fref)//2
psd_ref[ind_zero_freq] = mp.nan# np.mean(psd_ref[ind_zero_freq-1:ind_zero_freq+2])
psd_hop[ind_zero_freq] = mp.nan# np.mean(psd_hop[ind_zero_freq-1:ind_zero_freq+2])
csd[ind_zero_freq] = mp.nan# np.mean(csd[ind_zero_freq-1:ind_zero_freq+2])

spectral_coh = cF.compute_spectral_coherence(psd_ref, psd_hop, csd)

fig, ax = plot_1d([], [], grid=True)
ax.plot(fref/1e6, psd_ref, color='blue', label='zref')
ax.plot(fhop/1e6, psd_hop, color='red', label='zref+noise')
ax.plot(fcsd/1e6, np.abs(csd), color='green', label='csd')
ax.set_xlim(-0.5, 0.5)
ax.set_xlabel('Frequency [Hz]')




# for tests (correlation directly in real space using scipy.signal.correlate)
tcorr_scipy, corr_scipy = cF.scipy_correlation_function(zref_noise,zref, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True)
    

fig, ax = plot_1d([], [], grid=True)
ax.plot(tcorr_no_noise/1e-6, corr_no_noise, color='blue', label='zref')
ax.plot(tcorr_spec/1e-6, corr, color='red', label='zref+noise')
ax.set_xlabel(r'Time [$\mu$ s]')
my_legend(ax)
ax.set_xlim(-5, 5)
ax.set_ylabel('Correlation')

fig, ax = plot_1d([], [], grid=True)
ax.plot(fcsd_no_noise/1e6, spectral_coh_no_noise, color='blue', label='zref')
ax.plot(fcsd/1e6, spectral_coh, color='red', label='zref+noise')
ax.set_xlabel('Frequency [MHz]')
my_legend(ax)
ax.set_xlim(-2, 2)
ax.set_ylabel('Coherence')

from DBS.processing.fit_utils import perform_fits, gaussian, lorentzian
noise=0
xdata=fcsd
ydata=abs(spectral_coh)
fit_results = perform_fits(xdata, ydata, noise, dt, p0=None, verbose=False)


a_gauss = fit_results['gaussian']['params'][0]
x0 = fit_results['gaussian']['params'][1]
sigma = fit_results['gaussian']['params'][2]
ygauss = gaussian(xdata, a_gauss, x0, sigma, sigma_is_FWHM=True)

a_lorentz = fit_results['lorentzian']['params'][0]
x0 = fit_results['lorentzian']['params'][1]
FWHM = fit_results['lorentzian']['params'][2]
ylorentz = lorentzian(xdata, a_lorentz, x0, FWHM)

fig, ax = plot_1d([], [], grid=True)
ax.plot(xdata, spectral_coh, color='red')
ax.plot(xdata,ygauss, color='blue')
ax.plot(xdata, ylorentz, color='black')
 
#%% Test with real data
shot = 82607
isweep = 3

shotobj=CorrelationAnalysis(82607, numDemod=True)
shotobj.get_normalized_data_isweep(isweep)

#%%
dt = shotobj.processedData['sweep'+str(isweep)]['dt']
# zref = shotobj.processedData['sweep'+str(isweep)]['z_list_ref'][39,:]
# zhop = shotobj.processedData['sweep'+str(isweep)]['z_list_hop'][39,:]
zref = shotobj.processedData['sweep'+str(isweep)]['z_list_ref'][44,:]
zhop = shotobj.processedData['sweep'+str(isweep)]['z_list_hop'][58,:]

tcorr_old, corr_old, fcsd_old, spectral_coh_old, tcorr_scipy_old, corr_scipy_old= cF.full_coherence_analysis(zref, zhop, dt=dt, mode='amp', plot=False)

# tcorr_spec2, corr2, fcsd2, spectral_coh2, tcorr_scipy2, corr_scipy2= cF.full_coherence_analysis(zref, zref_noise, dt=dt, mode='amp', plot=False)
 # prepare signals
zref, zhop = cF.normalize_signals(zref, zhop)
zref, zhop = cF.choose_mode(zref, zhop, mode='amp')

# fourier space
fref, psd_ref, fhop, psd_hop, fcsd, csd = cF.compute_psd_csd(zref, zhop, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True)

# from fourier to real space
corr_from_csd, tcorr_spec = cF.compute_correlation_function(csd, dt, nperseg=1024)

corr = cF.shift_correlation_function(corr_from_csd)

import math as mp 
ind_zero_freq = len(fref)//2
psd_ref[ind_zero_freq] = mp.nan# np.mean(psd_ref[ind_zero_freq-1:ind_zero_freq+2])
psd_hop[ind_zero_freq] = mp.nan# np.mean(psd_hop[ind_zero_freq-1:ind_zero_freq+2])
csd[ind_zero_freq] = mp.nan# np.mean(csd[ind_zero_freq-1:ind_zero_freq+2])

spectral_coh = cF.compute_spectral_coherence(psd_ref, psd_hop, csd)

fig, ax = plot_1d([], [], grid=True)
ax.plot(fref/1e6, psd_ref, color='blue', label='zref')
ax.plot(fhop/1e6, psd_hop, color='red', label='zref+noise')
ax.plot(fcsd/1e6, np.abs(csd), color='green', label='csd')
ax.set_xlim(-0.5, 0.5)
ax.set_xlabel('Frequency [Hz]')



# for tests (correlation directly in real space using scipy.signal.correlate)
tcorr_scipy, corr_scipy = cF.scipy_correlation_function(zhop,zref, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True)
    

fig, ax = plot_1d([], [], grid=True)
ax.plot(tcorr_spec/1e-6, corr, color='red', label='')
ax.set_xlabel(r'Time [$\mu$ s]')
my_legend(ax)
ax.set_xlim(-5, 5)
ax.set_ylabel('Correlation')
my_text(ax, 0.2, 0.9, 'old max = {:.2f}'.format(np.max(abs(corr))), color='black', fontsize=12)
ax.set_ylim(-0.1, 1.1)

fig, ax = plot_1d([], [], grid=True)
ax.plot(fcsd/1e6, spectral_coh_old, color='red', label='')
# ax.plot(fcsd/1e6, spectral_coh, color='red', label='')
my_text(ax, 0.2, 0.9, 'old max = {:.2f}'.format(np.max(abs(spectral_coh_old))), color='black', fontsize=12)
ax.set_xlabel('Frequency [MHz]')
my_legend(ax)
ax.set_xlim(-2, 2)
ax.set_ylabel('Coherence')
ax.set_ylim(-0.1, 1.1)

from DBS.processing.fit_utils import perform_fits, gaussian, lorentzian, gauss_lorentz_fit_wrapper
noise=0
xdata=fcsd/1e6
ydata=abs(spectral_coh)
# ydata[xdata<-2]=0
# ydata[xdata>2]=0

bounds = ([0, -2, 0], [1, 2 , 5]) #bounds for the fit parameters: amplitude, center, FWHM
p0 = [1, 0, 2] #initial guess for the fit parameters: amplitude, center, FWHM

fitfunc, popt, pcov = custom_lorentz_fit_wrapper(xdata, ydata,p0=p0,bounds=bounds, verbose=False)

param_errors = np.sqrt(np.diag(pcov))
a_lorentz_err, x0_err, FWHM_err = param_errors

amp_err = max(a_lorentz_err, np.std(spectral_coh_old))

a_lorentz = popt[0]
x0 = popt[1]
FWHM = popt[2]
ylorentz = lorentzian(xdata, a_lorentz, x0, FWHM)

fig, ax = plot_1d([], [], grid=True)
ax.plot(xdata, ydata, color='red', label='spectral coherence')
ax.plot(xdata, ylorentz, color='black', label='lorentzian fit')
ax.set_ylabel('Coherence')
my_legend(ax, loc='upper right')
ax.set_xlim(-2, 2)
ax.set_xlabel('Frequency [MHz]')
ax.set_ylim(-0.1, 1.1)
my_text(ax, 0.2, 0.9, 'old max = {:.2f}'.format(np.max(abs(spectral_coh_old))), color='black', fontsize=12)
my_text(ax, 0.2, 0.75, 'new max = {:.2f} +/- {:.2f}'.format(a_lorentz, amp_err), color='black', fontsize=12)

#%%



# %%
#%% Method with both gaussian and lorentzian fits
dt = shotobj.processedData['sweep'+str(isweep)]['dt']
zref = shotobj.processedData['sweep'+str(isweep)]['z_list_ref'][54,:]
zhop = shotobj.processedData['sweep'+str(isweep)]['z_list_hop'][54,:]

tcorr_old, corr_old, fcsd_old, spectral_coh_old, tcorr_scipy_old, corr_scipy_old= cF.full_coherence_analysis(zref, zhop, dt=dt, mode='amp', plot=False)

# tcorr_spec2, corr2, fcsd2, spectral_coh2, tcorr_scipy2, corr_scipy2= cF.full_coherence_analysis(zref, zref_noise, dt=dt, mode='amp', plot=False)
 # prepare signals
zref, zhop = cF.normalize_signals(zref, zhop)
zref, zhop = cF.choose_mode(zref, zhop, mode='amp')

# fourier space
fref, psd_ref, fhop, psd_hop, fcsd, csd = cF.compute_psd_csd(zref, zhop, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True)

# from fourier to real space
corr_from_csd, tcorr_spec = cF.compute_correlation_function(csd, dt, nperseg=1024)
corr = cF.shift_correlation_function(corr_from_csd)

import math as mp 
ind_zero_freq = len(fref)//2
psd_ref[ind_zero_freq] = mp.nan# np.mean(psd_ref[ind_zero_freq-1:ind_zero_freq+2])
psd_hop[ind_zero_freq] = mp.nan# np.mean(psd_hop[ind_zero_freq-1:ind_zero_freq+2])
csd[ind_zero_freq] = mp.nan# np.mean(csd[ind_zero_freq-1:ind_zero_freq+2])

spectral_coh = cF.compute_spectral_coherence(psd_ref, psd_hop, csd)

fig, ax = plot_1d([], [], grid=True)
ax.plot(fref/1e6, psd_ref, color='blue', label='zref')
ax.plot(fhop/1e6, psd_hop, color='red', label='zref+noise')
ax.plot(fcsd/1e6, np.abs(csd), color='green', label='csd')
ax.set_xlim(-0.5, 0.5)
ax.set_xlabel('Frequency [Hz]')



# for tests (correlation directly in real space using scipy.signal.correlate)
tcorr_scipy, corr_scipy = cF.scipy_correlation_function(zhop,zref, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True)
    

fig, ax = plot_1d([], [], grid=True)
ax.plot(tcorr_spec/1e-6, corr, color='red', label='')
ax.set_xlabel(r'Time [$\mu$ s]')
my_legend(ax)
ax.set_xlim(-5, 5)
ax.set_ylabel('Correlation')
my_text(ax, 0.2, 0.9, 'old max = {:.2f}'.format(np.max(abs(corr))), color='black', fontsize=12)
ax.set_ylim(-0.1, 1.1)

fig, ax = plot_1d([], [], grid=True)
ax.plot(fcsd/1e6, spectral_coh_old, color='red', label='')
# ax.plot(fcsd/1e6, spectral_coh, color='red', label='')
my_text(ax, 0.2, 0.9, 'old max = {:.2f}'.format(np.max(abs(spectral_coh_old))), color='black', fontsize=12)
ax.set_xlabel('Frequency [MHz]')
my_legend(ax)
ax.set_xlim(-2, 2)
ax.set_ylabel('Coherence')
ax.set_ylim(-0.1, 1.1)

from DBS.processing.fit_utils import perform_fits, gaussian, lorentzian
noise=0
xdata=fcsd
ydata=abs(spectral_coh)
fit_results = perform_fits(xdata, ydata, noise, dt, p0=None, verbose=False)


a_gauss = fit_results['gaussian']['params'][0]
x0 = 0#fit_results['gaussian']['params'][1]
sigma = fit_results['gaussian']['params'][2]
ygauss = gaussian(xdata, a_gauss, x0, sigma, sigma_is_FWHM=True)

a_lorentz = fit_results['lorentzian']['params'][0]
x0 = 0#fit_results['lorentzian']['params'][1]
FWHM = fit_results['lorentzian']['params'][2]
ylorentz = lorentzian(xdata, a_lorentz, x0, FWHM)

fig, ax = plot_1d([], [], grid=True)
ax.plot(xdata/1e6, abs(spectral_coh), color='red', label='spectral coherence')
ax.plot(xdata/1e6,ygauss, color='blue', label='gaussian fit')
ax.plot(xdata/1e6, ylorentz, color='black', label='lorentzian fit')
ax.set_ylabel('Coherence')
my_legend(ax, loc='upper right')
ax.set_xlim(-2, 2)
ax.set_xlabel('Frequency [MHz]')
ax.set_ylim(-0.1, 1.1)
my_text(ax, 0.2, 0.9, 'old max = {:.2f}'.format(np.max(abs(spectral_coh_old))), color='black', fontsize=12)
my_text(ax, 0.2, 0.75, 'new max = {:.2f} +/- {:.2f}'.format(np.mean([a_lorentz]), np.std([a_gauss, a_lorentz])), color='black', fontsize=12)













#%% Doing my own fit
bounds = ([0, -np.inf, 0], [np.inf, np.inf, np.inf])


from scipy.optimize import curve_fit

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




#%% Test for a full radial correlation function
shot = 81069
isweep = 5
plat = 1

if plat == 0:
    indmin=0
    indmax=20
elif plat == 1:
    indmin=20
    indmax=40
elif plat == 2:
    indmin=40
    indmax=60

shotobj=CorrelationAnalysis(shot, numDemod=True)
shotobj.get_normalized_data_isweep(isweep)

dt = shotobj.processedData['sweep'+str(isweep)]['dt']
maxfitcoh = []
maxfitcoherr = []
maxfitcoh_old = []

for i in range(indmin,indmax):
    

    zref = shotobj.processedData['sweep'+str(isweep)]['z_list_ref'][i,:]
    zhop = shotobj.processedData['sweep'+str(isweep)]['z_list_hop'][i,:]


    zref, zhop = cF.normalize_signals(zref, zhop)
    zref, zhop = cF.choose_mode(zref, zhop, mode='amp')

    # fourier space
    fref, psd_ref, fhop, psd_hop, fcsd, csd = cF.compute_psd_csd(zref, zhop, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True)
    _, psd_ref_old, _, psd_hop_old, _, csd_old = cF.compute_psd_csd(zref, zhop, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True)

    # from fourier to real space
    corr_from_csd, tcorr_spec = cF.compute_correlation_function(csd, dt, nperseg=1024)

    corr = cF.shift_correlation_function(corr_from_csd)

    import math as mp 
    ind_zero_freq = len(fref)//2
    psd_ref[ind_zero_freq] = mp.nan# np.mean(psd_ref[ind_zero_freq-1:ind_zero_freq+2])
    psd_hop[ind_zero_freq] = mp.nan# np.mean(psd_hop[ind_zero_freq-1:ind_zero_freq+2])
    csd[ind_zero_freq] = mp.nan# np.mean(csd[ind_zero_freq-1:ind_zero_freq+2])

    spectral_coh = cF.compute_spectral_coherence(psd_ref, psd_hop, csd)
    spectral_coh_old = cF.compute_spectral_coherence(psd_ref_old, psd_hop_old, csd_old)

    # for tests (correlation directly in real space using scipy.signal.correlate)
    tcorr_scipy, corr_scipy = cF.scipy_correlation_function(zhop,zref, dt, nperseg=1024, noverlap=512, window=None, remove_mean=True)
        
    xdata=fcsd/1e6
    ydata=abs(spectral_coh)

    ydata[xdata<-1]=0
    ydata[xdata>1]=0

    bounds = ([0, -2, 0], [max(ydata), 2 , 5]) #bounds for the fit parameters: amplitude, center, FWHM
    p0 = [max(ydata), 0, 2] #initial guess for the fit parameters: amplitude, center, FWHM

    popt, pcov = custom_lorentz_fit_wrapper(xdata, ydata,p0=p0,bounds=bounds, verbose=False)
    # popt, pcov = custom_lorentz_fit_wrapper(xdata, ydata,p0=p0,bounds=None, verbose=False)

    param_errors = np.sqrt(np.diag(pcov))
    a_lorentz_err, x0_err, FWHM_err = param_errors

    amp_err = max(a_lorentz_err, np.std(spectral_coh_old))

    a_lorentz = popt[0]
    x0 = popt[1]
    FWHM = popt[2]
    ylorentz = lorentzian(xdata, a_lorentz, x0, FWHM)

    maxfitcoh.append(a_lorentz)
    maxfitcoherr.append(amp_err)
    maxfitcoh_old.append(np.max(abs(spectral_coh_old)))
    
    fig, ax = plot_1d([], [], grid=True)
    ax.plot(xdata, ydata, color='red', label='spectral coherence')
    ax.plot(xdata, ylorentz, color='black', label='lorentzian fit')
    ax.set_ylabel('Coherence')
    my_legend(ax, loc='upper right')
    ax.set_xlim(-2, 2)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylim(-0.1, 1.1)
    my_text(ax, 0.2, 0.9, 'old max = {:.2f}'.format(np.max(abs(spectral_coh_old))), color='black', fontsize=12)
    my_text(ax, 0.2, 0.75, 'new max = {:.2f} +/- {:.2f}'.format(a_lorentz, amp_err), color='black', fontsize=12)
    

#%%
delta = shotobj.get_delta(isweep, retdata=True)
fig, ax = plot_1d([], [], grid=True)
ax.errorbar(delta[indmin:indmax], maxfitcoh, yerr=maxfitcoherr, fmt='^', color='red')
ax.errorbar(delta[indmin:indmax], maxfitcoh_old, fmt='o', color='teal')

ax.set_yscale('log')
ax.set_xlabel('delta cm')
ax.set_ylabel('max coherence')
ax.set_title('#{} sweep{}'.format(shot, isweep))
# ax.set_ylim(0.1,1.2)
# %%
