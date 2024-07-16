#%% Mixed files technique for noise estimation


#1. Import the necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy

from diag_tcv.dbsAnalysis.correlationAnalysis.correlationAnalysis import CorrelationAnalysis
import diag_tcv.dbsAnalysis.correlationAnalysis.correlationFunctions as cF

from dataAnalysis.utils.plot_utils import plot_1d, my_text, my_legend

def correlation_mixed_files(zref_norm, zhop_norm, dt, plot=False):

    nperseg=1024
    noverlap=512
    window=None
    remove_mean=True
    # plot=True
    mode='amp'
    verbose=False

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

    fref, psd_ref, fhop, psd_hop, fcsd, csd = cF.compute_psd_csd(zref_norm, zhop_norm, dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean)
    spectral_coh = cF.compute_spectral_coherence(psd_ref, psd_hop, csd)
    corr_from_csd = cF.compute_correlation_function(csd, dt)
    corr = cF.shift_correlation_function(corr_from_csd)
    tcorr_spec = scipy.signal.correlation_lags(nperseg, nperseg, mode='same')*dt

    tcorr_scipy, corr_scipy = cF.scipy_correlation_function(zref_norm, zhop_norm, dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean)
        
    if plot:
        #4. Plot the correlation function and spectral coherence
        fig, ax = plot_1d([], [], grid=True)
        ax.plot(fref, abs(psd_ref), color='blue', label='psd ref', marker='')
        ax.plot(fhop, abs(psd_hop), color='red', label='psd hop', marker='')
        ax.plot(fcsd, abs(csd), color='green', label='csd', marker='')
        ax.set_yscale('log')
        my_legend(ax)

        fig, ax = plot_1d([], [], grid=True)
        ax.plot(fref, spectral_coh, color='darkorchid', label='spectral coherence', marker='')
        ax.set_xlabel(r'frequency [Hz]')
        ax.set_ylabel('coherence')
        my_legend(ax)
        plt.title('Spectral coherence')

        fig, ax = plot_1d([], [], grid=True)
        ax.plot(tcorr_spec*1e6,(corr), color='darkorchid', label='correlation function', markersize=2)
        ax.plot(tcorr_spec*1e6,np.sqrt(corr.real**2+corr.imag**2), color='indigo', label='amplitude', markersize=2)
        ax.plot(tcorr_scipy*1e6, abs(corr_scipy), color='red', label='scipy', markersize=2)
        ax.set_xlabel(r'delay $[\mu_s]$')
        ax.set_ylabel('correlation')
        my_legend(ax)
        plt.title('Spectral correlation function')


    return np.max(corr_scipy), np.max(corr), np.max(spectral_coh)




#%% Example of the mixed files technique for noise estimation
#2. Load zref and zhop for different discharges & isweep
shot1 = 80949
isweep1 = 6
ifreq1 = [1]
a=CorrelationAnalysis(shot1, numDemod=True)
z_list_ref1, z_list_hop1, t_reduced_list_ref1, t_reduced_list_hop1 = a.get_normalized_data_isweep(isweep1, ifreq_list=ifreq1, dtsart=400.e-6, dtend=100.e-6, ret=True)

shot2 = 80940
isweep2 = 9
ifreq2 = [8]
b=CorrelationAnalysis(shot2, numDemod=True)
z_list_ref2, z_list_hop2, t_reduced_list_ref2, t_reduced_list_hop2 = b.get_normalized_data_isweep(isweep2, ifreq_list=ifreq2, dtsart=400.e-6, dtend=100.e-6, ret=True)

#%% Plot the two signals
fig, ax = plot_1d([],[],grid=True)
# ax.plot(t_reduced_list_ref1[0], z_list_ref1[0], label='Shot {}, isweep {}, ifreq {}'.format(shot1, isweep1, ifreq1[0]), marker='')
ax.plot(z_list_ref1[0], label='Shot {}, isweep {}, ifreq {}'.format(shot1, isweep1, ifreq1[0]), marker='')
ax2=ax.twinx()
ax2.plot(z_list_hop2[0], label='Shot {}, isweep {}, ifreq {}'.format(shot2, isweep2, ifreq2[0]), marker='')
ax.set_title('#{} isweep {} ifreq {} vs #{} isweep {} ifreq {}'.format(shot1, isweep1, ifreq1[0], shot2, isweep2, ifreq2[0]))
ax.set_xlabel('samples')
ax.set_ylabel('zref1')
ax2.set_ylabel('zhop2')

zref_norm = z_list_ref1[0]
zhop_norm = z_list_hop2[0]
dt = t_reduced_list_ref1[0][1] - t_reduced_list_ref1[0][0]
correlation_mixed_files(zref_norm, zhop_norm, dt, plot=True)

# %% Perform mixed files for 100 random pairs of shots
shots = [80940, 80949, 81069, 81084]
isweep = [4,5,6,7,8,9]
ifreq = np.linspace(0,39, 40)

import random
maxlist=[]
for i in range(100):
    
    #shot1
    shot1 = random.choice(shots)
    isweep1 = random.choice(isweep)
    ifreq1 = [random.choice(ifreq)]
    
    #shot2
    shot2 = random.choice(shots)
    isweep2 = random.choice(isweep)
    ifreq2 = [random.choice(ifreq)]
    
    print(' === chosen pair of shots === ')    
    print('shot1: {}, isweep1: {}, ifreq1: {}'.format(shot1, isweep1, ifreq1[0]))
    print('shot2: {}, isweep2: {}, ifreq2: {}'.format(shot2, isweep2, ifreq2[0]))
    
    
    a=CorrelationAnalysis(shot1, numDemod=True)
    z_list_ref1, z_list_hop1, t_reduced_list_ref1, t_reduced_list_hop1 = a.get_normalized_data_isweep(isweep1, ifreq_list=ifreq1, dtsart=400.e-6, dtend=100.e-6, ret=True)
    b=CorrelationAnalysis(shot2, numDemod=True)
    z_list_ref2, z_list_hop2, t_reduced_list_ref2, t_reduced_list_hop2 = b.get_normalized_data_isweep(isweep2, ifreq_list=ifreq2, dtsart=400.e-6, dtend=100.e-6, ret=True)
    dt=t_reduced_list_ref1[0][1] - t_reduced_list_ref1[0][0]
    
    max1, max2, max3 = correlation_mixed_files(z_list_ref1[0], z_list_hop2[0], dt, plot=False)
    print('max1: {}, max2: {}, max3: {}'.format(max1, max2, max3))
    print(' ============================ ')

    maxlist.append(np.max([max1,max2,max3]))


# %%
fig, ax = plot_1d(maxlist, grid=True, label='max correlation', linestyle='')
ax.set_xlabel('random pairs')
ax.set_ylabel('max correlation')
plt.title('max correlation for 100 random pairs')
ax.axhline(np.mean(maxlist), color='red', label='mean')
ax.axhline(np.mean(maxlist)+np.std(maxlist), color='orange', label='mean+std')
ax.axhline(np.mean(maxlist)-np.std(maxlist), color='orange', label='mean-std')
my_legend(ax)
my_text(ax,0.7,0.5, r'mean: {:.2f} $\pm$ {:.2f}'.format(np.real(np.mean(maxlist)), np.real(np.std(maxlist))), color='red')


#%% mixed files at specific frequencies 
shots = [80940, 80949, 81069, 81084]
isweep = [4,5,6,7,8,9]
ifreq = np.linspace(2,37, 36)

import random
maxlist=[]
for i in range(100):
    
    #shot1
    shot1 = random.choice(shots)
    isweep1 = random.choice(isweep)
    ifreq1 = [random.choice(ifreq)]
    
    #shot2
    shot2 = random.choice(shots)
    isweep2 = random.choice(isweep)
    ifreq_list_constrained = np.linspace(ifreq1[0]-2, ifreq1[0]+2, 5)
    ifreq2 = [random.choice(ifreq_list_constrained)]
    
    print(' === chosen pair of shots === ')    
    print('shot1: {}, isweep1: {}, ifreq1: {}'.format(shot1, isweep1, ifreq1[0]))
    print('shot2: {}, isweep2: {}, ifreq2: {}'.format(shot2, isweep2, ifreq2[0]))
    
    
    a=CorrelationAnalysis(shot1, numDemod=True)
    z_list_ref1, z_list_hop1, t_reduced_list_ref1, t_reduced_list_hop1 = a.get_normalized_data_isweep(isweep1, ifreq_list=ifreq1, dtsart=400.e-6, dtend=100.e-6, ret=True)
    b=CorrelationAnalysis(shot2, numDemod=True)
    z_list_ref2, z_list_hop2, t_reduced_list_ref2, t_reduced_list_hop2 = b.get_normalized_data_isweep(isweep2, ifreq_list=ifreq2, dtsart=400.e-6, dtend=100.e-6, ret=True)
    dt=t_reduced_list_ref1[0][1] - t_reduced_list_ref1[0][0]
    
    max1, max2, max3 = correlation_mixed_files(z_list_ref1[0], z_list_hop2[0], dt, plot=False)
    print('max1: {}, max2: {}, max3: {}'.format(max1, max2, max3))
    print(' ============================ ')

    maxlist.append(np.max([max1,max2,max3]))
    
#%% 
fig, ax = plot_1d(maxlist, grid=True, label='max correlation', linestyle='')
ax.set_xlabel('random pairs')
ax.set_ylabel('max correlation')
plt.title('max correlation for 100 random pairs')
ax.axhline(np.mean(maxlist), color='red', label='mean')
ax.axhline(np.mean(maxlist)+np.std(maxlist), color='orange', label='mean+std')
ax.axhline(np.mean(maxlist)-np.std(maxlist), color='orange', label='mean-std')
my_legend(ax)
my_text(ax,0.7,0.5, r'mean: {:.2f} $\pm$ {:.2f}'.format(np.real(np.mean(maxlist)), np.real(np.std(maxlist))), color='red')
# %%
