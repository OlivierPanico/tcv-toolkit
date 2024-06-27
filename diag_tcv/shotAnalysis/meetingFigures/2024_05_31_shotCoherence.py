#%%
### ======= ###
### IMPORTS ###
### ======= ###
from DBS.beamtracing.DBSbeam import _DBSbeam
from diag_tcv.dbsAnalysis.correlationAnalysis.correlationAnalysis import CorrelationAnalysis


#%%
### ================ ###
### LOADING DBS BEAM ###
### ================ ###
# shotnb = 80940
# isweep_list = [4,5,6,7,8,9]

# for i, isweep in enumerate(isweep_list):
#     _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=3, verbose=True, plot=True, load_if_existing=True)
#     _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=4, verbose=True, plot=True, load_if_existing=True)
    
# shotnb = 81065
# isweep_list = [4,5,6,7,8,9]

# for i, isweep in enumerate(isweep_list):
#     _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=3, verbose=True, plot=True, load_if_existing=True)
#     _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=4, verbose=True, plot=True, load_if_existing=True)


# shotnb = 81084
# isweep_list = [4,5,6,7,8,9]
# for i, isweep in enumerate(isweep_list):
#     _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=3, verbose=True, plot=True, load_if_existing=True)
#     _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=4, verbose=True, plot=True, load_if_existing=True)


# shotnb = 81069
# isweep_list = [4,5,6,7,8,9]
# for i, isweep in enumerate(isweep_list):
#     _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=3, verbose=True, plot=True, load_if_existing=True)
#     _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=4, verbose=True, plot=True, load_if_existing=True)


# shotnb = 81087
# isweep_list = [4,5,6,7,8,9]
# for i, isweep in enumerate(isweep_list):
#     _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=3, verbose=True, plot=True, load_if_existing=True)
#     _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=4, verbose=True, plot=True, load_if_existing=True)

#%%
shotnb = 80338
isweep_list = [4,5,6,7,8,9]
for i, isweep in enumerate(isweep_list):
    _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=3, verbose=True, plot=True, load_if_existing=True)
    _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=4, verbose=True, plot=True, load_if_existing=True)


### =========================================== ###
### INFLUENCE OF HEATING ON CORRELATION LENGTHS ###
### =========================================== ###
# Hydrogen: 80940 ; 80949
# Deuterium: 81065, 81084, 81069, 81087
#%% 80940
list_isweep=[4,5,6,7,8,9]
a=CorrelationAnalysis(80940, numDemod=True, verbose=False, plot=True)
a.plot_coherence(list_isweep, plot_pearson=True, plot_spectral=False)
a.plot_coherence(list_isweep, plot_pearson=False, plot_spectral=True)
a.plot_heating_with_sweeps(list_isweep)

#%% 80949
list_isweep=[4,5,6,7,8,9]
a=CorrelationAnalysis(80949, numDemod=True, verbose=False, plot=False)
a.plot_coherence(list_isweep, plot_pearson=True, plot_spectral=False)
a.plot_coherence(list_isweep, plot_pearson=False, plot_spectral=True)
a.plot_heating_with_sweeps(list_isweep)

#%% 81065
list_isweep=[4,5,6,7,8,9]
a=CorrelationAnalysis(81065, numDemod=True, verbose=False, plot=False)
a.plot_coherence(list_isweep, plot_pearson=True, plot_spectral=False)
a.plot_coherence(list_isweep, plot_pearson=False, plot_spectral=True)
a.plot_heating_with_sweeps(list_isweep)

#%% 81084
list_isweep=[4,5,6,7,8,9]
a=CorrelationAnalysis(81084, numDemod=True, verbose=False, plot=False)
a.plot_coherence(list_isweep, plot_pearson=True, plot_spectral=False)
a.plot_coherence(list_isweep, plot_pearson=False, plot_spectral=True)
a.plot_heating_with_sweeps(list_isweep)

#%% 81069
ylog=True
list_isweep=[4,5,6,7,8,9]
a=CorrelationAnalysis(81069, numDemod=True, verbose=False, plot=False)
a.plot_coherence_delta(list_isweep, plot_pearson=True, plot_spectral=False, ylog=ylog)
a.plot_coherence_delta(list_isweep, plot_pearson=False, plot_spectral=True, ylog=ylog)
a.plot_heating_with_sweeps(list_isweep)

#%% 81087 (ECRH)
ylog=True
mode='amp'
plot_spectral=False 
plot_pearson=True
add_pt_zero=False

list_isweep=[4,5,6,7,8,9]

a=CorrelationAnalysis(81087, numDemod=True, verbose=False, plot=False)
a.plot_coherence_delta(list_isweep, plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero, plot_rho=True)
a.plot_coherence_delta(list_isweep, plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero, plot_rho=True)
a.plot_heating_with_sweeps(list_isweep)



#%% 80338
ylog=True
mode='amp'
plot_spectral=False 
plot_pearson=True
add_pt_zero=False

list_isweep=[2, 3, 4, 5]

a=CorrelationAnalysis(80338, numDemod=True, verbose=False, plot=False)
a.plot_coherence(list_isweep, plot_pearson=True, plot_spectral=False)
a.plot_coherence(list_isweep, plot_pearson=False, plot_spectral=True)
a.plot_heating_with_sweeps(list_isweep)

### =========================================== ###
### INFLUENCE OF ISOTOPE ON CORRELATION LENGTHS ###
### =========================================== ###
#%% 80940 vs 81084
ylog=True
mode='amp'
plot_spectral=False 
plot_pearson=True
add_pt_zero=False

a=CorrelationAnalysis(80940, numDemod=True, verbose=False, plot=False)
ax_low, ax_high,_ = a.plot_coherence_delta([4], plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero)

b=CorrelationAnalysis(81065, numDemod=True, verbose=False, plot=False)
ax_low, ax_high,_ = b.plot_coherence_delta([4],ax_low_f=ax_low, ax_high_f=ax_high, plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero)

c=CorrelationAnalysis(81084, numDemod=True, verbose=False, plot=False)
ax_low, ax_high,_ = c.plot_coherence_delta([4],ax_low_f=ax_low, ax_high_f=ax_high, plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero)

# ax_low.set_xlim(-0.02, 0.02)
# ax_high.set_xlim(-0.02, 0.02)
# ax_low.set_ylim(0.01, 1)
# ax_high.set_ylim(0.01, 1)
ax_low.set_title('80940 ; 81084')
ax_high.set_title('80940 ; 81084')


#%% 80940 vs 81069
ylog=True
mode='amp'
plot_spectral=False 
plot_pearson=True
add_pt_zero=False

a=CorrelationAnalysis(80940, numDemod=True, verbose=False, plot=False)
ax_low, ax_high,_ = a.plot_coherence_delta([8, 9], plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero)

b=CorrelationAnalysis(81069, numDemod=True, verbose=False, plot=False)
ax_low, ax_high,_ = b.plot_coherence_delta([4, 5],ax_low_f=ax_low, ax_high_f=ax_high,
                                           plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero)
# ax_low.set_xlim(-0.02, 0.02)
# ax_high.set_xlim(-0.02, 0.02)
# ax_low.set_ylim(0.01, 1)
# ax_high.set_ylim(0.01, 1)
ax_low.set_title('80940 ; 81069')
ax_high.set_title('80940 ; 81069')

#%% 80949 vs 81069
ylog=True
mode='amp'
plot_spectral=False 
plot_pearson=True
add_pt_zero=False

a=CorrelationAnalysis(80940, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _ = a.plot_coherence_delta([4, 5, 6],plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero)
                                            

b=CorrelationAnalysis(81069, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _ = b.plot_coherence_delta([8, 9],ax_low_f=ax_low, ax_high_f=ax_high,
                                            plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero)
# ax_low.set_xlim(-0.02, 0.02)
# ax_high.set_xlim(-0.02, 0.02)
# ax_low.set_ylim(0.01, 1)
# ax_high.set_ylim(0.01, 1)
ax_low.set_title('80949 ; 81069')
ax_high.set_title('80949 ; 81069')

#%% 80949 vs 81084
ylog=True
mode='amp'
plot_spectral=False 
plot_pearson=True
add_pt_zero=False

a=CorrelationAnalysis(80949, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _ = a.plot_coherence_delta([8,9], plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero)

b=CorrelationAnalysis(81084, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _ = b.plot_coherence_delta([8,9],ax_low_f=ax_low, ax_high_f=ax_high,
                                            plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero)
# ax_low.set_xlim(-0.03, 0.015)
# ax_high.set_xlim(-0.03, 0.015)
# ax_low.set_ylim(0.6, 1)
# ax_high.set_ylim(0.6, 1)
ax_low.set_title('80949 ; 81084')
ax_high.set_title('80949 ; 81084')

# %%
ylog=True
mode='amp'
plot_spectral=False 
plot_pearson=True
add_pt_zero=False

a=CorrelationAnalysis(80949, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _  = a.plot_coherence_delta([4,5,6], plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero)


b=CorrelationAnalysis(81069, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _ = b.plot_coherence_delta([8,9],ax_low_f=ax_low, ax_high_f=ax_high,
                                           plot_pearson=plot_pearson,
                                           plot_spectral=plot_spectral, ylog=ylog,
                                           mode=mode, add_pt_zero=add_pt_zero)



# %%
#%% Mixed files technique for noise measurement
from diag_tcv.dbsAnalysis.correlationAnalysis.correlationAnalysis import CorrelationAnalysis

a=CorrelationAnalysis(80940, numDemod=True, verbose=False, plot=False)
b=CorrelationAnalysis(81084, numDemod=True, verbose=False, plot=False)
c=CorrelationAnalysis(81069, numDemod=True, verbose=False, plot=False)
d=CorrelationAnalysis(81087, numDemod=True, verbose=False, plot=False)

a_z_ref_list, a_z_hop_list, _, _ = a.get_normalized_data_isweep(4, ifreq_list=None, dtsart=400.e-6, dtend=100.e-6, ret=True)
b_z_ref_list, b_z_hop_list, _, _ = b.get_normalized_data_isweep(6, ifreq_list=None, dtsart=400.e-6, dtend=100.e-6, ret=True)
c_z_ref_list, c_z_hop_list, _, _ = c.get_normalized_data_isweep(9, ifreq_list=None, dtsart=400.e-6, dtend=100.e-6, ret=True)
d_z_ref_list, d_z_hop_list, _, _ = d.get_normalized_data_isweep(5, ifreq_list=None, dtsart=400.e-6, dtend=100.e-6, ret=True)

dt=a.processedData['sweep'+str(4)]['dt']

full_coherence_analysis(a_z_ref_list[5], c_z_hop_list[8], dt, nperseg=1024, noverlap=512, window=None, remove_mean=True, plot=False, verbose=False)






#%% Test amplitude vs full

ylog=True

a=CorrelationAnalysis(80949, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _ = a.plot_coherence_delta([8], plot_pearson=False, plot_spectral=True, ylog=ylog,add_pt_zero=True, mode='full')

a=CorrelationAnalysis(80949, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _ = a.plot_coherence_delta([8],ax_low_f=ax_low, ax_high_f=ax_high, plot_pearson=False, plot_spectral=True, ylog=ylog,add_pt_zero=True, mode='amp')
# %%
shot = 81065
ylog=True
add_pt_zero=False
plot_pearson=True
plot_spectral=False
isweep_list = [4]

mode='amp'
a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _ = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson, plot_spectral=plot_spectral, ylog=ylog, mode=mode, add_pt_zero=add_pt_zero)

mode="full"
a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _ = a.plot_coherence_delta(isweep_list,ax_low_f=ax_low, ax_high_f=ax_high, plot_pearson=plot_pearson, plot_spectral=plot_spectral, ylog=ylog, mode=mode,add_pt_zero=add_pt_zero)

mode="real"
a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _ = a.plot_coherence_delta(isweep_list,ax_low_f=ax_low, ax_high_f=ax_high, plot_pearson=plot_pearson, plot_spectral=plot_spectral, ylog=ylog, mode=mode, add_pt_zero=add_pt_zero)

mode="imag"
a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _ = a.plot_coherence_delta(isweep_list,ax_low_f=ax_low, ax_high_f=ax_high, plot_pearson=plot_pearson, plot_spectral=plot_spectral, ylog=ylog, mode=mode, add_pt_zero=add_pt_zero)

mode="phase"
a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
ax_low, ax_high, _ = a.plot_coherence_delta(isweep_list,ax_low_f=ax_low, ax_high_f=ax_high, plot_pearson=plot_pearson, plot_spectral=plot_spectral, ylog=ylog, mode=mode, add_pt_zero=add_pt_zero)

# %%
import diag_tcv.dbsAnalysis.correlationAnalysis.correlationFunctions as cF
import numpy as np 
import matplotlib.pyplot as plt
import scipy

# from dataAnalysis.utils.utils import get_closest_ind, normalize_array_1d
from dataAnalysis.utils.plot_utils import plot_1d
# from dataAnalysis.spectral.spectralAnalysis import custom_csd, custom_time_coherence

shot = 81069
ylog=True
add_pt_zero=False
plot_pearson=True
plot_spectral=False
isweep_list = [4]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)
        
ax_low, ax_high, _ = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson, plot_spectral=plot_spectral, ylog=ylog, mode=mode, add_pt_zero=add_pt_zero)


#%%
ifreq_list = None
a.get_normalized_data_isweep(isweep=4, ifreq_list=ifreq_list, dtsart=400.e-6, dtend=100.e-6)

nperseg=1024
noverlap=512
window=None
remove_mean=True

list_freq= [5]

for i in list_freq:
    zref = a.processedData['sweep'+str(4)]['z_list_ref'][i]
    zhop = a.processedData['sweep'+str(4)]['z_list_hop'][i]
    dt=a.processedData['sweep'+str(4)]['dt']

    zref_norm, zhop_norm = cF.normalize_signals(zref, zhop)

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
    tcorr_spec = scipy.signal.correlation_lags(1024, 1024, mode='same')*dt
    tcorr_scipy, corr_scipy = cF.scipy_correlation_function(zref_norm, zhop_norm, dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean)
   
    
    fig, ax = plot_1d([], [], grid=True)
    ax.plot(tcorr_spec*1e6,(corr), color='darkorchid', label='correlation function', markersize=5, markevery=20)
    # ax.plot(tcorr_spec*1e6,np.sqrt(corr.real**2+corr.imag**2), color='indigo', label='amplitude', markersize=5, markevery=20)
    ax.plot(tcorr_spec*1e6, corr_scipy, color='red', label='scipy correlation function', markersize=5, markevery=20)
    ax.set_xlabel(r'delay $[\mu_s]$')
    ax.set_ylabel('correlation')
    plt.legend()
    plt.title('Time correlation function')
    # plt.xlim(-20, 20)
    
    
    fig, ax = plot_1d([], [], grid=True)
    ax.plot(fref,(spectral_coh), color='green', label='spec correlation function', markersize=5, markevery=20)
    ax.set_xlabel(r'f [Hz]')
    ax.set_ylabel('correlation')
    plt.legend()
    plt.title('Spectral correlation function')
    # plt.xlim(-20, 20)
            
#%%
import scipy.signal as signal
n=len(zref)
corr_array = signal.correlate(zref_norm, zhop_norm, mode='same') / np.sqrt(signal.correlate(zref_norm,zref_norm, mode='same')[int(n/2)] * signal.correlate(zhop_norm, zhop_norm, mode='same')[int(n/2)])

delay_arr = signal.correlation_lags(n, n, mode='same')*dt
plot_1d(corr_array, delay_arr, markevery=25, markersize=5)


      
# %%
