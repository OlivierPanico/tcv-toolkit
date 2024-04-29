#%% 2024/04/29 meeting: how to compute the coherence 


from diag_tcv.dbsAnalysis.correlationAnalysis.correlationAnalysis import CorrelationAnalysis
from dataAnalysis.utils.plot_utils import plot_1d
from dataAnalysis.utils.utils import get_closest_ind, normalize_array_1d
from dataAnalysis.spectral.spectralAnalysis import custom_csd, custom_time_coherence

import numpy as np
import matplotlib.pyplot as plt
import scipy

#%% Step 0: loading data from DBS


correlObject = CorrelationAnalysis(80257)
isweep_list=[2]
ifreq_list=np.linspace(20,39,20)
for i,sweeploc in enumerate(isweep_list):
    correlObject.wrapper_coherence_analysis(sweeploc, ifreq_list, method='time')
   
# fig, ax = plot_1d([],[], grid=True)
# ax.plot(correlObject.rho_list_hop[20:40], correlObject.corrSigDic['coh_max_list'], marker='+', color='blue')
# ax.axvline(np.mean(correlObject.rho_list_ref[20:40]), color='black')
# ax.set_xlabel(r'$\Delta_\rho$')
# ax.set_ylabel('coherence')
# plt.yscale('log')
# plt.title('#{} ; sweep {}'.format(80257, 2))


#%% Step 1: taking two complex signals


zref = correlObject.corrSigDic['z_list_ref'][5]
zhop = correlObject.corrSigDic['z_list_hop'][5]
t_reduced_ref = correlObject.corrSigDic['t_reduced_list_ref'][5]
t_reduced_hop = correlObject.corrSigDic['t_reduced_list_hop'][5]

fig, ax = plot_1d([], [], grid=True)
ax.plot(t_reduced_hop*1000, zhop, color='red', label='zhop')
ax.plot(t_reduced_ref*1000, zref, color='blue', label='zref')
ax.set_xlabel('t [ms]')
ax.legend()


def normalize_complex_arr(a):
    a_oo = a - a.real.min() - 1j*a.imag.min() # origin offsetted
    return a_oo/np.abs(a_oo).max()

def normalize_complex_1d(a):
    a_r = normalize_array_1d(a.real)
    a_i = normalize_array_1d(a.imag)
    return a_r + 1j*a_i

# zref_norm = normalize_complex_arr(zref)
# zhop_norm = normalize_complex_arr(zhop)
zref_norm = normalize_complex_1d(zref)
zhop_norm = normalize_complex_1d(zhop)

fig, ax = plot_1d([], [], grid=True)
ax.plot(t_reduced_hop*1000, zhop_norm, color='red', label='normalized zhop')
ax.plot(t_reduced_ref*1000, zref_norm, color='blue', label='normalized zref')
ax.set_xlabel('t [ms]')
ax.legend()

dt=np.diff(t_reduced_ref)[0]
# %% Step 2: compute psd and csd


fref, psd_ref = custom_csd(zref_norm, zref_norm, dt=dt, nperseg=1024, noverlap=512, window=None, remove_mean=True)
fhop, psd_hop = custom_csd(zhop_norm, zhop_norm,dt=dt, nperseg=1024, noverlap=512, window=None, remove_mean=True)
fcsd, csd = custom_csd(zref_norm, zhop_norm, dt=dt, nperseg=1024, noverlap=512, window=None, remove_mean=True)

fig, ax = plot_1d([], [], grid=True)
ax.plot(fref/1000, psd_ref, color='blue', label='ref')
ax.plot(fhop/1000, psd_hop, color='red', label='hop')
ax.plot(fcsd/1000, abs(csd), color='green', label='abs(csd)')
plt.legend()
plt.xlabel('f [MHz]')
plt.ylabel('spectral density')
plt.yscale('log')



# %% Step 3: normalization of csd by psd to obtain spectral coherence


spectral_coh = abs(csd)**2/(psd_ref*psd_hop)

fig, ax = plot_1d([], [], grid=True)
ax.plot(fcsd/1000, spectral_coh, color='green', label='spectral coherence')
plt.legend()
plt.xlabel('f [MHz]')
plt.ylabel('coherence')
# plt.xlim(-10, 10)
# plt.yscale('log')

# %% Step 4: ifftshift + inverse Fourier transform to obtain the correlation function 
#Careful: scaling depends on the code you are using => it can preserve either the energy or the amplitude of the signal



csd_ifftshift = np.fft.ifftshift(csd)
corr_from_csd = np.fft.ifft(csd_ifftshift)/dt/2

psd_ref_ifftshift = np.fft.ifftshift(psd_ref)
corr_from_psd_ref = np.fft.ifft(psd_ref_ifftshift)/dt/np.sqrt(2)

psd_hop_ifftshift = np.fft.ifftshift(psd_hop)
corr_from_psd_hop = np.fft.ifft(psd_hop_ifftshift)/dt/np.sqrt(2)


fig, ax = plot_1d([], [], grid=True)
ax.plot(corr_from_csd, color='darkorchid')
# ax.plot(corr_from_psd_ref, color='blue')
# ax.plot(corr_from_psd_hop, color='red')



# %% Step 5: fftshift on the correlation function to replace peak at the center


corr = np.fft.fftshift(corr_from_csd)
corr_psd_ref = np.fft.fftshift(corr_from_psd_ref)
corr_psd_hop = np.fft.fftshift(corr_from_psd_hop)

# tcorr_spec = np.linspace(-511, 512, 1024)*dt
tcorr_spec = scipy.signal.correlation_lags(1024, 1024, mode='same')*dt

fig, ax = plot_1d([], [], grid=True)
ax.plot(tcorr_spec*1e6,(corr), color='darkorchid', label='correlation function')
ax.plot(tcorr_spec*1e6,np.sqrt(corr.real**2+corr.imag**2), color='indigo', label='amplitude')
ax.set_xlabel(r'delay $[\mu_s]$')
ax.set_ylabel('correlation')
plt.legend()
plt.title('Spectral correlation function')
plt.xlim(-20, 20)
# plt.ylim(0.9,1.1)


# fig, ax = plot_1d([], [], grid=True)
# ax.plot(tcorr*1e6, corr_psd_ref, color='blue')
# ax.plot(tcorr*1e6, corr_psd_hop, color='red')
# ax.set_xlabel(r'delay $[\mu_s]$')
# ax.set_ylabel('correlation function')
# plt.xlim(-20, 20)


# fig, ax = plot_1d([], [], grid=True)
# ax.plot(corr**2/(corr_psd_ref*corr_psd_hop))
# plt.xlim(300, 700)



#%% Alternate method: using Pearson coefficient and staying in real time space
tcorr, pearsonCorr = custom_time_coherence(zhop, zref, nperseg=1024, noverlap=512)

fig, ax = plot_1d([], [], grid=True)
ax.plot(tcorr*dt*1e6, pearsonCorr, color='teal', label='correlation function')
ax.plot(tcorr*dt*1e6, abs(pearsonCorr), color='darkslategrey', label='amplitude')
plt.xlim(-20, 20)
plt.legend()
plt.xlabel(r'delay $[\mu_s]$')
plt.ylabel('correlation')
plt.title('Pearson correlation function')



#%% Compare time method & spectral method

fig, ax = plot_1d([], [], grid=True)
ax.plot(tcorr*dt*1e6, pearsonCorr, color='teal', label='Pearson correlation')
ax.plot(tcorr*dt*1e6, abs(pearsonCorr), color='darkslategrey', label='Pearson amplitude')
ax.plot(tcorr_spec*1e6, (corr), color='darkorchid', label='Spectral correlation')
ax.plot(tcorr_spec*1e6,np.sqrt(corr.real**2+corr.imag**2), color='indigo', label='Spectral amplitude')
plt.legend()
plt.xlabel(r'delay $[\mu_s]$')
plt.ylabel('correlation')
plt.title('Comparison between time and spectral correlation')
plt.xlim(-20, 20)

# %% Test using scipy csd

import scipy

f, csd_from_scipy = scipy.signal.csd(zref_norm, zhop_norm, fs=1/dt, nperseg=1024, noverlap=512, scaling='spectrum')

f=np.fft.fftshift(f)
csd_from_scipy = np.fft.fftshift(csd_from_scipy)

fig, ax = plot_1d([], [], grid=True)
ax.plot(f/1e6,abs(csd_from_scipy))
plt.yscale('log')


csd_from_scipy_ifftshift = np.fft.ifftshift(csd_from_scipy)

corr_from_scipy = np.fft.ifft(csd_from_scipy_ifftshift)


fig, ax = plot_1d([], [], grid=True)
ax.plot(corr_from_scipy)
# plt.yscale('log')



# %%
