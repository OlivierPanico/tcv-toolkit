# %%
import scipy
import numpy as np
import matplotlib.pyplot as plt



from dataAnalysis.utils.plot_utils import plot_1d, my_text, my_legend

from diag_tcv.dbsAnalysis.correlationAnalysis.correlationAnalysis import CorrelationAnalysis
import diag_tcv.dbsAnalysis.correlationAnalysis.correlationFunctions as cF

shot = 81069
isweep = 4
ifreq = [10]
a=CorrelationAnalysis(shot, numDemod=True)
z_list_ref, z_list_hop, t_reduced_list_ref, t_reduced_list_hop = a.get_normalized_data_isweep(isweep, ifreq_list=ifreq,
                                                                                              dtsart=400.e-6, dtend=100.e-6, ret=True)
dt = t_reduced_list_ref[0][1]-t_reduced_list_ref[0][0]

nperseg = 1024
noverlap = 512
window=None
remove_mean=True
mode='full'

zref=z_list_ref[0]
zhop=z_list_hop[0]

 # prepare signals
zref_norm, zhop_norm = cF.normalize_signals(zref, zhop)
zref_norm, zhop_norm = cF.choose_mode(zref_norm, zhop_norm, mode=mode)

# fourier space
fref, psd_ref, fhop, psd_hop, fcsd, csd = cF.compute_psd_csd(zref_norm, zhop_norm, dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean)
spectral_coh = cF.compute_spectral_coherence(psd_ref, psd_hop, csd)
spectral_coh_smooth = scipy.signal.savgol_filter(spectral_coh, 51, 3)

# from fourier to real space
corr_from_csd = cF.compute_correlation_function(csd, dt, mode=mode)
corr_spec = cF.shift_correlation_function(corr_from_csd).real
tcorr_spec = scipy.signal.correlation_lags(nperseg, nperseg, mode='same')*dt
corr_spec_amp = np.sqrt(corr_spec.real**2+corr_spec.imag**2)
tcorr_spec_mus = tcorr_spec*1e6
# for tests (correlation directly in real space using scipy.signal.correlate)
tcorr_scipy, corr_scipy = cF.scipy_correlation_function(zhop_norm,zref_norm, dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean)
corr_scipy=corr_scipy.real
tcorr_scipy_mus = tcorr_scipy*1e6

fig, ax = plot_1d([], [], grid=True)
ax.plot(fref, abs(psd_ref), color='blue', label='psd ref', marker='')
ax.plot(fhop, abs(psd_hop), color='red', label='psd hop', marker='')
ax.plot(fcsd, abs(csd), color='green', label='csd', marker='')
ax.set_yscale('log')
ax.set_title('#{} ; isweep {} ; ifreq {}'.format(shot, isweep, ifreq))
my_legend(ax)

fig,ax = plot_1d([], [], grid=True)
ax.plot(fref, abs(spectral_coh), color='darkorchid', label='spectral coh', marker='')
ax.plot(fref, abs(spectral_coh_smooth), color='black', label='spectral coh smooth', marker='')
my_text(ax, 0.2, 0.85, r'$\max$(sc) = {:.3f}'.format(np.max(abs(spectral_coh))), color='darkorchid', fontsize='12')
my_text(ax, 0.2, 0.70, r'$\max$(sc_smoo) = {:.3f}'.format(np.max(abs(spectral_coh_smooth))), color='black', fontsize='12')
ax.set_xlabel(r'frequency [Hz]')
ax.set_ylabel('coherence')
my_legend(ax, fontsize='12')
ax.set_xlim(-3e6, 3e6)
plt.title('Spectral coherence')

fig, ax = plot_1d([], [], grid=True)
ax.plot(tcorr_spec_mus, corr_spec, color='black', label='correlation', marker='')
ax.plot(tcorr_spec*1e6, corr_spec_amp, color='red', label='amplitude', marker='')
ax.plot(tcorr_scipy_mus, corr_scipy, color='blue', label='scipy correlation', marker='')
ax.set_xlabel(r'time [$\mu_s$]')
ax.set_ylabel('correlation')
ax.set_xlim(-15, 15)
my_text(ax, 0.2, 0.85, r'$\max$(sc corr) = {:.3f}'.format(np.max(abs(corr_scipy))), color='blue', fontsize='12')
my_text(ax, 0.2, 0.7, r'$\max$(corr) = {:.3f}'.format(np.max(abs(corr_spec))), color='black', fontsize='12')
my_text(ax, 0.8, 0.25, r'delay = {:.2f} $\mu_s$'.format(tcorr_scipy_mus[np.argmax(abs(corr_scipy))]), color='blue', fontsize='12')
my_text(ax, 0.8, 0.1, r'delay = {:.2f} $\mu_s$'.format(tcorr_spec_mus[np.argmax(abs(corr_spec))]), color='black', fontsize='12')
my_legend(ax, fontsize='12')
plt.title('Correlation function')



fig, ax = plot_1d([], [], grid=True)
ax.set_yscale('log')
ax.plot(tcorr_spec*1e6, abs(corr_spec_amp), color='red', label='amplitude', marker='')
ax.plot(tcorr_scipy_mus, abs(corr_scipy), color='blue', label='scipy correlation', marker='')
ax.set_ylabel('correlation')
ax.set_xlim(-10, 10)
ax.set_ylim
my_legend(ax, fontsize='12')
plt.title('Correlation amplitude (log)')
# %%
