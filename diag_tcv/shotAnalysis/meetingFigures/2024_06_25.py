#%%
from DBS.beamtracing.DBSbeam import _DBSbeam
from diag_tcv.dbsAnalysis.correlationAnalysis.correlationAnalysis import CorrelationAnalysis
import diag_tcv.dbsAnalysis.correlationAnalysis.correlationFunctions as cF
import numpy as np 
import matplotlib.pyplot as plt
import scipy
from dataAnalysis.utils.plot_utils import plot_1d, my_text, my_legend
from dataAnalysis.utils.utils import my_linearRegression
# from dataAnalysis.spectral.spectralAnalysis import custom_csd, custom_time_coherence
def plot_correlation_slopes(xdata, ydata, rho_s=None, ind_turb=None, ind_aval=None, ind_turb_pos=None, xunit='rho_s', ax=None, color='teal', ylog=True):
    '''
    xunit: 'rho_s' or 'cm'
    '''
    if xunit=='cm':
        xdata = xdata*rho_s
        
    if ax is None:
        fig, ax = plot_1d([], [], grid=True)
    
    plt.plot(xdata, (ydata), 'o', color=color)
    my_text(ax, 0.4, 0.9, r'$\rho_s$ = {:.2f} cm'.format(rho_s[0]), fontsize=12, color='k')

    if ind_turb is not None:
        ind_turb_min=ind_turb[0]    
        ind_turb_max=ind_turb[1]
        popt, perr, rsquared = my_linearRegression(xdata[ind_turb_min:ind_turb_max], np.log(ydata[ind_turb_min:ind_turb_max]), mode='affine')
        plt.plot(xdata[ind_turb_min:ind_turb_max], np.exp(popt[0]*xdata[ind_turb_min:ind_turb_max]+ popt[1]), 'r', marker='')
        if xunit=='rho_s':
            my_text(ax, 0.15, 0.9, r'$L_c \approx$ {:.2f} $\rho_s$'.format(1/popt[0]), fontsize=12, color='r')
        elif xunit=='cm':
            my_text(ax, 0.15, 0.9, r'$L_c \approx$ {:.2f} cm'.format(1/popt[0]), fontsize=12, color='r')
        
    if ind_aval is not None:
        ind_aval_min=ind_aval[0]
        ind_aval_max=ind_aval[1]
        popt, perr, rsquared = my_linearRegression(xdata[ind_aval_min:ind_aval_max], np.log(ydata[ind_aval_min:ind_aval_max]), mode='affine')
        plt.plot(xdata[ind_aval_min:ind_aval_max], np.exp(popt[0]*xdata[ind_aval_min:ind_aval_max] + popt[1]), 'g', marker='')
        if xunit=='rho_s':
            my_text(ax, 0.15, 0.75, r'$L_a \approx$ {:.2f} $\rho_s$'.format(1/popt[0]), fontsize=12, color='g')
        elif xunit=='cm':
            my_text(ax, 0.15, 0.75, r'$L_a \approx$ {:.2f} cm'.format(1/popt[0]), fontsize=12, color='g')

    if ind_turb_pos is not None:
        ind_turb_min=ind_turb_pos[0]
        ind_turb_max=ind_turb_pos[1]
        popt, perr, rsquared = my_linearRegression(xdata[ind_turb_min:ind_turb_max], np.log(ydata[ind_turb_min:ind_turb_max]), mode='affine')
        plt.plot(xdata[ind_turb_min:ind_turb_max], np.exp(popt[0]*xdata[ind_turb_min:ind_turb_max]+ popt[1]), 'b', marker='')
        if xunit=='rho_s':
            my_text(ax, 0.15, 0.6, r'$L_c \approx$ {:.2f} $\rho_s$'.format(1/popt[0]), fontsize=12, color='b')
        elif xunit=='cm':
            my_text(ax, 0.15, 0.6, r'$L_c \approx$ {:.2f} cm'.format(1/popt[0]), fontsize=12, color='b')

    if ylog:
        plt.yscale('log')
    
    plt.ylim(0.1,1.1)
    plt.ylabel('correlation')
    plt.axhline(1, color='black', linestyle='--')
    plt.axvline(0, color='black', linestyle='--')
    
    if xunit=='rho_s':
        plt.xlabel(r'$\Delta$ $[\rho_s]$')
    elif xunit=='cm':
        plt.xlabel(r'$\Delta$ $[cm]$')
    
    return ax


#%%
shot = 81069
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [4]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[8,16]
ind_aval_low=[0,10]
ind_turb_pos_low=[15,19]
ind_turb_high=[9,15]
ind_aval_high=[3,9]
ind_turb_pos_high=[16,21]

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))

#%%
shot = 81069
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [5]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
# xunit='rho_s'
xunit='cm'

ind_turb_low=[8,16]
ind_aval_low=[0,10]
ind_turb_pos_low=[15,19]
ind_turb_high=[9,16]
ind_aval_high=[3,9]
ind_turb_pos_high=[15,21]

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))

#%%
shot = 81069
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [6]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=True)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[7,16]
ind_aval_low=[0,9]
ind_turb_pos_low=[15,19]
ind_turb_high=[9,15]
ind_aval_high=[3,9]
ind_turb_pos_high=[16,21]

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))

#%%
shot = 81069
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [7]
mode='amp'

ifreq_list=np.linspace(0,19,20)

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=ifreq_list, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[7,16]
ind_aval_low=[0,9]
ind_turb_pos_low=[15,19]
ind_turb_high=None
ind_aval_high=[5,15]
ind_turb_pos_high=[15,21]

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))



#%%
shot = 81069
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [8]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[8,16]
ind_aval_low=[2,7]
ind_turb_pos_low=[15,21]
ind_turb_high=[10,16]
ind_aval_high=[3,10]
ind_turb_pos_high=[15,21]

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))

#%%
shot = 81069
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [9]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[7,16]
ind_aval_low=[1,7]
ind_turb_pos_low=[15,21]
ind_turb_high=[8,16]
ind_aval_high=[3,8]
ind_turb_pos_high=[15,21]

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))


#%%
### ============= ###
### SUMMARY 81069 ###
### ============= ###

isweep_list = [4,5,6,7,8,9]

Lc_low_list = [3.92, 4.58, 4.15, 3.96, 4.7, 4.15]
La_low_list = [9.08, 8.84, 10.26, 9.42, 12.92, 16.22]
rho_s_low_list = [0.09, 0.08, 0.09, 0.09, 0.09, 0.07]
rho_low_list = [0.962, 0.969, 0.931, 0.915, 0.948, 0.963]
rho_high_list = [0.848, 0.864, 0.8, 0.783, 0.932, 0.857]

Lc_high_list = [3.83, 3.56, np.nan, np.nan, 1.48, 1.85]
La_high_list = [5.09, 5.39, np.nan, 10.07, 7.79, 6.30]
rho_s_high_list = [0.14, 0.13, 0.13, 0.14, 0.13, 0.13]

fig, ax = plot_1d([], [], grid=True)
ax.plot(isweep_list, Lc_low_list, 'o', label='Lc low f', color='teal')

fig, ax = plot_1d([], [], grid=True)
ax.plot(isweep_list, La_low_list, 'o', label='La low f', color='firebrick')

fig, ax = plot_1d([], [], grid=True)
ax.plot(isweep_list, Lc_high_list, 'o', label='Lc high f', color='teal')
fig, ax = plot_1d([], [], grid=True)
ax.plot(isweep_list, La_high_list, 'o', label='La high f', color='firebrick')

### ============= ###


#%%
shot = 81065
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [4]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[8,15]
ind_aval_low=[0,8]
ind_turb_pos_low=[15,21]
ind_turb_high=[10,16]
ind_aval_high=[3,10]
ind_turb_pos_high=[15,21]

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))


#%%
shot = 81065
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [5]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[8,15]
ind_aval_low=[0,8]
ind_turb_pos_low=[15,21]
ind_turb_high=[10,16]
ind_aval_high=[3,10]
ind_turb_pos_high=[15,21]

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))


#%%
shot = 81065
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [6]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[11,16]
ind_aval_low=[1,8]
ind_turb_pos_low=[15,21]
ind_turb_high=[10,16]
ind_aval_high=[3,10]
ind_turb_pos_high=[15,19]

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))

#%%
shot = 81065
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [7]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[8,16]
ind_aval_low=[0,8]
ind_turb_pos_low=[15,21]
ind_turb_high=None
ind_aval_high=[4,15]
ind_turb_pos_high=None

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))


#%%
shot = 81065
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [8]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[8,16]
ind_aval_low=[0,8]
ind_turb_pos_low=[15,19]
ind_turb_high=[10,16]
ind_aval_high=[4,12]
ind_turb_pos_high=[15,19]

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))

#%%
shot = 81065
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [9]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[8,15]
ind_aval_low=[0,8]
ind_turb_pos_low=[15,21]
ind_turb_high=[10,16]
ind_aval_high=[3,10]
ind_turb_pos_high=[15,21]

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))

#%%
### ============= ###
### SUMMARY 81065 ###
### ============= ###

isweep_list = [4,5,6,7,8,9]

Lc_low_list = [5.08, 4.27, 5.34, 3.94, 5.15, 3.77]
La_low_list = [12, 10.81, 14.45, 16.02, 10.51, 8.04]
rho_s_low_list = [0.09, 0.09, 0.09, 0.09, 0.08, 0.08]
rho_low_list = [0.943,0.946, 0.925, 0.928, 0.97, 0.983]

Lc_high_list = [3.69, 3.13, np.nan, np.nan, 1.97, 2.29]
La_high_list = [5.01, 4.8, np.nan, 7.77, 6.09, 6.03]
rho_s_high_list = [0.14, 0.14, 0.14, 0.14, 0.13, 0.11]
rho_high_list = [0.815, 0.82, 0.789, 0.795, 0.86, 0.901]

#%% Correlation lengths as a function of heating 


heating = np.array([180, 180, 260, 260, 360, 360, 500, 500])

Lc_low = np.array([4.7, 4.15, 5.08, 4.27, 5.15, 3.77, 3.92, 4.58])
La_low = np.array([12.92, 16.22, 12, 10.81, 10.51, 8.04, 9.08, 8.84])
rho_s_low = np.array([0.09, 0.07, 0.09, 0.09, 0.08, 0.08, 0.09, 0.08])
rho_low = [0.948, 0.963, 0.943,0.946, 0.97, 0.983, 0.962, 0.969]

Lc_high = np.array([1.48, 1.85, 3.69, 3.13, 1.97, 2.29, 3.83, 3.56])
La_high = np.array([7.79, 6.30, 5.01, 4.8, 6.09, 6.03, 5.09, 5.39])
rho_s_high = np.array([0.13, 0.13, 0.14, 0.14, 0.13, 0.11, 0.14, 0.13])
rho_high = [0.932, 0.857, 0.815, 0.82, 0.86, 0.901, 0.848, 0.864]

fig, ax = plot_1d([], [], grid=True)
ax.plot(heating, Lc_low, 'o', label=r'$L_c$', color='teal')
ax.plot(heating, La_low, '^', label=r'$L_a$', color='firebrick')
ax.set_ylabel(r'Correlation length $[\rho_s]$')
ax.set_xlabel('Heating [kW]')
ax.set_title('Correlation lengths as a function of NBI-1')
my_legend(ax)

fig, ax = plot_1d([], [], grid=True)
ax.plot(heating, Lc_high, 'o', label=r'$L_c$', color='teal')
ax.plot(heating, La_high, '^', label=r'$L_a$', color='firebrick')
ax.set_ylabel(r'Correlation length $[\rho_s]$')
ax.set_xlabel('Heating [kW]')
ax.set_title('Correlation lengths as a function of NBI-1')
my_legend(ax)

# %%
                ### ====================================================== ###
                ### ================= HYDROGEN =========================== ###
                ### ====================================================== ###

#%%
shot = 80949
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [6]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[12,15]
ind_aval_low=[0,10]
ind_turb_pos_low=[15,19]
ind_turb_high=[12,16]
ind_aval_high=[3,12]
ind_turb_pos_high=[15,19]

ax=plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))

# %%
shots=[80949, 80940, 81065, 81069]
from DBS.io.DBSsync import push_processed_data
for shot in shots:
    push_processed_data('tcv', shot, data='beamtracing', dest_hostname='altair1')
    push_processed_data('tcv', shot, data='fDop_estimation', dest_hostname='altair1')