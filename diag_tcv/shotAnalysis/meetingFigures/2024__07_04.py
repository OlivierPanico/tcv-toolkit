#%%

from DBS.beamtracing.DBSbeam import _DBSbeam
from diag_tcv.dbsAnalysis.correlationAnalysis.correlationAnalysis import CorrelationAnalysis
import diag_tcv.dbsAnalysis.correlationAnalysis.correlationFunctions as cF
import numpy as np 
import matplotlib.pyplot as plt
import scipy
from dataAnalysis.utils.plot_utils import plot_1d, my_text, my_legend



# %%
                ### ====================================================== ###
                ### ================= HYDROGEN =========================== ###
                ### ====================================================== ###

shot = 80949
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=False
isweep_list = [4,5,6,7,8,9]
mode='amp'
a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
#%%
shot = 80940
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=False
isweep_list = [4,5,6,7,8,9]
mode='amp'
a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
#%%
shot = 80949
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
ind_aval_low=[0,8]
ind_turb_pos_low=[15,19]
ind_turb_high=[9,16]
ind_aval_high=None
ind_turb_pos_high=[15,19]

ax=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))


#%%
shot = 80949
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

ind_turb_low=[8,16]
ind_aval_low=[0,8]
ind_turb_pos_low=[15,19]
ind_turb_high=[9,16]
ind_aval_high=None
ind_turb_pos_high=[15,19]

ax=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))


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

ind_turb_low=[8,16]
ind_aval_low=[0,8]
ind_turb_pos_low=[15,19]
ind_turb_high=[9,16]
ind_aval_high=None
ind_turb_pos_high=[15,19]

ax=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))




#%%
shot = 80949
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
ind_turb_pos_low=[15,19]
ind_turb_high=[9,16]
ind_aval_high=None
ind_turb_pos_high=[15,19]

ax=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))



#%%
shot = 80949
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=False
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
ind_turb_high=[9,16]
ind_aval_high=None
ind_turb_pos_high=[15,19]

ax=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))


#%%

shot = 80949
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=False
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

ind_turb_low=[8,16]
ind_aval_low=[0,8]
ind_turb_pos_low=[15,19]
ind_turb_high=[9,16]
ind_aval_high=None
ind_turb_pos_high=[15,19]

ax=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))



#%%
### ========== ###
### SHOT 80940 ###
### ========== ###

shot = 80940
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
ind_aval_low=[0,8]
ind_turb_pos_low=[15,19]
ind_turb_high=[10,16]
ind_aval_high=[3,10]
ind_turb_pos_high=[15,19]

ax=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))

ax=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))


#%%
shot = 80940
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

ind_turb_low=[8,16]
ind_aval_low=[0,8]
ind_turb_pos_low=[15,19]
ind_turb_high=[9,16]
ind_aval_high=[3,9]
ind_turb_pos_high=[15,19]

ax=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))



#%%
shot = 80940
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

ind_turb_low=[8,16]
ind_aval_low=None
ind_turb_pos_low=[15,19]
ind_turb_high=[9,16]
ind_aval_high=[3,9]
ind_turb_pos_high=[15,19]

ax=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))



#%%
shot = 80940
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
ind_turb_pos_low=[15,19]
ind_turb_high=[9,16]
ind_aval_high=[3,9]
ind_turb_pos_high=[15,19]

ax=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))


#%%
shot = 80940
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
# xunit='rho_s'
xunit='cm'

ind_turb_low=[6,16]
ind_aval_low=[0,6]
ind_turb_pos_low=[15,19]
ind_turb_high=[9,16]
ind_aval_high=[3,9]
ind_turb_pos_high=[15,19]

ax=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))



#%%
shot = 80940
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [9]
mode='amp'

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=True)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)
xunit='rho_s'
# xunit='cm'

ind_turb_low=[6,16]
ind_aval_low=[0,6]
ind_turb_pos_low=[15,19]
ind_turb_high=[9,16]
ind_aval_high=[3,9]
ind_turb_pos_high=[15,19]

ax=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))






#%% 
### =============================== ###
### COMPARISON HYDROGEN / DEUTERIUM ###
### =============================== ###

### 80940 (end) vs 81069 (beginning) ###
mode='amp'
xunit='rho_s'
# xunit='cm'



shot = 80940
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [8]

a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)

ind_turb_low=None
ind_aval_low=None
ind_turb_pos_low=None
ind_turb_high=None
ind_aval_high=None
ind_turb_pos_high=None

ax_low=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True)
ax_low.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax_high=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True)
ax_high.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))


# shot = 80940
# ylog=True
# add_pt_zero=True
# plot_pearson=True
# plot_spectral=True
# isweep_list = [9]


# a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
# a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
# _, _, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
#                                                           plot_spectral=plot_spectral, ylog=ylog,
#                                                           mode=mode, add_pt_zero=add_pt_zero,
#                                                           retdata=True, plot_rho=False)

# ind_turb_low=None
# ind_aval_low=None
# ind_turb_pos_low=None
# ind_turb_high=None
# ind_aval_high=None
# ind_turb_pos_high=None

# ax_low=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
#                            ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True, ax=ax_low, color='green')
# # ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
# ax_high=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
#                            ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True, ax=ax_high, color='green')
# # ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))



shot = 81069
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list = [4]


a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
_, _, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)

ind_turb_low=None
ind_aval_low=None
ind_turb_pos_low=None
ind_turb_high=None
ind_aval_high=None
ind_turb_pos_high=None

ax_low=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True, ax=ax_low, color='firebrick')
# ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax_high=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True, ax=ax_high, color='firebrick')
# ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))



# shot = 81069
# ylog=True
# add_pt_zero=True
# plot_pearson=True
# plot_spectral=True
# isweep_list = [4]

# a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
# a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
# _, _, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
#                                                           plot_spectral=plot_spectral, ylog=ylog,
#                                                           mode=mode, add_pt_zero=add_pt_zero,
#                                                           retdata=True, plot_rho=False)


# ind_turb_low=None
# ind_aval_low=None
# ind_turb_pos_low=None
# ind_turb_high=None
# ind_aval_high=None
# ind_turb_pos_high=None

# ax_low=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
#                            ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True, ax=ax_low, color='purple')
# # ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
# ax_high=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
#                            ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True, ax=ax_high, color='purple')
# # ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))





#%%
#80940 (beginning) vs 810-ç (beginning)

mode='amp'
# xunit='rho_s'
xunit='cm'



shot_h = 80940
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list_h = [8]

a=CorrelationAnalysis(shot_h, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list_h[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
ax_low_f, ax_high_f, delta, data_h = a.plot_coherence_delta(isweep_list_h, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)

ind_turb_low=None
ind_aval_low=None
ind_turb_pos_low=None
ind_turb_high=None
ind_aval_high=None
ind_turb_pos_high=None

ax_low=cF.plot_correlation_slopes(data_h['delta_low_f'], data_h['val_low_f_pear'],data_h['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True, caption=False)
# ax_low.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax_high=cF.plot_correlation_slopes(data_h['delta_high_f'], data_h['val_high_f_pear'],data_h['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True, caption=False)
# ax_high.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))


# shot = 80940
# ylog=True
# add_pt_zero=True
# plot_pearson=True
# plot_spectral=True
# isweep_list = [5]


# a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
# a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
# _, _, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
#                                                           plot_spectral=plot_spectral, ylog=ylog,
#                                                           mode=mode, add_pt_zero=add_pt_zero,
#                                                           retdata=True, plot_rho=False)

# ind_turb_low=None
# ind_aval_low=None
# ind_turb_pos_low=None
# ind_turb_high=None
# ind_aval_high=None
# ind_turb_pos_high=None

# ax_low=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
#                            ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True, ax=ax_low, color='green')
# # ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
# ax_high=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
#                            ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True, ax=ax_high, color='green')
# # ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))



# shot = 81069
# ylog=True
# add_pt_zero=True
# plot_pearson=True
# plot_spectral=True
# isweep_list = [4]


# a=CorrelationAnalysis(shot, numDemod=True, verbose=False, plot=False)
# a.get_correlation_isweep(isweep_list[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
# _, _, delta, data = a.plot_coherence_delta(isweep_list, plot_pearson=plot_pearson,
#                                                           plot_spectral=plot_spectral, ylog=ylog,
#                                                           mode=mode, add_pt_zero=add_pt_zero,
#                                                           retdata=True, plot_rho=False)

# ind_turb_low=None
# ind_aval_low=None
# ind_turb_pos_low=None
# ind_turb_high=None
# ind_aval_high=None
# ind_turb_pos_high=None

# ax_low=cF.plot_correlation_slopes(data['delta_low_f'], data['val_low_f_pear'],data['rho_s_low_f'],
#                            ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True, ax=ax_low, color='firebrick')
# # ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
# ax_high=cF.plot_correlation_slopes(data['delta_high_f'], data['val_high_f_pear'],data['rho_s_high_f'],
#                            ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True, ax=ax_high, color='firebrick')
# # ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))



shot_d = 81069
ylog=True
add_pt_zero=True
plot_pearson=True
plot_spectral=True
isweep_list_d = [5]

a=CorrelationAnalysis(shot_d, numDemod=True, verbose=False, plot=False)
a.get_correlation_isweep(isweep_list_d[0], ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False)        
_, _, delta, data_d = a.plot_coherence_delta(isweep_list_d, plot_pearson=plot_pearson,
                                                          plot_spectral=plot_spectral, ylog=ylog,
                                                          mode=mode, add_pt_zero=add_pt_zero,
                                                          retdata=True, plot_rho=False)


ind_turb_low=None
ind_aval_low=None
ind_turb_pos_low=None
ind_turb_high=None
ind_aval_high=None
ind_turb_pos_high=None

ax_low=cF.plot_correlation_slopes(data_d['delta_low_f'], data_d['val_low_f_pear'],data_d['rho_s_low_f'],
                           ind_turb_low, ind_aval_low, ind_turb_pos_low, xunit=xunit, ylog=True, ax=ax_low,caption=False, color='purple')
# ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:3f}'.format(shot, isweep_list[0],data['rho_low_f']))
ax_high=cF.plot_correlation_slopes(data_d['delta_high_f'], data_d['val_high_f_pear'],data_d['rho_s_high_f'],
                           ind_turb_high, ind_aval_high, ind_turb_pos_high, xunit=xunit, ylog=True, ax=ax_high,caption=False, color='purple')
# ax.set_title(r'#{} ; isweep={} ; $\rho$ = {:.3f}'.format(shot, isweep_list[0], data['rho_high_f']))

ax_low.plot([], [], color='teal', label='hydrogen')
ax_low.plot([], [], color='purple', label='deuterium')
my_legend(ax_low, fontsize='12')
my_text(ax_low, 0.4, 0.9, r'$\rho_s$ = {:.2f} cm'.format(data_h['rho_s_low_f'][0]), fontsize=12, color='teal')
my_text(ax_low, 0.4, 0.75, r'$\rho_s$ = {:.2f} cm'.format(data_d['rho_s_low_f'][0]), fontsize=12, color='purple')
ax_low.set_title( r'#{} isweep {} vs #{} isweep {}'.format(shot_h, isweep_list_h[0], shot_d, isweep_list_d[0]))

ax_high.plot([], [], color='teal', label='hydrogen')
ax_high.plot([], [], color='purple', label='deuterium')
my_legend(ax_high, fontsize='12')
my_text(ax_high, 0.4, 0.9, r'$\rho_s$ = {:.2f} cm'.format(data_h['rho_low_f'][0]), fontsize=12, color='teal')
my_text(ax_high, 0.4, 0.75, r'$\rho_s$ = {:.2f} cm'.format(data_d['rho_low_f'][0]), fontsize=12, color='purple')
ax_high.set_title( r'#{} isweep {} vs #{} isweep {}'.format(shot_h, isweep_list_h[0], shot_d, isweep_list_d[0]))
# %%
shots=[80949, 80940, 81065, 81069]
from DBS.io.DBSsync import push_processed_data
for shot in shots:
    push_processed_data('tcv', shot, data='beamtracing', dest_hostname='altair1')
    push_processed_data('tcv', shot, data='fDop_estimation', dest_hostname='altair1')