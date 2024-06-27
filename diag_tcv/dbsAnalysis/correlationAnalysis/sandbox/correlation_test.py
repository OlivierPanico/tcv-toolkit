
#%%
import numpy as np
import matplotlib.pyplot as plt

#If choc is not loaded from altair1
#DBSsync numchoc --machine tcv

from DBS.io.interface import DataInterface, DataInterface_ifreq
#isweep:index declenche (commence a 1)
#channelval: index DBS (1 or 2)
#shot, isweep, channelval, machine = 78541, 5, 1, 'tcv'

# shot, isweep, ifreq, channelval, machine = 78549, 1, 2, 1, 'tcv'
# data = DataInterface(shot, isweep=isweep,ifreq=ifreq, channelval=channelval, machine=machine)
# tmaster, xmaster, ymaster = data.get_signal()

# shot, isweep, channelval, machine = 78549, [], 2, 'tcv'
# data = DataInterface(shot, isweep=isweep, channelval=channelval, machine=machine)
# tslave, xslave, yslave = data.get_signal()

# shot = 78549
# isweep = 2
# machine = 'tcv'
# ifreq = 1

shot = 79797
isweep = 1
machine = 'tcv'
ifreq = 2

### ================== ###
### CHANNEL 1 : MASTER ###
### ================== ###
channelval = 1
data_master = DataInterface_ifreq(shot, isweep=isweep,ifreq=ifreq, channelval=channelval, machine=machine)
# print(data_master.params)
NbStep_master = data_master.params.NbStep
dtStep_master = data_master.params.dtStep

t0F_master = data_master.params.t0F
t0F_seq1_master = t0F_master[NbStep_master:2*NbStep_master]

# for i in range(1,6):
#     ifreq=i
#     print(ifreq)
#     data = DataInterface(shot, isweep=isweep,ifreq=ifreq, channelval=channelval, machine=machine)
#     tmaster, xmaster, ymaster = data.get_signal()
#     print(tmaster[0], tmaster[-1])

### ================= ###
### CHANNEL 2 : SLAVE ###
### ================= ###
channelval = 2
data_slave = DataInterface_ifreq(shot, isweep=isweep,ifreq=ifreq, channelval=channelval, machine=machine)
# print(data_slave.params)
NbStep_slave = data_slave.params.NbStep
dtStep_slave = data_slave.params.dtStep

t0F_slave = data_slave.params.t0F
t0F_seq1_slave = t0F_slave[NbStep_slave:2*NbStep_slave]

### ==================== ###
### FREQUENCY COMPARISON ###
### ==================== ###
Nb_per_step = int(dtStep_master//dtStep_slave) #Actually seems better to use dtStep_master / dtStep_slave
print('Nb_per_step = ', Nb_per_step)

F_master = np.zeros((NbStep_slave))

for i in range(NbStep_master):
    F_master[i*Nb_per_step:(i+1)*Nb_per_step] = data_master.params.F[i]

plt.figure()
plt.plot(t0F_seq1_slave, F_master, marker='+', label='f master')
plt.plot(t0F_seq1_slave, data_slave.params.F, marker='+', label='f slave')
plt.legend()

# plt.figure()
# plt.plot(t0F_seq1_master, data_master.params.F, marker='+')
# plt.plot(t0F_seq1_slave, data_slave.params.F, marker='+')

# %%
import numpy as np
import matplotlib.pyplot as plt

plt.figure();
plt.plot(tmaster,xmaster, color='blue')
plt.plot(tmaster,ymaster, color='red')

# plt.figure();
# plt.plot(tslave,xslave, color='blue')
# plt.plot(tslave,yslave, color='red')

# plt.figure()
# plt.hist(x, bins=1000, color='blue')
# plt.hist(y, bins=1000, color='red')



#%% Check whether we can go down to smaller analysis time

import numpy as np
import matplotlib.pyplot as plt

#If choc is not loaded from altair1
#DBSsync numchoc --machine tcv

from DBS.io.interface import DataInterface, DataInterface_ifreq


machine='tcv'
channelval=2
shot = 78549
isweep=1
ifreq=4
print(ifreq)
data = DataInterface_ifreq(shot, isweep=isweep,ifreq=ifreq, channelval=channelval, machine=machine)
par = data.params

t, x, y = data.get_signal()
print(t[0], t[-1])

plt.figure()
plt.plot(t, x, color='blue')
plt.plot(t,y, color='red')

#dt = 40 ns
#1ms = 25 000 dt 
#Je peux calculer la vitesse pour 1, 2, 3 ..., 10 ms et voir à quel point je me rapproche

# from DBS.processing.sigprocessing import remove_noisy_signal_edges, get_normalized_complex_signal
# remove noisy edges 0.4 ms at the beginning and 0.2 ms at the end:
# _t, _x = remove_noisy_signal_edges(t, x)
# _t, _y = remove_noisy_signal_edges(t, y)


# plt.figure()
# plt.plot(_t, _x, color='blue')
# plt.plot(_t,_y, color='red')

# z = get_normalized_complex_signal(_x,_y, par.phase_cor)
        
# plt.figure()
# plt.plot(_t, (z))

from DBS.processing.sigprocessing import init_specobjs_specific_time, perform_specobj_fits, show_spec, get_fDop_from_fit_results, make_title
sobjs = init_specobjs_specific_time(data, [ifreq], tinit=0.03, tfin=0.031)

for s in sobjs[:]:
    perform_specobj_fits(s)
    fig, ax = plt.subplots()
    plot_dict = show_spec(s, ax=ax)
    make_title(ax, data, s.header.ifreq)
    ax.set_xlim(-2,2)

    dic_fDop_fit_results = get_fDop_from_fit_results(s.fit_params, s.xscale) 
    ax.text(0.8, 0.8, dic_fDop_fit_results['fDop'], transform=ax.transAxes)

    z = s.z
    fig, ax = plt.subplots()
    ax.plot(z.real, color='blue')
    ax.plot(z.imag, color='red')

# 
# get_fDop_from_fit_results

print(get_fDop_from_fit_results(sobjs[0].fit_params, sobjs[0].xscale))

# %%

from DBS.io.interface import DataInterface_ifreq

machine='tcv'
channelval=2
shot = 78549
isweep=2
ifreq=8#4
fDop_list, dfDop_list, analysis_time = [], [], []

data = DataInterface_ifreq(shot, isweep=isweep,ifreq=ifreq, channelval=channelval, machine=machine)
par = data.params
t, x, y = data.get_signal()
for i in range(1,101):
    tinit=t[0]
    tfin = tinit + 0.0005 + 0.0001*i
    print(tinit, tfin)
    
    sobjs = init_specobjs_specific_time(data, [ifreq], tinit=tinit, tfin=tfin)
    for j, s in enumerate(sobjs):
        print(j)
        perform_specobj_fits(s)
        dic_fDop_fit_results = get_fDop_from_fit_results(s.fit_params, s.xscale) 
    fDop = dic_fDop_fit_results['fDop']
    dfDop = dic_fDop_fit_results['dfDop']

    analysis_time.append(tfin-tinit)
    fDop_list.append(fDop)
    dfDop_list.append(dfDop)

analysis_time = np.array(analysis_time)
fDop_list = np.array(fDop_list)
dfDop_list = np.array(dfDop_list)


#%% Idem mais sur master

from DBS.io.interface import DataInterface_ifreq


machine='tcv'
channelval=1
shot = 78549
isweep=2
ifreq=1#4
fDop_list, dfDop_list, analysis_time = [], [], []

data = DataInterface_ifreq(shot, isweep=isweep,ifreq=ifreq, channelval=channelval, machine=machine)
par = data.params
t, x, y = data.get_signal()
for i in range(1,21):
    tinit=t[0]
    tfin = tinit + 0.01*i
    print(tinit, tfin)
    
    sobjs = init_specobjs_specific_time(data, [ifreq], tinit=tinit, tfin=tfin)
    for j, s in enumerate(sobjs):
        print(j)
        perform_specobj_fits(s)
        dic_fDop_fit_results = get_fDop_from_fit_results(s.fit_params, s.xscale) 
    fDop = dic_fDop_fit_results['fDop']
    dfDop = dic_fDop_fit_results['dfDop']

    analysis_time.append(tfin-tinit)
    fDop_list.append(fDop)
    dfDop_list.append(dfDop)

analysis_time = np.array(analysis_time)
fDop_list = np.array(fDop_list)
dfDop_list = np.array(dfDop_list)



#%%
plt.figure()
plt.scatter(analysis_time, fDop_list)
plt.errorbar(analysis_time, fDop_list, yerr=dfDop_list, fmt="o")
# plt.figure()
# plt.scatter(analysis_time, dfDop_list)
# 
# 
# %% ====================================================
### ================================== ###
### LOADING DATA FOR CORRELATION TESTS ###
### ================================== ###

from DBS.processing.sigprocessing import remove_noisy_signal_edges, get_normalized_complex_signal

from DBS.io.interface import DataInterface_ifreq
def get_DBS_signal(shot, isweep, ifreq, channelval, machine='tcv'):
    data = DataInterface_ifreq(shot, isweep=isweep,ifreq=ifreq, channelval=channelval, machine=machine)
    par = data.params
    t, x, y = data.get_signal()
    z = get_normalized_complex_signal(x, y, par.phase_cor)

    f = par.F[ifreq-1]

    return t, z, f


def get_master_slave(shot, isweep, reference_ifreq):
    t_master, z_master, f_master = get_DBS_signal(shot, isweep=isweep,ifreq=reference_ifreq, channelval=1, machine='tcv')
    

machine='tcv'
channelval=1
shot = 78549
isweep=2
ifreq=1
# ifreq=2

data_master = DataInterface_ifreq(shot, isweep=isweep,ifreq=ifreq, channelval=channelval, machine=machine)
par = data_master.params
t_master, x_master, y_master = data_master.get_signal()

z_master = get_normalized_complex_signal(x_master,y_master, par.phase_cor)

channelval=2
# ifreq=8
ifreq=9
# ifreq=24
data_slave = DataInterface_ifreq(shot, isweep=isweep,ifreq=ifreq, channelval=channelval, machine=machine)
par = data_slave.params
t_slave, x_slave, y_slave = data_slave.get_signal()

z_slave = get_normalized_complex_signal(x_slave,y_slave, par.phase_cor)




#%%

import sys
import numpy as np
import matplotlib.pyplot as plt
1

sys.path.append("/home/panico")
from dataAnalysis.spectral.spectralAnalysis import custom_csd, custom_coherence, custom_time_coherence
from dataAnalysis.utils.array_splitting import split_array_1d #, normalize_signal
from dataAnalysis.utils.utils import normalize_array_1d, normalize_array_2d, get_closest_ind

t_master, z_master, f_master = get_DBS_signal(shot=78549, isweep=2,ifreq=2, channelval=1, machine='tcv')

t_master_split = split_array_1d(t_master, 250000)
z_master_split = split_array_1d(z_master, 250000)

nperseg=1024
noverlap=nperseg//2

dt = t_master[1]-t_master[0]

allcoh = []
allf = []
for i, ifreq in enumerate(range(20,40)):
    # ifreq=i+1
    ifreq=ifreq+1
    print("ifreq =", ifreq, "i = ", i)
    t_slave, z_slave, f_slave = get_DBS_signal(shot=78549, isweep=2, ifreq=ifreq, channelval=2, machine='tcv')
    print('time matching : ', t_slave[0], t_master_split[i,0])

    x= np.array(z_slave)
    y= np.array(z_master_split[i,:])

    # x=normalize_signal(x)
    # y = normalize_signal(y)
    # if i==0:
    #     plt.figure()
    #     plt.plot(z_slave)
    #     plt.plot(z_master_split[i,:])

    # x = np.sqrt(z_slave.real**2 + z_slave.imag**2)
    # y = np.sqrt(z_master_split[i,:].real**2 + z_master_split[i,:].imag**2)

    # timecorr, corr = custom_time_coherence(x, y, nperseg=nperseg, noverlap=noverlap)
    f, coh = custom_coherence(x.imag, y, nperseg=nperseg, noverlap=noverlap, window='hanning', remove_mean=True, dt=dt)
    
    indmin = get_closest_ind(f, -8E6)
    indmax = get_closest_ind(f, +8E6)
    # allcoh.append(np.max(coh[1556:2540]))
    allcoh.append(np.max(coh[indmin:indmax]))
    allf.append(f_slave)
    plt.plot(f[indmin:indmax], coh[indmin:indmax], label = str(ifreq))

plt.legend()
# plt.xlim(-3E6, 3E6)

plt.figure()
plt.scatter(allf, allcoh)
plt.axvline(f_master, color='red')
#Master needs to be splitted in 20 different windows 


# %% =============================================================
### ================ ###
### TIME CORRELATION ###
### ================ ###


from scipy.signal import correlation_lags,windows, detrend, correlate



def correlate_time_signals(t, sig1, sig2):
    '''
    Calculate the correlation between two time signals 
    Positive delay implies that signal 1 is leading in front of signal 2
    returns: 
    - correlation (time delays) 
    - maximum of correlation
    - corresponding array of the time delays 
    - delays corresponding to maximum of correlation
    '''
    
    #Considering the timesteps are constant
    dt = t[1] - t[0]

    n = len(sig1)

    corr_array = np.zeros((n))
    delay_array = np.zeros((n))

    # delay_maxcorr_array = []
    # corr_max_array = []
    delay_arr = correlation_lags(n, n, mode='same')*dt
    
    # Windowing
    # window = windows.hann(n)
    
    # Detrending
    # _y1 = sig1
    # _y2 = sig2
    _y1 = detrend(sig1, type='constant')
    _y2 = detrend(sig2, type='constant')
    
    # if all(item == 0 for item in y1) or all(item == 0 for item in y2):
    #     corr_array[:] = np.zeros((np.shape(corr_array)))
    #     continue
    
    # _y1 *= window       # reduce effect of finite length 
    # _y2 *= window     


    corr_array[:] = correlate(_y2, _y1, mode='same') / np.sqrt(correlate(_y1, _y1, mode='same')[int(n/2)] * correlate(_y2, _y2, mode='same')[int(n/2)])
    
    delay_maxcorr =  delay_arr[np.where(abs(corr_array[:]) == np.max(abs(corr_array[:])))[0][0]]
    # delay_maxcorr = delay_arr[np.argmax(corr_array[:,i])] # delay corresponding to the maximum of correlation
    # delay_array[:,i] = delay_arr_1D 
    
    # delay_maxcorr_array.append(delay_maxcorr) 
    # corr_max_array.append(np.max(corr_array[:,i])) # maximum of the correlation for a given position
    
    return corr_array, delay_arr, delay_maxcorr


from numpy import split
nperseg=2000

x_split = split(x, nperseg)
y_split = split(y, nperseg)
t_split = split(t_slave, nperseg)
x_split=np.array(x_split)
y_split=np.array(y_split)
t_split=np.array(t_split)
n,p = np.shape(x_split)
print(n,p)
corr = np.zeros((nperseg))
for j in range(p):

    xloc = x_split[:,j]
    yloc = y_split[:,j]
    tloc = t_split[:,j]

    corr_array, delay_arr, delay_maxcorr = correlate_time_signals(tloc, xloc, yloc)
    corr += corr_array/p




# %% ============================================================

plt.figure()
plt.plot(delay_arr, corr)
plt.xlim(-0.0001, 0.0001)
# %% ============================================================
### ==================================== ###
### Cross spectral density and coherence ###
### ==================================== ###

import sys

sys.path.append("/home/panico")
from dataAnalysis.spectral.spectralAnalysis import custom_csd, custom_coherence, custom_time_coherence
from dataAnalysis.utils.array_splitting import split_array_1d #, normalize_signal
from dataAnalysis.utils.utils import normalize_array_1d, normalize_array_2d, get_closest_ind

from scipy.signal import welch, csd, coherence, correlate
import numpy as np
import matplotlib.pyplot as plt

dt = t_slave[1]-t_slave[0]
fs = 1/dt

t_master_split = split_array_1d(t_master, 250000)
z_master_split = split_array_1d(z_master, 250000)



x = z_master_split[7, :]
# x = z_master_split[8, :]
y = z_slave[:]

# x = (z_master[10000: 500000])
# y = (z_master[10300:500300])

# x = normalize_signal(x)
# y = normalize_signal(y)

nperseg = 4096

f, pxy = csd(x, y, fs=fs, nperseg=nperseg, nfft=None, noverlap=nperseg//2)
f, pxx = welch(x, fs=fs, nperseg=nperseg, nfft=None, noverlap=nperseg//2)
f, pyy = welch(y, fs=fs, nperseg=nperseg,nfft=None, noverlap=nperseg//2)

coh = abs(pxy)**2/(pxx*pyy)

# plt.figure()
# plt.semilogy(f, abs(pyy)**2)
# plt.semilogy(f, abs(pxx)**2)

# plt.figure()
# plt.plot(f, coh)

f, coh2 = coherence(x, y, fs=fs, nperseg=nperseg, noverlap=nperseg//2)
# plt.figure()
# plt.plot(f, coh2)

tcorr, mycorr = custom_time_coherence(x, y, nperseg=nperseg, noverlap=nperseg//2)

plt.figure()
plt.plot(tcorr, mycorr)
#%%
# plt.figure()
# plt.plot(timecoh)
# plt.xlim(0,10)
# plt.xlim(0.49E7,0.51E7)

# f,cxy = coherence(x, y, fs=fs, window='hann', nperseg=nperseg, nfft=nperseg)

# plt.figure()
# plt.plot(f,cxy)

# corr = correlate(x, y, mode='same')
# plt.figure()
# plt.plot(corr)


frq, spec_master = custom_csd(x, x, nfft=nperseg, noverlap=None, dt=dt, norm=False)

# plt.figure()
# plt.semilogy(frq, abs(spec))

frq, spec_slave = custom_csd(y, y, nfft=nperseg, noverlap=None, dt=dt, norm=False)
plt.figure()
plt.semilogy(frq, abs(spec_slave)**2)
plt.semilogy(frq, abs(spec_master)**2)
plt.ylim(1E-18, 1E-11)

frq, mycoh = custom_coherence(x, y, nfft=nperseg, noverlap=None, dt=dt, coherence=True, norm=False)

plt.figure()
plt.plot(frq, mycoh)


plt.figure()
plt.semilogy(frq, abs(abs(spec_slave)**2-abs(pyy)**2)/(abs(pyy)**2))
plt.semilogy(frq, abs(abs(spec_master)**2-abs(pxx)**2)/(abs(pxx)**2))

# print(t_corr[np.argmax(corr)])
# %%
# %%