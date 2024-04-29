#%%
#=== general info ===#
#Author: Olivier Panico
#contact: olivier.panico@free.fr

#Goal: performs the correlation of two Doppler channels

#=== imports ===#
#General imports
import numpy as np
import matplotlib.pyplot as plt
#Local imports
from DBS.io.interface import DataInterface
from dataAnalysis.utils.plot_utils import plot_1d
from dataAnalysis.utils.utils import get_closest_ind
from dataAnalysis.spectral.spectralAnalysis import custom_coherence, custom_time_coherence


def get_raw_signals(shot, isweep, ifreq_list, numDemod=False, verbose=False):
    '''
    WRAPPER FUNCTION
    Get the raw DBS signals (t, I, Q) for both channels
    Check the number of points per frequency step then calls get_raw_signals_diff_nbpts or get_raw_signals_same_nbpts
    depending on the result
    '''
    if numDemod:
        data_ref = DataInterface(shot=shot, isweep=isweep, channelval=3, machine='tcv', verbose=1)
        data_hop = DataInterface(shot=shot, isweep=isweep, channelval=4, machine='tcv', verbose=1)
    else:
        data_ref = DataInterface(shot=shot, isweep=isweep, channelval=1, machine='tcv', verbose=1)
        data_hop = DataInterface(shot=shot, isweep=isweep, channelval=2, machine='tcv', verbose=1)
        
    
    dtStep_ref = data_ref.params.dtStep #size of a frequency step
    dtAcq_ref = data_ref.params.dtAcq #acquisition time step

    dtStep_hop = data_hop.params.dtStep
    dtAcq_hop = data_hop.params.dtAcq
    
    assert(dtAcq_ref == dtAcq_hop)
    #We check whether both channels have the same nbts of pts per frequency step
    nbpts_per_freq_ref = int(dtStep_ref*1e-3/dtAcq_ref) #nb of points on a freq plateau
    nbpts_per_freq_hop = int(dtStep_hop*1e-3/dtAcq_hop)
    
    if nbpts_per_freq_ref==nbpts_per_freq_hop:
        if verbose:
            print(' === both ref & hop have same nbpts per freq === ')
        rawSigDic = get_raw_signals_same_nbpts(shot=shot, isweep=isweep, ifreq_list=ifreq_list, numDemod=numDemod, verbose=verbose)
    
    else:
        if verbose:
            print(' === ref & hop have different nbpts per freq === ')
        rawSigDic = get_raw_signals_diff_nbpts(shot=shot, isweep=isweep, ifreq_list=ifreq_list, numDemod=numDemod, verbose=verbose)

    return rawSigDic


def get_raw_signals_diff_nbpts(shot, isweep, ifreq_list, numDemod=False, verbose=False):
    '''
    Get the raw DBS signals (t, I, Q) for both channels
    case if channels have different number of points per frequency steps 
    => In this case we assume that reference channel is on long times and hopping is on short times
        The reference is time splitted so that it is decomposed on time windows of the same time as the hopping channel
    '''
    
    rawSigDic = dict()
    
    if numDemod:
        data_ref = DataInterface(shot=shot, isweep=isweep, channelval=3, machine='tcv', verbose=1)
        data_hop = DataInterface(shot=shot, isweep=isweep, channelval=4, machine='tcv', verbose=1)
    else:
        data_ref = DataInterface(shot=shot, isweep=isweep, channelval=1, machine='tcv', verbose=1)
        data_hop = DataInterface(shot=shot, isweep=isweep, channelval=2, machine='tcv', verbose=1)
  
    dtStep_ref = data_ref.params.dtStep #size of a frequency step
    dtAcq_ref = data_ref.params.dtAcq #acquisition time step

    dtStep_hop = data_hop.params.dtStep
    dtAcq_hop = data_hop.params.dtAcq
    
    assert(dtAcq_ref == dtAcq_hop)
    
    nbpts_per_freq_ref = int(dtStep_ref*1e-3/dtAcq_ref) #nb of points on a freq plateau
    nbpts_per_freq_hop = int(dtStep_hop*1e-3/dtAcq_hop)
    
    NbStep_ref = data_ref.params.NbStep #nb of frequency steps
    NbStep_hop = data_hop.params.NbStep
    nbhop_per_ref = NbStep_hop/NbStep_ref #nb of different hop steps in a single ref step
    
    #Initialize arrays
    nfreq = len(ifreq_list)
    nt = nbpts_per_freq_hop
    allfreq_ref, t_allfreq, I_allfreq_ref, Q_allfreq_ref = np.zeros((nfreq)), np.zeros((nfreq, nt)), np.zeros((nfreq, nt)) , np.zeros((nfreq, nt)) 
    allfreq_hop, t_allfreq_hop, I_allfreq_hop, Q_allfreq_hop = np.zeros((nfreq)), np.zeros((nfreq, nt)), np.zeros((nfreq, nt)) , np.zeros((nfreq, nt)) 

        
    for i, ifreq_loc in enumerate(ifreq_list):
        ifreq_loc = int(ifreq_loc)
        
        #We find in the reference the starting point of the hoping signal
        ifreq_ref = int(np.ceil(ifreq_loc / nbhop_per_ref)) # ifreq_ref that corresponds to ifreq_hop 
        print('ifreq_ref : ', ifreq_ref)
        _t_ref, _I_ref, _Q_ref = data_ref.get_signal(ifreq_ref) #loading ifreq_ref
        _t_hop, _I_hop, _Q_hop = data_hop.get_signal(ifreq=ifreq_loc) #loading ifreq_hop

        ind_timematch_ref_hop = get_closest_ind(_t_ref, _t_hop[0]) #finding the matching time condition
        print('time match condition : ', ind_timematch_ref_hop)
        
        print('nbpts_per_freq_hop : ', nbpts_per_freq_hop)

        #selects the time window that corresponds between the two signals
        _t_ref_match = _t_ref[ind_timematch_ref_hop: ind_timematch_ref_hop+ nt]
        _I_ref_match = _I_ref[ind_timematch_ref_hop: ind_timematch_ref_hop+ nt]
        _Q_ref_match = _Q_ref[ind_timematch_ref_hop: ind_timematch_ref_hop+ nt]
        
        #Assert the matching condition is respected
        assert( (_t_ref_match - _t_hop).all() == 0)
        
        #Store into arrays (ifreq_list, time)
        allfreq_ref[i] = data_ref.params.F[ifreq_ref]
        t_allfreq[i,:] = _t_ref_match
        I_allfreq_ref[i, :] = _I_ref_match
        Q_allfreq_ref[i, :] = _Q_ref_match

        allfreq_hop[i] = data_hop.params.F[ifreq_loc]
        t_allfreq_hop[i,:] = _t_hop
        I_allfreq_hop[i, :] = _I_hop
        Q_allfreq_hop[i, :] = _Q_hop

    assert((t_allfreq_hop-t_allfreq).all() == 0)
    
    rawSigDic['freq_list_ref'] = allfreq_ref
    rawSigDic['freq_list_hop'] = allfreq_hop
    rawSigDic['t'] = t_allfreq
    
    rawSigDic['I_list_ref'] = I_allfreq_ref
    rawSigDic['Q_list_ref'] = Q_allfreq_ref
    
    rawSigDic['I_list_hop'] = I_allfreq_hop
    rawSigDic['Q_list_hop'] = Q_allfreq_hop
    
    rawSigDic['params_ref'] = data_ref.params
    rawSigDic['params_hop'] = data_hop.params
    
    return rawSigDic




def get_raw_signals_same_nbpts(shot, isweep, ifreq_list, numDemod=False, verbose=False):
    '''
    Get the raw DBS signals (t, I, Q) for both channels
    case if both channels have the same number of points per frequency steps 
    => simple case as they are both on the same trigger / time frames
    '''
    
    rawSigDic = dict()
    
    nfreq = len(ifreq_list)
    
    if numDemod:
        data_ref = DataInterface(shot=shot, isweep=isweep, channelval=3, machine='tcv', verbose=1)
        data_hop = DataInterface(shot=shot, isweep=isweep, channelval=4, machine='tcv', verbose=1)
    else:
        data_ref = DataInterface(shot=shot, isweep=isweep, channelval=1, machine='tcv', verbose=1)
        data_hop = DataInterface(shot=shot, isweep=isweep, channelval=2, machine='tcv', verbose=1)
  
    for i, ifreq_loc in enumerate(ifreq_list):
        ifreq_loc = int(ifreq_loc)
        _t_ref, _I_ref, _Q_ref = data_ref.get_signal(ifreq_loc)
        _t_hop, _I_hop, _Q_hop = data_hop.get_signal(ifreq_loc)
        
        if i==0:
            nt = len(_t_ref)
            allfreq_ref, t_allfreq, I_allfreq_ref, Q_allfreq_ref = np.zeros((nfreq)), np.zeros((nfreq, nt)), np.zeros((nfreq, nt)) , np.zeros((nfreq, nt)) 
            allfreq_hop, t_allfreq_hop, I_allfreq_hop, Q_allfreq_hop = np.zeros((nfreq)), np.zeros((nfreq, nt)), np.zeros((nfreq, nt)) , np.zeros((nfreq, nt)) 
        # print(len(_t_ref))
        allfreq_ref[i] = data_ref.params.F[ifreq_loc]
        t_allfreq[i,:] = _t_ref
        I_allfreq_ref[i, :] = _I_ref
        Q_allfreq_ref[i, :] = _Q_ref

        allfreq_hop[i] = data_hop.params.F[ifreq_loc]
        t_allfreq_hop[i,:] = _t_hop
        I_allfreq_hop[i, :] = _I_hop
        Q_allfreq_hop[i, :] = _Q_hop
        
        assert((t_allfreq_hop-t_allfreq).all() == 0)

    rawSigDic['freq_list_ref'] = allfreq_ref
    rawSigDic['freq_list_hop'] = allfreq_hop
    rawSigDic['t'] = t_allfreq
    
    rawSigDic['I_list_ref'] = I_allfreq_ref
    rawSigDic['Q_list_ref'] = Q_allfreq_ref
    
    rawSigDic['I_list_hop'] = I_allfreq_hop
    rawSigDic['Q_list_hop'] = Q_allfreq_hop
    
    rawSigDic['params_ref'] = data_ref.params
    rawSigDic['params_hop'] = data_hop.params
        
    return rawSigDic



def get_signal(datainterface, I, Q, remove_edges=False):
    
    if remove_edges: #maye use the dedicated method instead
        x = I[5:-5]
        y = Q[5:-5]  
    else:
        x = I[:]
        y = Q[:]

    z = get_normalized_complex_signal(x, y, datainterface.params.phase_cor)


def correlation_channels(z_ref, z_hop, ifreq_list, mode='spectral'):
    '''
    if the shot is in correlation mode => correlation z_ref with np.conjugate(z_hop)
    if shot in normal mode => correlation z_ref with z_hop
    '''
    
    coh_max = []
    for i in range(len(ifreq_list)):
        if mode=='spectral':
            f,coh = custom_coherence((z_hop), (z_ref), nperseg=1024,noverlap=512, dt=dt, window='hanning',remove_mean=True)
            coh_max.apped(np.max(coh))
        elif mode == 'time':
            tcorr, corr = custom_time_coherence((z_hop), z_ref, nperseg=1024, noverlap=512)
            coh_max.append(np.max(corr))

    return np.array(coh_max)

#%%
from DBS.processing.sigprocessing import remove_signal_edges, get_normalized_complex_signal

if __name__ == '__main__':
    
    maxcoh = []
    for i in range(len(ifreq_list)):
        
        dt = t_allfreq[i,1] - t_allfreq[i,0]
        
        x_ref = I_allfreq_ref[i,5:-5]
        y_ref = Q_allfreq_ref[i,5:-5]
        
        x_hop = I_allfreq_hop[i, 5:-5]
        y_hop = Q_allfreq_hop[i, 5:-5]
        
        z_ref = get_normalized_complex_signal(x_ref, y_ref, 0)
        z_hop = get_normalized_complex_signal(x_hop, y_hop, 0)
        
        # print(len(z_ref), len(z_hop))
        # amp_ref = np.sqrt(I_allfreq_ref[i,:]**2 + Q_allfreq_ref[i,:]**2)
        # amp_hop = np.sqrt(I_allfreq_hop[i,:]**2 + Q_allfreq_hop[i,:]**2)
        
        # f,coh = custom_coherence((z_hop), (z_ref), nperseg=1024,noverlap=512, dt=dt, window='hanning',remove_mean=True)
        
        # f,coh = custom_coherence(np.conjugate(z_hop), z_ref, nperseg=1024,noverlap=512, dt=dt, window='hanning',remove_mean=True)
        # maxcoh.append(np.max(coh))
        
        tcorr,corr = custom_time_coherence((z_hop), z_ref, nperseg=4096, noverlap=2048)
        maxcoh.append(np.max(corr))
        
        


    #%% SHOT 79797
    data_ref = DataInterface(shot=79797, isweep=7, channelval=1, machine='tcv', verbose=1)
    data_hop = DataInterface(shot=79797, isweep=7, channelval=2, machine='tcv', verbose=1)
    
    # [20,21,22,23,24, 25, 26, 27, 28, 29 , 30]
    ifreq_list = np.linspace(20,39, 20)
    allfreq_ref, allfreq_hop, t_allfreq, I_allfreq_ref, Q_allfreq_ref, I_allfreq_hop, Q_allfreq_hop = get_raw_signals(79797, 7, ifreq_list)

    '''
    Pour le shot 79797 => bien mieux de correler z_hop et z_ref ensemble
    '''


    #%% SHOT 78549
    ifreq_list = np.linspace(1,20, 20)
    allfreq_ref, allfreq_hop, t_allfreq, I_allfreq_ref, Q_allfreq_ref, I_allfreq_hop, Q_allfreq_hop = get_raw_signals(78549, 2, ifreq_list)
    '''
    Pour le shot 78549 => bien mieux de correler np.conjugate(z_hop) et z_ref ensemble
'''