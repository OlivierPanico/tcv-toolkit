#%%
#=== general info ===#
#Author: Olivier Panico
#contact: olivier.panico@free.fr

#Goal: provides quick way of knowing the general information of
#       the studied shot. Gives infos on DBS channels frequency / positions 

#=== contains ===#



#=== imports ===#
#General imports
import numpy as np
import matplotlib.pyplot as plt
#Local imports
from DBS.io.interface import DataInterface
from dataAnalysis.utils.plot_utils import plot_1d

#=== remarks ===#
'''
If choc is not loaded from altair1
DBSsync numchoc --machine tcv

signals ordered by
isweep (nb of time the same pattern is repeated)
ifreq: frequency inside the pattern

example of utilisation: 
get_frequencies(shot=78549, isweep=2, ifreq=1, machine='tcv', fig=False)

'''


def plot_general_infos(shot, time, machine='tcv'):
    data_ref = DataInterface.from_time(shot=shot, time=time, channelval=1, machine='tcv', verbose=True)
    data_hop = DataInterface.from_time(shot=shot, time=time, channelval=2, machine='tcv', verbose=True)
    
    ### ===================== ###
    ### CHANNEL 1 : REFERENCE ###
    ### ===================== ###
    channelval = 1
    NbStep_ref = data_ref.params.NbStep
    dtStep_ref = data_ref.params.dtStep

    t0F_ref = data_ref.params.t0F
    t0F_seq1_ref = t0F_ref[NbStep_ref:2*NbStep_ref]

    ### =================== ###
    ### CHANNEL 2 : HOPPING ###
    ### =================== ###
    channelval = 2
    NbStep_hop = data_hop.params.NbStep
    dtStep_hop = data_hop.params.dtStep

    t0F_hop = data_hop.params.t0F
    t0F_seq1_hop = t0F_hop[NbStep_hop:2*NbStep_hop]
    
    ### ==================== ###
    ### FREQUENCY COMPARISON ###
    ### ==================== ###
    Nb_per_step = int(dtStep_ref//dtStep_hop) #Actually seems better to use dtStep_master / dtStep_slave
    # print('Nb_per_step = ', Nb_per_step)

    F_ref = np.zeros((NbStep_hop))

    for i in range(NbStep_ref):
        F_ref[i*Nb_per_step:(i+1)*Nb_per_step] = data_ref.params.F[i]


    fig, ax = plot_1d([],[], grid=True)
    # ax.plot(t0F_seq1_hop, F_ref, marker='+', label='f reference', color='blue')
    # ax.plot(t0F_seq1_hop, data_hop.params.F, marker='+', label='f hopping', color='red')
    ax.plot(F_ref, marker='+', label='f reference', color='blue')
    ax.plot(data_hop.params.F, marker='+', label='f hopping', color='red')
    plt.legend()


def plot_frequencies(shot, isweep):
    data_ref = DataInterface(shot=shot, isweep=isweep, channelval=1, machine='tcv', verbose=True)
    data_hop = DataInterface(shot=shot, isweep=isweep, channelval=2, machine='tcv', verbose=True)
    
    t0F_structured_ref = data_ref.params.t0F_structured[isweep-1,:]
    freq_ref = data_ref.params.F
    
    t0F_structured_hop = data_hop.params.t0F_structured[isweep-1,:]
    freq_hop = data_hop.params.F
    
    fig, ax = plot_1d([],[], grid=True)
    ax.plot(t0F_structured_ref, freq_ref, color='blue', marker='+')
    ax.plot(t0F_structured_hop, freq_hop, color='red', marker='+')
    










