#%% 

'''
Objective:
    - Functions to calculate velocity profiles
'''

#General imports
import numpy as np
import matplotlib.pyplot as plt
# from collections import Counter

from DBS.analysis import DBS_Profile

from dataAnalysis.utils.plot_utils import plot_1d, my_legend

def get_velocity_prof_object(shot, isweep_list, xmode=1, channelvals=[2]):
    '''
    Objective:
        - Get velocity profile for a given shot, isweep, xmode and channel
    Note: 
        - You have to perform the beamtracing before  
        - You have to save the velocity profile from the GUI before 
    Inputs:
        - shot: int
        - isweep_list: list of int
        - xmode: int
        - channelvals: list of int
    Outputs:
        - prof: DBS_Profile object => contains position of measure fDop and velocity
    '''
    
    machine = 'tcv'
    args = (machine, shot, isweep_list, xmode, channelvals)
    prof = DBS_Profile(*args)
    
    # prof.sort_values(by='rho_psi', ascending=True, inplace=True)
    
    return prof



def plot_velocity_prof(prof, ax=None, ifreq_list=None, sort=False, validated=False):
    '''
    Objective:
        - Plot velocity profile
    '''    
    
    if sort:
        prof.sort_values(by='rho_psi', ascending=True, inplace=True)
    
    nbsweep = len(prof.isweep.unique())
    colors = plt.cm.Dark2(np.linspace(0,1,nbsweep+2))
    marker_list = ['o', '^', 's', 'D', 'v',  '<', '>', 'p', 'h', 'H', 'd', 'P', 'X']
    
    if ifreq_list is None:
        ifreq_list = prof.ifreq.unique()
    
    
    #PLOT
    if ax is None:
        fig ,ax = plot_1d([], [], grid=True)
    for i in range(nbsweep):
        prof_i = prof[prof.isweep == prof.isweep.unique()[i]]
        
        ifreqlow = int(ifreq_list[0])
        ifreqhigh = int(ifreq_list[-1])
        
        if validated:
            rho_psi = prof_i[prof_i.validated==1].rho_psi.values[ifreqlow:ifreqhigh+1]
            v_perp = prof_i[prof_i.validated==1].v_perp.values[ifreqlow:ifreqhigh+1]
            
            err = np.zeros((2, len(prof_i[prof_i.validated==1].v_perp.values[ifreqlow:ifreqhigh+1])))
            err[1] = prof_i[prof_i.validated==1].dv_perp_up.values[ifreqlow:ifreqhigh+1]
            err[0] = prof_i[prof_i.validated==1].dv_perp_low.values[ifreqlow:ifreqhigh+1]
        
        else:
                
            rho_psi = prof_i.rho_psi.values[ifreqlow:ifreqhigh+1]
            v_perp = prof_i.v_perp.values[ifreqlow:ifreqhigh+1]
            
            err = np.zeros((2, len(prof_i.v_perp.values[ifreqlow:ifreqhigh+1])))
            err[1] = prof_i.dv_perp_up.values[ifreqlow:ifreqhigh+1]
            err[0] = prof_i.dv_perp_low.values[ifreqlow:ifreqhigh+1]
        
        
        ax.errorbar(rho_psi, v_perp, err, color=colors[i], marker=marker_list[i], linestyle='', label='isweep {}'.format(prof_i.isweep.unique()[0]))
    
    my_legend(ax, loc='upper right')
    ax.set_xlabel(r'$\rho_\psi$')
    ax.set_ylabel(r'$v_\perp$ [km/s]')
    ax.set_title('#{}'.format(prof.shot.values[0]))
    fig.tight_layout()
    
    
# %%
