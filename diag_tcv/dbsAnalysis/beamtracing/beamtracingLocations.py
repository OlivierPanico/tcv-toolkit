
#%% Compare Beamtracing 

import numpy as np
import matplotlib.pyplot as plt

from DBS.beamtracing.DBSbeam import _DBSbeam
from DBS.analysis import DBS_Profile
from DBS.beamtracing import DBSbeam

from diag_tcv.dbsAnalysis.correlationAnalysis.correlationAnalysis import CorrelationAnalysis
from diag_tcv.shotAnalysis.dischargeInfoMdsObject import plot_profiles_comparison_single_column

from dataAnalysis.utils.utils import get_closest_ind, normalize_array_1d, my_linearRegression
from dataAnalysis.utils.plot_utils import plot_1d, my_text, my_legend, prep_multiple_subplots

import diag_tcv.utils.myMplStyle



def plot_beam_on_prof(shot, isweep, ifreq_list, ax=None, cprof='black', cref='blue', chop='red', plot_k_perp=False):
    a = CorrelationAnalysis(shot, numDemod=True)
    a.get_raytracing_isweep(isweep, ifreq_list=ifreq_list)

    output_ref = a.processedData['sweep'+str(isweep)]['output_ref']
    output_hop = a.processedData['sweep'+str(isweep)]['output_hop']
    # ------------------------------------
    # locate the beam in the density profile
    ne = output_ref.plasma.neprofile
    rho_ne = output_ref.plasma.neprofilerho
    rho_beam_ref = output_ref.rho
    rho_beam_hop = output_hop.rho
    # densityprof.ne
    # rho_ne = densityprof.rho_psi
    # rho_beam = outp.rho

    rho_freq_ind_hop = []
    rho_freq_ind_ref = []
    for i in range(len(rho_beam_hop)):
        rho_freq_ind_hop.append(get_closest_ind(rho_ne, rho_beam_hop[i]))

    for i in range(len(rho_beam_ref)):
        rho_freq_ind_ref.append(get_closest_ind(rho_ne, rho_beam_ref[i]))


    rhoinit = 0.7
    rhofin = 1.2

    rhoinit_ind = get_closest_ind(rho_ne, rhoinit)
    rhofin_ind = get_closest_ind(rho_ne, rhofin)


    
    if ax is None:
        # fig, ax = plt.subplots()
        fig, ax = plot_1d([], [], grid=True)
    
    ax.plot(rho_ne[rhoinit_ind:rhofin_ind], ne[rhoinit_ind:rhofin_ind], color=cprof, marker='', label=r'#{} ; isweep {}'.format(shot, isweep))
    ax.scatter(rho_ne[rho_freq_ind_hop], ne[rho_freq_ind_hop], color=chop, marker='s', s=40)
    ax.scatter(rho_ne[rho_freq_ind_ref], ne[rho_freq_ind_ref], color=cref, marker='*', s=60)
    ax.axvline(1, color='black', linestyle='--')
    
    if plot_k_perp:
        for i in range(len(output_hop.k_perp)):
            ax.text(0.9, 0.95-0.05*i, r'$k_\perp =$ {:.2f}'.format(output_hop.k_perp[i]), horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)
    
    ax.set_ylabel('density (m$^{-3}$)')
    ax.set_xlabel(r'$\rho$')
    ax.set_title('beamtracing location on density profile')
    # ------------------------------------

    return ax


def plot_velocity(shot, isweep, xmode=1, channelvals=[4],machine='tcv', ax=None, cvel='blue', legend=True):
    
    
    print(' ---- DEBUG  #{} isweep {} ---- '.format(shot, isweep))
    
    args = (machine, shot, isweep, xmode, channelvals)
    prof = DBS_Profile(*args)

    prof.sort_values(by='rho_psi', ascending=True, inplace=True)
    if ax is None:
        fig ,ax = plot_1d([], [], grid=True)
    # err = np.zeros((2, len(prof.v_perp.values[:20])))
    # err[1] = prof.dv_perp_up.values[:20]
    # err[0] = prof.dv_perp_low.values[:20]
    # ax.errorbar(prof.rho_psi.values[:20], prof.v_perp.values[:20], err, linewidth=0.8,color=cvel)
    # ax.set_xlabel(r'$\rho_\psi$')
    # ax.set_ylabel(r'$v_\perp$ [km/s]')

    # # fig ,ax = plot_1d([], [], grid=True)
    # err = np.zeros((2, len(prof.v_perp.values[20:])))
    # err[1] = prof.dv_perp_up.values[20:]
    # err[0] = prof.dv_perp_low.values[20:]
    # ax.errorbar(prof.rho_psi.values[20:], prof.v_perp.values[20:], err, linewidth=0.8,color=cvel)
    # ax.set_xlabel(r'$\rho_\psi$')
    # ax.set_ylabel(r'$v_\perp$ [km/s]')

    
    err = np.zeros((2, len(prof.v_perp.values[:])))
    err[1] = prof.dv_perp_up.values[:]
    err[0] = prof.dv_perp_low.values[:]
    ax.errorbar(prof.rho_psi.values[:], prof.v_perp.values[:], err, linewidth=0.8,color=cvel, label=r'#{} ; isweep {}'.format(shot, isweep))
    if legend:
        my_legend(ax)
    ax.set_xlabel(r'$\rho_\psi$')
    ax.set_ylabel(r'$v_\perp$ [km/s]')
    ax.set_title('Perpendicular velocity')
    
    # rhoinit = 0.7
    # rhofin = 1.2
    # ax.set_xlim(rhoinit, rhofin)
    
    return ax


def wrapper_beam_vel_prof(shot_list, isweep_list, time_list, rhomin, rhomax, tavg, xmode, channelvals):
    ifreq_list = np.linspace(1,20,20, dtype=int)
    plot_k_perp = False
    
    cprof_list = ['blue', 'red', 'green', 'purple']
    cref_list = ['cyan', 'orange', 'lime', 'magenta']
    chop_list = ['blue', 'red', 'green', 'purple']
    
    fig, axs = prep_multiple_subplots(2,1,axgrid=[0,1])
    for i, shot in enumerate(shot_list):
        isweep = isweep_list[i]
        axs[0] = plot_beam_on_prof(shot, isweep, ifreq_list, ax=axs[0],cprof=cprof_list[i], cref=cref_list[i], chop=chop_list[i], plot_k_perp=False)
        axs[1] = plot_velocity(shot, isweep, xmode=xmode, channelvals=channelvals,machine='tcv', ax=axs[1], cvel=cprof_list[i], legend=False)

    my_legend(axs[0])
    axs[0].set_xlim(rhomin, rhomax)
    axs[1].set_xlim(rhomin, rhomax)

    plot_profiles_comparison_single_column(shot_list, time_list, plot_err=None, rhomin=rhomin, rhomax=rhomax, tavg=tavg)




#%% interactive plotting
from DBS.io.utils import run_line_magic
# run_line_magic('matplotlib', 'widget')

run_line_magic('matplotlib', 'inline')

#%% 81069 vs 80940
shot=81069
isweep = 4
ifreq_list = np.linspace(1,20,20, dtype=int)
plot_k_perp = False
ax = plot_beam_on_prof(shot, isweep, ifreq_list, ax=None,cprof='blue', cref='cyan', chop='blue', plot_k_perp=False)

shot=80940
isweep = 8
ifreq_list = np.linspace(1,20,20, dtype=int)
plot_k_perp = False
ax = plot_beam_on_prof(shot, isweep, ifreq_list, ax=ax,cprof='purple', cref='orange', chop='purple', plot_k_perp=False)
my_legend(ax)

shot=81069
isweep = 5
ifreq_list = np.linspace(1,20,20, dtype=int)
plot_k_perp = False
ax = plot_beam_on_prof(shot, isweep, ifreq_list, ax=None,cprof='blue', cref='cyan', chop='blue', plot_k_perp=False)

shot=80940
isweep = 9
ifreq_list = np.linspace(1,20,20, dtype=int)
plot_k_perp = False
ax = plot_beam_on_prof(shot, isweep, ifreq_list, ax=ax,cprof='purple', cref='orange', chop='purple', plot_k_perp=False)
my_legend(ax)

# %% 80940 vs 81084 : bad matching
shot=81084
isweep = 5
ifreq_list = np.linspace(1,20,20, dtype=int)
plot_k_perp = False
ax = plot_beam_on_prof(shot, isweep, ifreq_list, ax=None,cprof='blue', cref='cyan', chop='blue', plot_k_perp=False)

shot=80940
isweep = 4
ifreq_list = np.linspace(1,20,20, dtype=int)
plot_k_perp = False
ax = plot_beam_on_prof(shot, isweep, ifreq_list, ax=ax,cprof='purple', cref='orange', chop='purple', plot_k_perp=False)
my_legend(ax)




# %% 80949 vs 81069: bad match
shot=81069
isweep = 9
ifreq_list = np.linspace(1,20,20, dtype=int)
plot_k_perp = False
ax = plot_beam_on_prof(shot, isweep, ifreq_list, ax=None,cprof='blue', cref='cyan', chop='blue', plot_k_perp=False)

shot=80949
isweep = 6
ifreq_list = np.linspace(1,20,20, dtype=int)
plot_k_perp = False
ax = plot_beam_on_prof(shot, isweep, ifreq_list, ax=ax,cprof='purple', cref='orange', chop='purple', plot_k_perp=False)
my_legend(ax)



# %% 80949 vs 81065: good match
shot=81065
isweep = 9
ifreq_list = np.linspace(1,20,20, dtype=int)
plot_k_perp = False
ax = plot_beam_on_prof(shot, isweep, ifreq_list, ax=None,cprof='blue', cref='cyan', chop='blue', plot_k_perp=False)

shot=80949
isweep = 8
ifreq_list = np.linspace(1,20,20, dtype=int)
plot_k_perp = False
ax = plot_beam_on_prof(shot, isweep, ifreq_list, ax=ax,cprof='purple', cref='orange', chop='purple', plot_k_perp=False)
my_legend(ax)

# %% 80949 vs 81065: bad match

shot=81084
isweep = 8
ifreq_list = np.linspace(1,20,20, dtype=int)
plot_k_perp = False
ax = plot_beam_on_prof(shot, isweep, ifreq_list, ax=None,cprof='blue', cref='cyan', chop='blue', plot_k_perp=False)

shot=80949
isweep = 8
ifreq_list = np.linspace(1,20,20, dtype=int)
plot_k_perp = False
ax = plot_beam_on_prof(shot, isweep, ifreq_list, ax=ax,cprof='purple', cref='orange', chop='purple', plot_k_perp=False)
my_legend(ax)

# %%
from diag_tcv.shotAnalysis.dischargeInfoMdsObject import plot_profiles_comparison

shot_list=[81069, 80940]
time_list=[0.9, 1.7]

plot_profiles_comparison(shot_list, time_list, plot_err=None, rhomin=None, rhomax=None, tavg=None)



#%% plot velocity profile
# machine, shot, isweep_list, xmode, channelvals = 'tcv', 81069, [4], 1, [4]
# plot_velocity(shot, isweep_list, xmode=1, channelvals=[4],machine='tcv', ax=None, cvel='blue')




# %% 80940 (beg) VS 81065 (beg)
wrapper_beam_vel_prof(shot_list=[80940, 81065], isweep_list=[4, 4], time_list=[0.7, 0.7],
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])

#%% 80949 (end) VS 81065 (end)
wrapper_beam_vel_prof(shot_list=[80949, 81065], isweep_list=[8, 9], time_list=[1.5, 1.7],
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])


#%% 80949 (end) VS 81084 (end)
wrapper_beam_vel_prof(shot_list=[80949, 81084], isweep_list=[8, 9], time_list=[1.5, 1.7],
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])

#%% 80949 (beg) VS 81069 (beg)
wrapper_beam_vel_prof(shot_list=[80949, 81069], isweep_list=[4, 5], time_list=[0.7, 0.9],
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])

#%% 80940 (end) vs 81069 (end)
wrapper_beam_vel_prof(shot_list=[80940, 81069], isweep_list=[8, 9], time_list=[1.5, 1.7],
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])
# %% 80940 (beg) vs 81084 (beg)
wrapper_beam_vel_prof(shot_list=[80940, 81084], isweep_list=[4, 5], time_list=[0.7, 0.9],
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])


# %% 80949 (end) vs 81084 (end)
wrapper_beam_vel_prof(shot_list=[80949, 81084], isweep_list=[8, 9], time_list=[1.5, 1.7],
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])








# %% 
wrapper_beam_vel_prof(shot_list=[81087], isweep_list=[4], time_list=[0.7],
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])



#%% Week CW29
###############################################################
########################## WEEK CW29 ##########################
###############################################################
#%%
shot=82549
isweep = 3
ifreq_list='all'

plot_beam_on_prof(shot, isweep, ifreq_list, ax=None, cprof='black', cref='blue', chop='red', plot_k_perp=False)
   
wrapper_beam_vel_prof(shot_list=[shot], isweep_list=[isweep], time_list=[0.7],
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])

#%%
shot_list=[82549, 81087]
isweep_list = [3, 4]
time_list = [0.7, 0.7]
ifreq_list='all'

wrapper_beam_vel_prof(shot_list, isweep_list, time_list,
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])
#%%
shot_list=[82549, 82549]
isweep_list = [5, 6]
time_list = [1.4, 1.7]
ifreq_list='all'

wrapper_beam_vel_prof(shot_list, isweep_list, time_list,
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])


#%%
shot_list=[82547, 81084]
isweep_list = [5, 9]
time_list = [1.3, 1.7]
ifreq_list='all'

wrapper_beam_vel_prof(shot_list, isweep_list, time_list,
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])


#%%
shot_list=[82556, 81084]
isweep_list = [3, 5]
time_list = [0.7, 0.7]
ifreq_list='all'

wrapper_beam_vel_prof(shot_list, isweep_list, time_list,
                    rhomin=0.75, rhomax=1.05, tavg=0.2, xmode=1, channelvals=[4])
#%%
shot=82556
isweep = 3
ifreq_list='all'

plot_beam_on_prof(shot, isweep, ifreq_list, ax=None, cprof='black', cref='blue', chop='red', plot_k_perp=False)
   
    
   
#%% PREPARE BEAMTRACING
shotnb = 82556
isweep_list = [3,4,5,6]

for i, isweep in enumerate(isweep_list):
    _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=3, verbose=True, plot=True, load_if_existing=True)
    _DBSbeam('tcv',shot=shotnb, isweep=isweep, xmode=1, channelval=4, verbose=True, plot=True, load_if_existing=True)
    



