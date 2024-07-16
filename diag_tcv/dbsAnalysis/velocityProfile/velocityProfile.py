#%%

import matplotlib.pyplot as plt
from DBS.analysis import DBS_Profile
from DBS.beamtracing import DBSbeam


# machine, shot, isweep_list, xmode, channelvals = 'tcv', 80745, [9,11,13], 1, [1,2]

#to do: shot 80322 (ECRH / MIX (actually should be only NBI) / NBI)

# machine, shot, isweep_list, xmode, channelvals = 'tcv', 80376, [17], 1, [1]
machine, shot, isweep_list, xmode, channelvals = 'tcv', 81069, [4], 1, [4]

# if not done already: run beamtracing:
for ch in channelvals:
    for isweep in isweep_list:
        outp14, interface = DBSbeam(machine, shot, xmode=xmode, isweep=isweep, channelval=ch, verbose=True, plot=False, load_if_existing=True)

#%% plot 2 ch of 1 sweep

isweep=17
channelvals = [1]
color_list = ['xkcd:blue', 'xkcd:red']
fig, ax = plt.subplots()
for i,ch in enumerate(channelvals):
    args = (machine, shot, isweep, xmode, ch)
    prof = DBS_Profile(*args)
    plots_kwargs = {
        'sorted_along_rho': True,
        'show_fit': False,    
    }
    prof.plot(ax=ax, ls='', color=color_list[i], marker='+', lw=1,errorbars=True, **plots_kwargs)


ax.legend()
plt.grid()
plt.axvline(1.00, color='black')
plt.axhline(0, color='black')
# ax.set_xlim(0.92, 0.98)
# ax.set_ylim(-3.2, -0.9)
# %% plot same shot for different sweeps

# machine, shot, isweep_list, xmode, channelvals = 'tcv', 80745, [9,11, 13], 1, [1,2]
machine, shot, isweep_list, xmode, channelvals = 'tcv', 80376, [17,18], 1, [1]

color_list = ['blue', 'red', 'green', 'purple']
marker_list = ['+', 'o', 'x', 'd', '8']
isweep_to_plot = isweep_list

fig ,ax = plt.subplots(figsize=(9,6))

for i in range(len(isweep_to_plot)):
    args = (machine, shot, isweep_to_plot[i], xmode, channelvals)
    prof = DBS_Profile(*args)
    plots_kwargs = {
        'sorted_along_rho': True,
        'show_fit': False,    
    }
    prof.plot(ax=ax, color=color_list[i], marker=marker_list[i],ls='--', lw=1,errorbars=True, **plots_kwargs)
    
ax.grid()
plt.axvline(1.00, color='black')
plt.axhline(0, color='black')
# ax.set_xlim(0.65, 0.9)
# ax.set_ylim(-3.2, -0.9)
ax.legend(loc='upper left')

plt.tight_layout()



#%% PLOT ONE CH
import numpy as np
import matplotlib.pyplot as plt
from DBS.analysis import DBS_Profile
from DBS.beamtracing import DBSbeam

from dataAnalysis.utils.plot_utils import plot_1d, my_legend, my_text
import diag_tcv.utils.myMplStyle

machine, shot, isweep_list, xmode, channelvals = 'tcv', 81069, [4], 1, [3]

for i in range(len(isweep_list)):
    args = (machine, shot, isweep_list[i], xmode, channelvals)
    prof = DBS_Profile(*args)

fig ,ax = plot_1d([], [], grid=True)
err = np.zeros((2, len(prof.v_perp.values[:20])))
err[1] = prof.dv_perp_up.values[:20]
err[0] = prof.dv_perp_low.values[:20]
plt.errorbar(prof.t.values[:20], prof.v_perp.values[:20],yerr=err,linewidth=0.8,color='teal')
my_text(ax, 0.2, 0.9, r'$\rho$ = {:.2f} $\pm$ {:.2f}'.format(np.mean(prof.rho_psi.values[:20]), np.std(prof.rho_psi.values[:20])))
ax.set_xlabel('t [s]')
ax.set_ylabel(r'$v_\perp$ [km/s]')
ax.set_title('#{} ; isweep {} ; ch {}'.format(shot, isweep_list, channelvals))
plt.tight_layout()


fig ,ax = plot_1d([], [], grid=True)
err = np.zeros((2, len(prof.v_perp.values[20:])))
err[1] = prof.dv_perp_up.values[20:]
err[0] = prof.dv_perp_low.values[20:]
plt.errorbar(prof.t.values[20:], prof.v_perp.values[20:],yerr=err,linewidth=0.8,color='teal')
my_text(ax, 0.2, 0.9, r'$\rho$ = {:.2f} $\pm$ {:.2f}'.format(np.mean(prof.rho_psi.values[20:]), np.std(prof.rho_psi.values[20:])))
ax.set_xlabel('t [s]')
ax.set_ylabel(r'$v_\perp$ [km/s]')
ax.set_title('#{} ; isweep {} ; ch {}'.format(shot, isweep_list, channelvals))
plt.tight_layout()

#%% PLOT BOTH CH
import numpy as np
import matplotlib.pyplot as plt
from DBS.analysis import DBS_Profile
from DBS.beamtracing import DBSbeam

from dataAnalysis.utils.plot_utils import plot_1d, my_legend, my_text
import diag_tcv.utils.myMplStyle

machine, shot, isweep_list, xmode, channelvals = 'tcv', 81069, [4], 1, [3,4]

for i in range(len(isweep_list)):
    args = (machine, shot, isweep_list[i], xmode, channelvals)
    prof = DBS_Profile(*args)

prof3 = prof[prof.channel == channelvals[0]]
prof4 = prof[prof.channel == channelvals[1]]

fig ,ax = plot_1d([], [], grid=True)
err = np.zeros((2, len(prof.v_perp.values[:20])))
err[1] = prof.dv_perp_up.values[:20]
err[0] = prof.dv_perp_low.values[:20]

plt.errorbar(prof3.t.values[:20], prof3.v_perp.values[:20],yerr=err,linewidth=0.8,color='teal')
# my_text(ax, 0.2, 0.9, r'$\rho$ = {:.2f} $\pm$ {:.2f}'.format(np.mean(prof3.rho_psi.values[:20]), np.std(prof3.rho_psi.values[:20])))
my_text(ax, 0.2, 0.9, r'$\rho$ = {:.2f} - {:.2f} '.format(prof3.rho_psi.values[0], prof3.rho_psi.values[19]))
plt.errorbar(prof4.t.values[:20], prof4.v_perp.values[:20],yerr=err,linewidth=0.8,color='firebrick')
# my_text(ax, 0.2, 0.6, r'$\rho$ = {:.2f} $\pm$ {:.2f}'.format(np.mean(prof4.rho_psi.values[:20]), np.std(prof4.rho_psi.values[:20])))
my_text(ax, 0.2, 0.7, r'$\rho$ = {:.2f} - {:.2f} '.format(prof4.rho_psi.values[0], prof4.rho_psi.values[19]))
ax.set_xlabel('t [s]')
ax.set_ylabel(r'$v_\perp$ [km/s]')
ax.set_title('#{} ; isweep {}'.format(shot, isweep_list))
plt.tight_layout()



fig ,ax = plot_1d([], [], grid=True)
err3 = np.zeros((2, len(prof3.v_perp.values[20:])))
err3[1] = prof3.dv_perp_up.values[20:]
err3[0] = prof3.dv_perp_low.values[20:]

err4 = np.zeros((2, len(prof4.v_perp.values[20:])))
err4[1] = prof4.dv_perp_up.values[20:]
err4[0] = prof4.dv_perp_low.values[20:]


plt.errorbar(prof3.t.values[20:], prof3.v_perp.values[20:],yerr=err3,linewidth=0.8,color='teal')
# my_text(ax, 0.2, 0.9, r'$\rho$ = {:.2f} $\pm$ {:.2f}'.format(np.mean(prof3.rho_psi.values[:20]), np.std(prof3.rho_psi.values[:20])))
my_text(ax, 0.2, 0.9, r'$\rho$ = {:.2f} - {:.2f} '.format(prof3.rho_psi.values[20], prof3.rho_psi.values[39]))
# prof4 = prof[prof.channel == channelvals[1]]
plt.errorbar(prof4.t.values[20:], prof4.v_perp.values[20:],yerr=err4,linewidth=0.8,color='firebrick')
# my_text(ax, 0.2, 0.6, r'$\rho$ = {:.2f} $\pm$ {:.2f}'.format(np.mean(prof4.rho_psi.values[:20]), np.std(prof4.rho_psi.values[:20])))
my_text(ax, 0.2, 0.7, r'$\rho$ = {:.2f} - {:.2f} '.format(prof4.rho_psi.values[20], prof4.rho_psi.values[39]))
ax.set_xlabel('t [s]')
ax.set_ylabel(r'$v_\perp$ [km/s]')
ax.set_title('#{} ; isweep {}'.format(shot, isweep_list))
plt.tight_layout()
# %% HOPPING PROFILE

machine, shot, isweep_list, xmode, channelvals = 'tcv', 81069, [4], 1, [4]

for i in range(len(isweep_list)):
    args = (machine, shot, isweep_list[i], xmode, channelvals)
    prof = DBS_Profile(*args)

    prof.sort_values(by='rho_psi', ascending=True, inplace=True)

fig ,ax = plot_1d([], [], grid=True)
err = np.zeros((2, len(prof.v_perp.values[:20])))
err[1] = prof.dv_perp_up.values[:20]
err[0] = prof.dv_perp_low.values[:20]
ax.errorbar(prof.rho_psi.values[:20], prof.v_perp.values[:20], err, linewidth=0.8,color='teal')
ax.set_xlabel(r'$\rho_\psi$')
ax.set_ylabel(r'$v_\perp$ [km/s]')
ax.set_title('#{} ; isweep {}'.format(shot, isweep_list))
fig.tight_layout()

# fig ,ax = plot_1d([], [], grid=True)
err = np.zeros((2, len(prof.v_perp.values[20:])))
err[1] = prof.dv_perp_up.values[20:]
err[0] = prof.dv_perp_low.values[20:]
ax.errorbar(prof.rho_psi.values[20:], prof.v_perp.values[20:], err, linewidth=0.8,color='teal')
ax.set_xlabel(r'$\rho_\psi$')
ax.set_ylabel(r'$v_\perp$ [km/s]')
ax.set_title('#{} ; isweep {}'.format(shot, isweep_list))
fig.tight_layout()
# %%
