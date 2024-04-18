#%%

import matplotlib.pyplot as plt
from DBS.analysis import DBS_Profile
from DBS.beamtracing import DBSbeam


# machine, shot, isweep_list, xmode, channelvals = 'tcv', 80745, [9,11,13], 1, [1,2]

#to do: shot 80322 (ECRH / MIX (actually should be only NBI) / NBI)

machine, shot, isweep_list, xmode, channelvals = 'tcv', 80376, [17], 1, [1]

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
machine, shot, isweep_list, xmode, channelvals = 'tcv', 80745, [17,18], 1, [1]

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
# ax.set_xlim(0.65, 0.9)
ax.set_ylim(-3.2, -0.9)
ax.legend(loc='upper left')

plt.tight_layout()