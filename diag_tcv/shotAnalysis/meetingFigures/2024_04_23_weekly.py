
import matplotlib.pyplot as plt
import numpy as np

from dataAnalysis.utils.plot_utils import plot_1d, prep_multiple_subplots

from diag_tcv.shotAnalysis.dischargeInfoMdsObject import TCVShot, plot_profiles_comparison


# %% CW16 resumé slides
# a=TCVShot(80946)
# plt.rcParams.update({'font.size': 16, 'figure.figsize': [8, 5],})
# a.get_FIR()
# a.get_nbi()
# fig, ax = plot_1d([],[], grid=True)
# ax.plot(a.fir_time, a.fir_int_ne/10**19, color='blue', label='nel')
# ax.plot(a.ref_int_ne_time, a.ref_int_ne, color='black', label='ref')
# ax.plot([], [], color='red', label='NBI-1')
# ax.set_xlim(-0.2, 2.5)
# ax2 = ax.twinx()
# ax2.plot(a.nbi_time, a.nb1, color='red')
# ax.legend(loc='upper left')
# ax.set_ylabel(r'nel $[10^{19}]$')
# ax2.set_ylabel(r'NBI-1 $[MW]$')
# ax.set_xlabel(r'T $[s]$')
# ax.set_title('#80946')


plot_profiles_comparison([80322, 80931], [1.55, 0.8])

b=TCVShot(80322)
fig,axs = prep_multiple_subplots(2,1, axgrid=[0,1], sharex=True)
b.plot_heating(ax=axs[0])
axs[0].axvline(1.75, color='black')
a=TCVShot(80931)
a.plot_heating(ax=axs[1])
axs[1].axvline(0.8, color='black')



# plot_profiles_comparison([80257, 80931], [1.65, 0.8])

# b=TCVShot(80257)
# fig,axs = prep_multiple_subplots(2,1, axgrid=[0,1], sharex=True)
# b.plot_heating(ax=axs[0])
# axs[0].axvline(1.65, color='black')
# a=TCVShot(80931)
# a.plot_heating(ax=axs[1])
# axs[1].axvline(0.8, color='black')


# plt.xlim(0.5, 1.2)
# %%




#%% Figures meeting 2024 / 04 / 22


# fig, ax=plot_1d([],[])
# a=TCVShot(80257)
# a.get_collisionality_prof([1.55, 1.65], ax=ax, rhomin=0.3)
# b=TCVShot(80931)
# b.get_collisionality_prof([0.7, 0.8], ax=ax, rhomin=0.3)

fig, ax = plot_1d([],[], grid=True)
ax.plot(a.nu_star_rho, a.nu_star,marker='+', color='blue', label='#{}'.format(a.shot))
plt.axhline(a.zeff/4, color='blue', label=r'$z_{eff}/4$')
ax.plot(b.nu_star_rho, b.nu_star,marker='+', color='red',  label='#{}'.format(b.shot))
plt.axhline(b.zeff/4, color='red', label=r'$z_{eff}/4$')
plt.legend()
ax.set_ylabel(r'$\nu_*$')
ax.set_xlabel(r'$\rho$')
plt.yscale('log')


#%% Profile comparison 

# a = TCVShot(80930)
# a.plot_summary()
# plot_profiles_comparison([80930, 80930, 80930], [0.8, 1.2, 1.8])

# a=TCVShot(80931)
# a.plot_summary()
# plot_profiles_comparison([80931, 80931, 80931], [0.8, 1.2, 1.8])


# a=TCVShot(80935)
# a.plot_summary()
# plot_profiles_comparison([80935, 80935, 80935], [0.8, 1.2, 1.8])

# a=TCVShot(80939)
# a.plot_summary()
# plot_profiles_comparison([80939, 80939, 80939], [0.8, 1.2, 1.8])


### Same pattern repeated 3 times

# a=TCVShot(80940)
# a.plot_summary()
# plot_profiles_comparison([80940, 80940, 80940], [0.8, 1.2, 1.8])

# a=TCVShot(80942)
# a.plot_summary()
# plot_profiles_comparison([80942, 80942, 80942], [0.8, 1.2, 1.8])

# a=TCVShot(80945)
# a.plot_summary()
# plot_profiles_comparison([80945, 80945, 80945], [0.8, 1.2, 1.8])

# a=TCVShot(80946)
# a.plot_summary()
# plot_profiles_comparison([80946, 80946, 80946], [0.8, 1.2, 1.8])


# a=TCVShot(80947)
# a.plot_summary()
# plot_profiles_comparison([80947, 80947], [1, 1.4])

# a=TCVShot(80949)
# a.plot_summary()
# plot_profiles_comparison([80949, 80949], [1, 1.4])


# a=TCVShot(80951)
# a.plot_summary()
# plot_profiles_comparison([80951, 80951], [1, 1.4])

# %%
# plot_profiles_comparison([80322, 80946], [1, 0.8])
# plot_profiles_comparison([80257, 80946], [1.1, 0.8])


a=TCVShot(80322)
fig,axs = prep_multiple_subplots(2,1, axgrid=[0,1], sharex=True)
a.plot_heating(ax=axs[0])
axs[0].axvline(1, color='black')
b=TCVShot(80946)
b.plot_heating(ax=axs[1])
axs[1].axvline(0.8, color='black')


fig, ax=plot_1d([],[])
a.get_collisionality_prof([0.9, 1.05], ax=ax, rhomin=0.3)
b.get_collisionality_prof([0.65, 0.85], ax=ax, rhomin=0.3)

fig, ax = plot_1d([],[], grid=True)
ax.plot(a.nu_star_rho, a.nu_star,marker='+', color='blue', label='#{}'.format(a.shot))
plt.axhline(a.zeff/4, color='blue', label=r'$z_{eff}/4$')
ax.plot(b.nu_star_rho, b.nu_star,marker='+', color='red',  label='#{}'.format(b.shot))
plt.axhline(b.zeff/4, color='red', label=r'$z_{eff}/4$')
plt.legend()
ax.set_ylabel(r'$\nu_*$')
ax.set_xlabel(r'$\rho$')
plt.yscale('log')


#%% Match at high power

plot_profiles_comparison([80163, 80946], [1, 0.8])