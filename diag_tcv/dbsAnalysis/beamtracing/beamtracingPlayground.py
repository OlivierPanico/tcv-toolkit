
# %% Plot cutoff on density profile
#libraries and parameters

import numpy as np 
import matplotlib.pyplot as plt


from dataAnalysis.utils.utils import get_closest_ind
from dataAnalysis.utils.plot_utils import plot_1d



import warnings
from DBS.beamtracing.src.tcv.io import retrieve_tcv_eq, load_TS_data
from DBS.beamtracing import Beam3dInterface, TCVDUALV, Equilibrium2d, DensityProf1d
from DBS import definitions

machine = 'tcv'

# shot, twindow = 80935, [1.6,1.7]
# shot, twindow = 80162, [1,1.2]
# shot, twindow = 80839, [0.8,1]
shot, twindow = 80957, [1, 1.2]

modex, anglepol, zshift = 1, -50, 0.0

# FreqGHz = np.linspace(48.5, 75, 6)
# FreqGHz = [51, 56, 58, 61, 65, 70, 72.5, 66, 62, 55, 50, 53, 59, 63, 67, 71, 73, 64, 60, 54]#np.linspace(48, 70, 6)
FreqGHz = [59, 69]

# time to chose for the equilibrium snapshot:
time = np.mean(twindow)


#%% beamtracing part => do not change this


eq_npz_data = retrieve_tcv_eq(shot) # pulls the equilibrium data from the TCV database and stores it to file if not already present
equilibrium = Equilibrium2d.from_shot(machine, shot, time)
equilibrium.plot()

# load the density profile:
TS_data = load_TS_data(shot) # pulls the TS data from the TCV database and stores it to file if not already present

densityprof = DensityProf1d.from_shot(machine,shot,twindow, freeze_ped_fitparam='all') 
densityprof.plot(detailed=True)


launcher  = [TCVDUALV(freqGHz=f, anglepol=anglepol, modex=modex, zshift=zshift) for f in FreqGHz]
mode = 'X' if launcher[0].modex == 1 else 'O'
outpath = definitions.DATA_TMP_DIR / f'beam3d_tests/{shot}_{time}s_{mode}mode_{definitions.USER_SHORTNAME}.mat'
if not outpath.parent.exists():
    outpath.parent.mkdir(parents=True)

interface = Beam3dInterface(shot, launcher, densityprof, equilibrium ,outpath=outpath)

interface.run_beam3d()
outp = interface.fetch_result()


#%% plot ray paths

np.set_printoptions(threshold=20, precision=3)
from DBS.io.utils import run_line_magic
run_line_magic('matplotlib', 'widget')

from DBS.beamtracing.src.visualize import plot_beam
from DBS.io.utils import run_line_magic
import matplotlib.pyplot as plt
    
fig, ax = plt.subplots(figsize=(4, 6))
try:
    import tcv
    tcv.tcvview(shot, time, col='k', ports=True, linewidths=.6, levels_sol=[])
except ModuleNotFoundError:
    warnings.warn('Could not import tcv module, so not plotting the TCV vessel. If this was not intended, make sure you are on lac10 and have added /usr/local/mdsplus/python to your PYTHONPATH.')
    equilibrium.plot(ax=ax)    
   
for ifreq, freq in enumerate(FreqGHz):
    beam = outp.beam[ifreq]
    dif = outp.dif[ifreq]
    beami = outp.beami[ifreq]
    beami.diagnm = '' # 'tcvvx' or 'difdop', but not sure this is actually needed
    integopt = outp.integopt[ifreq]
    # fig, ax = plt.subplots()
    plot_beam(beam, dif, beami, ax=ax, central_ray=True,
                other_rays=False, label=f'{freq:.1f} GHz', lw=1)
    
    # indicate turning point:
    ax.plot(dif.x / 100, dif.z / 100, 'x', color=f'C{ifreq}', markersize=5)

ax.set_aspect('equal')
mode = 'X' if launcher[0].modex == 1 else 'O'
ax.text(0.99, 0.99, r'$\theta$ = {:.1f}$\deg$'.format(anglepol) + '\n' + r'$\Delta z = {:.1f}$cm'.format(zshift*100) + '\n' + f'{mode} mode', transform=ax.transAxes, ha='right', va='top')
ax.legend(loc='lower right')

print('anglepol = {}, zshift = {}:'.format(anglepol, zshift))
print(f'freq = {FreqGHz},\nrho = {outp.rho},\nkperp = {outp.k_perp}\n')









#%% plot the results
def get_closest_ind(L,val):
    ires=0
    diff_min=abs(L[0]-val)
    for i in range(len(L)):
        diff=abs(L[i]-val)
        if diff<diff_min:
            ires = i
            diff_min = diff
    return ires

ne = densityprof.ne
rho_ne = densityprof.rho_psi

rho_beam = outp.rho

rho_freq_ind = []
for i in range(len(rho_beam)):
    rho_freq_ind.append(get_closest_ind(rho_ne, rho_beam[i]))

rhoinit = 0.7
rhofin = 1.2

rhoinit_ind = get_closest_ind(rho_ne, rhoinit)
rhofin_ind = get_closest_ind(rho_ne, rhofin)


# fig, ax = plot_1d([], [], grid=True)
fig, ax = plt.subplots()
ax.plot(rho_ne[rhoinit_ind:rhofin_ind], ne[rhoinit_ind:rhofin_ind], color='black')
ax.scatter(rho_ne[rho_freq_ind], ne[rho_freq_ind], color='red', marker='s', s=40)
ax.axvline(1, color='black', linestyle='--')
for i in range(len(outp.k_perp)):
    ax.text(0.9, 0.95-0.05*i, r'$k_\perp =$ {:.2f}'.format(outp.k_perp[i]), horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes)
ax.grid()
ax.set_title('#{}'.format(shot))
# %%
