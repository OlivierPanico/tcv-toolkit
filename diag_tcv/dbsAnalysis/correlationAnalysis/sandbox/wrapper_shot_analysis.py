# %% Wrapper file to perform the whole analysis of a shot at a given time

#General import
import numpy as np
import matplotlib.pyplot as plt

#Local import
from DBS import definitions as defs


shot, time = 78549, 1 #For plasma data
modex, anglepol, isweep, channelval = 1, -50, 2, 1 #For beamtracing
ifreq = 1 #For raw data and spectrums 
machine = 'tcv' #for spectrum fit

### ============================= ###
### RETRIEVE DATA FROM THE PLASMA ###
### ============================= ###
from DBS.beamtracing.src.tcv.io import retrieve_plasma_data
from DBS.beamtracing.src.tcv.show_example import show_density_prof
outpath = defs.DATA_PLASMA_DIR / 'tcv' / f'{shot}_t{time:0.2f}s.mat'
plasma = retrieve_plasma_data(shot, time, verbose=True, outpath=outpath)
#Show the density profile used for beamtracing
show_density_prof(shot,time, plasmaobj=plasma)

#%%
### ==================== ###
### PERFORM BEAM TRACING ###
### ==================== ###
from DBS.beamtracing.src.launcher import DIFDOP, TCVDUALV
from DBS.beamtracing.src.equilibrium import Equilibrium2d
from DBS.beamtracing.src.densityprof import DensityProf1d

from DBS.processing.fDop_estimation import InterfaceFDopEstimation
interf = InterfaceFDopEstimation(shot, isweep, channelval,machine='tcv')

from DBS.beamtracing.interface import Beam3dInterface

FreqGHz = interf.dop.F
launcher  = [TCVDUALV(freqGHz=f, anglepol=anglepol, modex=modex, zshift=0.0) for f in FreqGHz]
# load the density profile:
densityprof = DensityProf1d.from_shot(machine='tcv', shot=shot, time=time)
# load the equilibrium:
equilibrium = Equilibrium2d.from_shot(machine='tcv', shot=shot, time=time)
# neprof1d = DensityProf1d.from_shot(shot=66260, machine='tcv')

# beami,plasma,integopt = Beam3dInterface(shot, launcher, densityprof, equilibrium)._prepare_input()

# beami.to_mat('./src/tcv/tmp/beami.mat')# %%
# plasma.to_mat('./src/tcv/tmp/plasmai.mat')# %%
# integopt.to_mat('./src/tcv/tmp/integopt.mat')# %%

mode = 'X' if launcher[0].modex == 1 else 'O'

outpath = defs.get_path_for('beam3d_output', machine='tcv',
                        shot=shot, mode=mode, xmode=modex,channelval=channelval, isweep=isweep)


if not outpath.parent.exists():
    outpath.parent.mkdir(parents=True)

#Preparing the data for the beam tracing
interface = Beam3dInterface(shot, launcher, densityprof, equilibrium ,outpath=outpath)
#Performing the beamtracing
interface.run_beam3d()
#Fetching results from beamtracing 
outp = interface.fetch_result()
#Plotting data from beamtracing
from DBS.beamtracing.src.visualize import plot_gola, plot_beam
from DBS.io.utils import run_line_magic
# run_line_magic('matplotlib', 'inline')
run_line_magic('matplotlib', 'widget')
fig, ax = plt.subplots()
import tcv
tcv.tcvview(shot, time, col='k', ports=True, linewidths=.6, levels_sol=[])

ax= plt.gca()
for ifreq, freq in enumerate(FreqGHz):
    beam = outp.beam[ifreq]
    dif = outp.dif[ifreq]
    beami = outp.beami[ifreq]
    beami.diagnm = '' # 'tcvvx' or 'difdop', but not sure this is actually needed
    integopt = outp.integopt[ifreq]
    # fig, ax = plt.subplots()
    plot_beam(beam, dif, beami, ax=ax, central_ray=True, other_rays=False)#, color='r')

# plot_gola(beami, ax=ax, color='k', lw=2)
ax.set_aspect('equal')

#%%
### ==================== ###
### RAW SIGNALS FROM DBS ###
### ==================== ### 
from DBS.io.interface import DataInterface_ifreq

data_interface = DataInterface_ifreq(shot, isweep=isweep,ifreq=ifreq, channelval=channelval, machine=machine)
par = data_interface.params

t, x, y = data_interface.get_signal()

plt.figure()
plt.plot(t, x, color='blue')
plt.plot(t,y, color='red')

#%%
### ========================= ###
### FDOP ESTIMATION WITH FITS ###
### ========================= ###
from DBS.processing.sigprocessing import init_specobjs, perform_specobj_fits, show_spec, make_title, get_fDop_from_fit_results

sobjs = init_specobjs(data_interface, ifreqs=[10])
        
for s in sobjs[:]:
    perform_specobj_fits(s)
    fig, ax = plt.subplots()
    plot_dict = show_spec(s, ax=ax)
    make_title(ax, data_interface, s.header.ifreq)
    
ax.set_xlim(-5,5)

print(get_fDop_from_fit_results(sobjs[0].fit_params, sobjs[0].xscale))



# %%