#%% Get the values of g, C for a given discharge in TCV

#General import
import numpy as np
import matplotlib.pyplot as plt

#Local import 
from DBS.beamtracing.src.tcv.io import retrieve_plasma_data #for density and mag eq
from DBS import definitions as defs
from DBS.io.utils import get_closest_ind

#tcv import
import tcv




def calculate_g_C_from_plasma_parameter(nsep, Tsep, Bsep, qsep):
    
    ### PHYSICAL CONSTANTS ###
    k_b = 1.380649E-23   #Boltzmann Constant
    epsilon_0 = 8.854E-12 #Void permittivity
    m_i = 1.6726E-27    #Mass proton
    m_e = 9.1094E-31    #Mass electron
    e = 1.6022E-19        #Charge electron


    ### VARIABLES TCV ###
    R = 0.88
    a = 0.25
    Z = 1

    ### NORMALISATION ###
    n_0 = 1  #10^19
    T_0 = 50 #eV

    cs = np.sqrt(e*Tsep/m_i)  #sound speed in m/s 
    cs_norm = np.sqrt(e*T_0/m_i)
    # print("cs : " + "{:e}".format(cs) + " in m/s")
    omega_cs = (e*Bsep)/(m_i) #frequency in s**-1
    # print("omega_cs : " + "{:e}".format(omega_cs) + " in 1/s")
    rho_s = cs/omega_cs #sound larmor radius in m
    rho_s_norm = cs_norm/omega_cs
    # print("rho_s : " + "{:e}".format(rho_s) + " m")

    inverse_aspect_ratio = a/(R)
    lnLambda = 15 - 0.5*np.log(10*nsep) + np.log(1E-3 * Tsep) #Coulomb logarithm from X.Garbet magic booklet 
    nu_ei = (1/omega_cs) * ((4*np.sqrt(2*np.pi))/3) * ((e**2/(4*np.pi*epsilon_0))**2) * np.mean(lnLambda) * ((nsep*10**19)/(np.sqrt(m_e)*((e*Tsep)**(3/2)))) #electron-ion collision frequency (from X.Garbet magic booklet)
    # print("nu_ei : " + "{:e}".format(nu_ei) + " normalized to omega_cs")²
    
    ### ADIABATICITY PARAMETER ###
    #kparallel: we choose the length of a field line: L = pi*q*R
    kpar = ((2*np.pi*rho_s_norm)/(qsep*np.pi*R))
    sigma = (e*Bsep)/(m_e * nu_ei*omega_cs)
    C = sigma*kpar**2
    # print("C : " + "{:.2e}".format(C))

    ### INTERCHANGE PARAMETER ###
    g = 2*rho_s/R
    # print("g : " + "{:.2e}".format(g)) 

    return g, C

def get_g_C_from_shot(shot, time, rho=0.8):
    #thomson values are 
    t_thomson, rho_thomson, ne, ne_err, Te, Te_err = get_thomson_data(shot=shot)
    #Values from the raytracing code are avged on -100ms/+100ms
    B0, ne_raytracing, rho_psi_ne_raytracing, q_psi, rho_psi = get_ne_B_from_raytracing(shot=shot, time=time)
    
    #We average the density and time profiles on -100 ms / +100 ms
    tinit = time - 0.1
    tfin = time + 0.1
    ind_tinit_thomson = get_closest_ind(t_thomson, tinit)
    ind_tfin_thomson = get_closest_ind(t_thomson, tfin)

    #Sanity check: both rho should be equal
    ind_rho_thomson = get_closest_ind(rho_thomson, rho)
    ind_rho_psi = get_closest_ind(rho_psi, rho)

    Tsep = np.array(np.mean(Te[ind_tinit_thomson:ind_tfin_thomson, ind_rho_thomson]))
    nsep = np.array(np.mean(ne[ind_tinit_thomson:ind_tfin_thomson, ind_rho_thomson]))/1e19
    Bsep = np.array(B0)
    qsep = np.array(q_psi[ind_rho_psi])
    

    g, C = calculate_g_C_from_plasma_parameter(nsep, Tsep, Bsep, qsep)

    return g, C, Tsep, nsep, Bsep, qsep

### ======== ###
### PLOTTING ###
### ======== ###
def plot_g_C_prof(shot, time):
    #thomson values are 
    t_thomson, rho_thomson, ne, ne_err, Te, Te_err = get_thomson_data(shot=shot)
    #Values from the raytracing code are avged on -100ms/+100ms
    B0, ne_raytracing, rho_psi_ne_raytracing, q_psi, rho_psi = get_ne_B_from_raytracing(shot=shot, time=time)
    
    #We average the density and time profiles on -100 ms / +100 ms
    tinit = time - 0.1
    tfin = time + 0.1
    ind_tinit_thomson = get_closest_ind(t_thomson, tinit)
    ind_tfin_thomson = get_closest_ind(t_thomson, tfin)

    Tsep = np.array( np.mean(Te[ind_tinit_thomson:ind_tfin_thomson, :], axis=0))
    nsep = np.array(np.mean(ne[ind_tinit_thomson:ind_tfin_thomson, :], axis=0))/1e19
    Bsep = np.array(B0)
    qsep = np.array(q_psi)

    g, C = calculate_g_C_from_plasma_parameter(nsep, Tsep, Bsep, qsep)

    plt.figure()
    plt.plot(rho_thomson, g)

    plt.figure()
    plt.plot(rho_thomson, C)

    return rho_thomson, g, C




#%% 

# shot_list = [78922, 78923, 78924, 78925, 78926, 78927, 78928, 78929, 78931, 78932, 78933,
#             78934, 78936, 78939, 78942]
# shot_list = [78966, 78967, 78969, 78970]

shot_list = [78936, 78939, 78942, 78966, 78967, 78969, 78970]


shot_list_L_mode = [78922, 78923, 78924, 78926, 78927, 78931, 78932, 78933, 78934]
time_list_L_mode = [1.3, 1.2, 0.3, 0.7, 0.3, 1.4, 0.3, 1.2, 0.3]

comm_L_mode = ['ok', 'ok', 'ohmic', 'ohmic', 'ohmic', 'no good flat top', 'ohmic', 'good L mode',
                'ohmic']

shot_list_H_mode = [78924, 78926, 78927, 78929, 78932, 78934]
time_list_H_mode = [1.3, 1.1, 1.2, 0.9, 1.2, 1.1]
comm_H_mode = ['high ne', 'is it qce?', 'qce?', 'qce', 'qce']

shot_list_intermediate = [78925, 78928, 78930]
time_list_intermediate = [1.2, 0.8, 1.2]
comm_intermediate = ['qce approx L mode', 'no good flat top', 'good qce']

g_L, C_L, Tsep_L, nsep_L, Bsep_L, qsep_L = np.zeros((len(shot_list_L_mode))), np.zeros((len(shot_list_L_mode))), np.zeros((len(shot_list_L_mode))), np.zeros((len(shot_list_L_mode))), np.zeros((len(shot_list_L_mode))), np.zeros((len(shot_list_L_mode)))


for i, shot in enumerate(shot_list_L_mode):
    print(i, shot)
    time = time_list_L_mode[i]
    g_L[i], C_L[i], Tsep_L[i], nsep_L[i], Bsep_L[i], qsep_L[i] = get_g_C_from_shot(shot, time, 0.8)
    print(g_L[i], C_L[i], Tsep_L[i], nsep_L[i], Bsep_L[i], qsep_L[i] )


g_H, C_H, Tsep_H, nsep_H, Bsep_H, qsep_H = np.zeros((len(shot_list_H_mode))), np.zeros((len(shot_list_H_mode))), np.zeros((len(shot_list_H_mode))), np.zeros((len(shot_list_H_mode))), np.zeros((len(shot_list_H_mode))), np.zeros((len(shot_list_H_mode)))

for i, shot in enumerate(shot_list_H_mode):
    print(i, shot)
    time = time_list_H_mode[i]
    g_H[i], C_H[i], Tsep_H[i], nsep_H[i], Bsep_H[i], qsep_H[i] = get_g_C_from_shot(shot, time, 0.8)
    print(g_H[i], C_H[i], Tsep_H[i], nsep_H[i], Bsep_H[i], qsep_H[i] )


#%%
fig, ax = plt.subplots(figsize=(9,5))
ax.scatter(C_L, g_L, color='blue', label='L mode')
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.grid()
# plt.xlim(1.2,2.5)
# plt.ylim(0.0024,0.0032)
ax.set_xlabel(r'$C$')
ax.set_ylabel(r'$g$')
for i, shot in enumerate(shot_list_L_mode):
    ax.text(C_L[i], g_L[i], str(shot))

ax.scatter(C_H, g_H, color='red', label = 'H mode')

for i, shot in enumerate(shot_list_H_mode):
    ax.text(C_H[i], g_H[i], str(shot))

ax.legend()

# plt.grid()


plt.subplots(sharex=True)
plt.subplot(2,2,1)
plt.scatter(C_L, Tsep_L, color='blue')
plt.scatter(C_H, Tsep_H, color='red')

# plt.scatter(g_L/C_L, Tsep_L, color='blue')
# plt.scatter(g_H/C_H, Tsep_H, color='red')
plt.ylabel(r'$T$ [$eV$]')
# plt.xscale('log')
plt.grid()
plt.xlabel(r'C')

plt.subplot(2,2,2)
plt.scatter(C_L, nsep_L, color='blue')
plt.scatter(C_H, nsep_H, color='red')

# plt.scatter(g_L/C_L, nsep_L, color='blue')
# plt.scatter(g_H/C_H, nsep_H, color='red')
plt.ylabel(r'$n_e$ [$10^{19}$ $m^{-3}$]')
# plt.xscale('log')
plt.grid()
plt.xlabel(r'C')

plt.subplot(2,2,3)
plt.scatter(C_L, Bsep_L, color='blue')
plt.scatter(C_H, Bsep_H, color='red')

# plt.scatter(g_L/C_L, Bsep_L, color='blue')
# plt.scatter(g_H/C_H, Bsep_H, color='red')
plt.ylabel(r'$B_{axis}$')
# plt.xscale('log')
plt.grid()
plt.xlabel(r'C')

plt.subplot(2,2,4)
plt.scatter(C_L, qsep_L, color='blue')
plt.scatter(C_H, qsep_H, color='red')

# plt.scatter(g_L/C_L, qsep_L, color='blue')
# plt.scatter(g_H/C_H, qsep_H, color='red')
plt.ylabel(r'$q$')
# plt.xscale('log')
plt.grid()
plt.xlabel(r'C')

plt.tight_layout()



#%%
# shot = 78549
shot = 78922 
time = 1.2

with tcv.shot(shot) as dis:

    # NOTE: Not all data is available from any Lac#, connect to standard LAC to be shure
    #data = MDSplus.Data.execute(query)

    ##We load the fitted data from Thomson scattering diagnostic
    #electron density shape(t, rho)
    key=r'\tcv_shot::top.results.thomson.profiles.auto:ne'
    ne = dis.tdi(key)

    #electron temperature shape(t,rho)
    key = r'\tcv_shot::top.results.thomson.profiles.auto:te'
    Te = dis.tdi(key)

    #time
    key=r'\tcv_shot::top.results.thomson.profiles.auto:time'
    t_thomson = dis.tdi(key)

    #rho 
    key = r'\tcv_shot::top.results.thomson.profiles.auto:rho'
    rho_thomson = dis.tdi(key)

    #The magnetic seems to not be accessible in python
    # key = r'\magnetics::rbphi' 
    # rbphi = dis.tdi(tcv_eq("bzero","fbte"))

    key = r'\results::q_95'
    q95 = dis.tdi(key)
  
    key = r'\results::r_contour'
    r_contour = dis.tdi(key)

    key = r'\results::r_axis'
    print(dis.tdi(key))

    key = r'\results::z_axis'
    print(dis.tdi(key))

#%%
ttime = 40
xinit =  30
xfin = 40

plt.figure()
plt.plot(rho_thomson[xinit:xfin], ne[ttime,xinit:xfin]/(1E19))

plt.figure()
plt.plot(rho_thomson[xinit:xfin],Te[ttime,xinit:xfin])

plt.figure()
plt.plot(rho_thomson[xinit:xfin], ne[ttime,xinit:xfin]/(1E19*Te[ttime,xinit:xfin]**3/2))

#%% Value from matlab code (using wrapper)
#Local import 
from DBS.beamtracing.src.tcv.io import retrieve_plasma_data #for density and mag eq
from DBS import definitions as defs


outpath = defs.DATA_PLASMA_DIR / 'tcv' / f'{shot}_t{time:0.2f}s.mat'

plasma = retrieve_plasma_data(shot=shot, time=time, verbose=True, outpath=None)

#%%
Btor_tot = plasma.BtotStruct.B.phi
R_B = plasma.BtotStruct.coord.R  # 1d array

nr, nz = np.shape(Btor_tot)
B0 = abs(np.mean(Btor_tot[nr//2, :]))

#position of LCFS: r_contour
rsep = np.mean(r_contour[0]) #For now an avg on the whole shot but should be taken at the right time

ind_sep = get_closest_ind(R_B, rsep)

Btorsep = abs(np.mean(Btor_tot[ind_sep, :]))

ne_raytracing = plasma.neStruct.ne
rho_psi_ne_raytracing = plasma.neStruct.rho_psi



# %%

import matplotlib.pyplot as plt


ttime = 20

# fig, ax = plt.subplots()
# ax.plot(t_thomson, ne)
fig, ax = plt.subplots()
ax.plot(t_thomson, Te)


fig, ax = plt.subplots()
ax.plot(rho_thomson, ne[ttime,:], marker='+', color='blue')
fig, ax = plt.subplots()
ax.plot(rho_thomson, Te[ttime,:], marker='o', color='red')







#%% Estimation of g and C


### PHYSICAL CONSTANTS ###
k_b = 1.380649E-23   #Boltzmann Constant
epsilon_0 = 8.854E-12 #Void permittivity
m_i = 1.6726E-27    #Mass proton
m_e = 9.1094E-31    #Mass electron
e = 1.6022E-19        #Charge electron


### VARIABLES TCV ###
R = 0.88
a = 0.25
Z = 1



ind_rho_thomson = get_closest_ind(rho_thomson, 0.9)
ind_t_thomson = get_closest_ind(t_thomson, time)

Tsep = np.array(Te[ind_t_thomson, ind_rho_thomson])
nsep = np.array(ne[ind_t_thomson, ind_rho_thomson])/1e19
Bsep = np.array(Btorsep)
qsep = 10#np.array(np.mean(q95, axis=0))

print(Tsep, nsep, Bsep)

# ind_rho_thomson = np.where(rho_thomson == 0.9)[0][0]
# ind_t_thomson = np.where(t_thomson == time)[0][0]

print("\n ### Normalization parameters ### ")
cs = np.sqrt(e*Tsep/m_i)  #sound speed in m/s 
# print("cs : " + "{:e}".format(cs) + " in m/s")
omega_cs = (e*Bsep)/(m_i) #frequency in s**-1
# print("omega_cs : " + "{:e}".format(omega_cs) + " in 1/s")
rho_s = cs/omega_cs #sound larmor radius in m
print("rho_s : " + "{:e}".format(rho_s) + " m")



inverse_aspect_ratio = a/(R)
print ("\n ### Collision frequencies ###")
lnLambda = 15 - 0.5*np.log(10*nsep) + np.log(1E-3 * Tsep) #Coulomb logarithm from X.Garbet magic booklet 
nu_ei = (1/omega_cs) * ((4*np.sqrt(2*np.pi))/3) * ((e**2/(4*np.pi*epsilon_0))**2) * lnLambda * ((nsep*10**19)/(np.sqrt(m_e)*((e*Tsep)**(3/2)))) #electron-ion collision frequency (from X.Garbet magic booklet)
print("nu_ei : " + "{:e}".format(nu_ei) + " normalized to omega_cs")

print("\n ### Adiabaticity and Interchange parameter ###")
### ADIABATICITY PARAMETER ###
#kparallel: we choose the length of a field line: L = pi*q*R
kpar = (2*np.pi*rho_s)/(qsep*np.pi*R)
sigma = (e*Bsep)/(m_e * nu_ei*omega_cs)
C = sigma*kpar**2
print("C : " + "{:.2e}".format(C))

### INTERCHANGE PARAMETER ###
g = 2*rho_s/R
print("g : " + "{:.2e}".format(g)) 


print(g/C)









# %%