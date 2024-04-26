#%%
### Reference shot for mix ECRH / NBI
# a=TCVShot(80257)


### Mix ECRH / NBI 
## USN
# a=TCVShot(80322) 
# a=TCVShot(80338)
## LSN
#a=TCVShot(80335) #Disrupted at 0.4s
# a=TCVShot(80336)  #Data available but disrupts at 1.2s


### Mix NB1 / NB2
# a=TCVShot(80324) #good 
# a=TCVShot(80328) #good

### Density ramp
# a=TCVShot(80376) #Some nans but overall ok ?

a.get_cxrs_fit(set_nan_to_zero=False)
a.get_thomson_fit()

rho_plot = 0.9

thrho09 = get_closest_ind(a.th_rho, rho_plot)
cxrsrho09 = get_closest_ind(a.cxrs_rho, rho_plot)

fig, ax = plot_1d(a.th_te[thrho09,:], a.th_time,marker='+', color='blue', label=r'$T_e ; \rho = ${:.3f}'.format(a.th_rho[thrho09]), grid=True)
ax.fill_between(a.th_time, a.th_te[thrho09,:]- a.th_te_err[thrho09,:], a.th_te[thrho09,:]+a.th_te_err[thrho09,:], color='blue', alpha=0.2)
ax.plot(a.cxrs_time[:], a.cxrs_ti[:,cxrsrho09],marker='x', color='red', label=r'$T_i ; \rho = ${:.3f}'.format(a.cxrs_rho[cxrsrho09]))
ax.fill_between(a.cxrs_time[:], a.cxrs_ti[:,cxrsrho09]-a.cxrs_ti_err[:,cxrsrho09], a.cxrs_ti[:,cxrsrho09]+a.cxrs_ti_err[:,cxrsrho09], alpha=0.2, color='red')
ax.set_title('#{}'.format(a.shot))
ax.legend()
ax.set_xlabel('time [s]')
ax.set_ylabel(r'$T_i, T_e$')
# ax.axvline(0.9, color='green')

# ax.axvline(1.1, color='black')
# ax.axvline(1.3, color='black')
# ax.axvline(1.5, color='green')
# ax.axvline(1.7, color='black')

#%%


'''
Has been replaced by results from  raytracing because could not get  array for mag eq 
def get_mag_eq_data(shot):
    
    with tcv.shot(shot) as dis:
        #The magnetic seems to not be accessible in python
        # key = r'\magnetics::rbphi' 
        # rbphi = dis.tdi(tcv_eq("bzero","fbte"))

        key = r'\results::q_95'
        q95 = dis.tdi(key)
    
        key = r'\results::r_contour'
        r_contour = dis.tdi(key)

    return q95, r_contour
'''

from DBS.beamtracing.src.tcv.io import _retrieve_plasma_data
def deprecated_get_ne_B_from_raytracing(shot,time):

    plasma = _retrieve_plasma_data(shot=shot, time=time, verbose=False, outpath=None)
    
    Btor_tot = plasma.BtotStruct.B.phi
    
    R_B = plasma.BtotStruct.coord.R  # 1d array
    nr, nz = np.shape(Btor_tot)
    B0 = abs(np.mean(Btor_tot[nr//2, :]))

    return Btor_tot
    
    #position of LCFS: r_contour
    # rsep = np.mean(r_contour[0]) #For now an avg on the whole shot but should be taken at the right time
    # ind_sep = get_closest_ind(R_B, rsep)
    # Btorsep = abs(np.mean(Btor_tot[ind_sep, :]))

    ne = plasma.neStruct.ne
    rho_psi_ne = plasma.neStruct.rho_psi


    q_psi = plasma.eq.q_psi
    # psi_values = plasma.eq.psi_values
    # psi_axis = plasma.eq.psi_axis
    psi_LCFS = 0
    psi_normalized = (plasma.eq.psi_values - plasma.eq.psi_axis) / (psi_LCFS - plasma.eq.psi_axis)
    rho_psi = np.sqrt(abs(psi_normalized))


    return B0, ne, rho_psi_ne, q_psi, rho_psi



#%% ANALYZED DISCHARGES IN HYDROGEN

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