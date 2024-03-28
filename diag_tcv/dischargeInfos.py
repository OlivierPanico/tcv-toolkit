#%%
#=== general info ===#
#Author: Olivier Panico
#contact: olivier.panico@free.fr

#Goal: provides quick way of knowing the general information of
#       the studied shot. Gives infos on gradients / DBS channels frequency / positions 


#=== imports ===#
#General imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#Local imports

#DBS
# from DBS.beamtracing.src.tcv.io import retrieve_plasma_data #for density and mag eq
from DBS import definitions as defs

#dataAnalysis
from dataAnalysis.utils.plot_utils import plot_1d, plot_2d, get_cmap_list
from dataAnalysis.utils.utils import get_closest_ind
#TCV
import tcv
import MDSplus as mds


### ================= ###
### get data from tcv ###
### ================= ###
def get_thomson_data(shot):
    
    dthomson = {}
    
    with tcv.shot(shot) as dis:

        # NOTE: Not all data is available from any Lac#, connect to standard LAC to be shure
        # data = MDSplus.Data.execute(query)

        ##We load the fitted data from Thomson scattering diagnostic
        #electron density shape(t, rho)
        key=r'\tcv_shot::top.results.thomson.profiles.auto:ne'
        ne = dis.tdi(key)
        key=r'\tcv_shot::top.results.thomson.profiles.auto:ne:error_bar'
        ne_err = dis.tdi(key)

        #electron temperature shape(t,rho)
        key = r'\tcv_shot::top.results.thomson.profiles.auto:te'
        Te = dis.tdi(key)
        key=r'\tcv_shot::top.results.thomson.profiles.auto:te:error_bar'
        Te_err = dis.tdi(key)

        #time
        key=r'\tcv_shot::top.results.thomson.profiles.auto:time'
        t = dis.tdi(key)

        #rho 
        key = r'\tcv_shot::top.results.thomson.profiles.auto:rho'
        rho = dis.tdi(key)

    dthomson['t'] = t
    dthomson['rho'] = rho
    dthomson['ne'] = ne
    dthomson['ne_err'] = ne_err
    dthomson['Te'] = Te
    dthomson['Te_err'] = Te_err        
    return dthomson


def get_NBI_data(shot):
    
    dNBI = {}
    
    with tcv.shot(shot) as dis:
        
        try:
            key=r'\tcv_shot::top.results.NB1.POWR_TCV'
            NB1 = dis.tdi(key)
            dNBI['NB1'] = NB1
            dNBI['tagNB1'] = True
        except:
            dNBI['tagNB1'] = False
            print('No data for NB1')
        try:
            key=r'\tcv_shot::top.results.NB2.POWR_TCV'
            NB2 = dis.tdi(key)
            dNBI['NB2'] = NB2
            dNBI['tagNB2'] = True
        except:
            dNBI['tagNB2'] = False
            print('No data for NB2')
        
        try:
            key=r'\ATLAS::NBH.DATA.MAIN_ADC:DATA'
            t_NBI = dis.tdi(key).dim_0.values
            dNBI['t'] = t_NBI
        except:
            print('No tcv time')

    return dNBI


def get_ECRH_data(shot):
    '''
    ECRH data structure:
        12 colums: 1 to 6 => power from lines 1 to 6 (X2 launchers usually)
                   7 to 9 => power from lines 7 to 9 (X3 launchers usually)
                   10 to 11 => power from launchers 10 / 11
                   last column => total power
    to access the total power :
        ecrh_tot = ecrh.values[:,11]
    to access time:
        ecrh_time = ecrh.dim_0.values
    '''
    
    dECRH = {}
    
    with tcv.shot(shot) as dis:
        try:
            key=r'\RESULTS::TORAY.INPUT:P_GYRO'
            ecrh = dis.tdi(key)
            dECRH['ECRH'] = ecrh
            dECRH['tagECRH'] = True
        except:
            dECRH['tagECRH'] = False
            print('No ECRH data')
            
        
    return dECRH



def get_cxrs_data(shot, set_Nan_zero=True):
    
    a=mds.Tree('tcv_shot', shot)
    
    dcxrs={}
    
    try:
        tiNode = a.getNode('\RESULTS::CXRS.PROFFIT:TI')
        rho_ti = tiNode.getDimensionAt(0).data()[0,:]
        time_ti = tiNode.getDimensionAt(1).data()
        ti = tiNode.data()
        ti_err = a.getNode('\RESULTS::CXRS.PROFFIT:TI:ERR').data()
        dcxrs['tagCXRSti'] = True
    except:
        dcxrs['tagCXRSti'] = False
        print('No CXRS Ti data')
        return dcxrs
    
    if set_Nan_zero:
        tidf = pd.DataFrame(ti)
        tidf.fillna(0, inplace=True)
        ti = tidf.to_numpy()

        ti_errdf = pd.DataFrame(ti_err)
        ti_errdf.fillna(0, inplace=True)
        ti_err = ti_errdf.to_numpy()

    if dcxrs['tagCXRSti']:
        dcxrs['ti'] = ti 
        dcxrs['ti_err'] = ti_err
        dcxrs['rho'] = rho_ti
        dcxrs['t'] = time_ti
        
    return dcxrs

def plot_thomson_profile(shot, time, rhomin=0, rhomax=1):
    dthomson = get_thomson_data(shot=shot)
    
    t=dthomson['t']
    rho=dthomson['rho']
    ne=dthomson['ne']
    ne_err=dthomson['ne_err']
    Te=dthomson['Te']
    Te_err=dthomson['Te_err']       
    
    tinit = time - 0.1
    tfin = time + 0.1
    ind_tinit_thomson = get_closest_ind(t, tinit)
    ind_tfin_thomson = get_closest_ind(t, tfin)

    ind_rhomin = get_closest_ind(rho, rhomin)
    ind_rho_max = get_closest_ind(rho, rhomax)+1
    print(ind_rho_max)

    Te = np.array(np.mean(Te[ind_tinit_thomson:ind_tfin_thomson, ind_rhomin:ind_rho_max], axis=0))
    ne = np.array(np.mean(ne[ind_tinit_thomson:ind_tfin_thomson, ind_rhomin:ind_rho_max], axis=0))
    
    Te_err = np.array( np.mean(Te_err[ind_tinit_thomson:ind_tfin_thomson, ind_rhomin:ind_rho_max], axis=0))
    ne_err = np.array(np.mean(ne_err[ind_tinit_thomson:ind_tfin_thomson, ind_rhomin:ind_rho_max], axis=0))

    rho = rho[ind_rhomin:ind_rho_max]

    fig, ax = plot_1d([], [], grid=True)
    ax.plot(rho, Te, color='red', marker='+', label=r'T_e')
    ax.fill_between(rho, Te-Te_err, Te+Te_err, color='red', alpha=0.2)
    ax.plot([],[], color='blue', marker='x', label = r'n_e')
    # ax.grid()
    ax.legend()
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$T_e$ $[eV]$', color='red')
    # fig, ax2 = plt.subplots()
    ax2 = ax.twinx()
    ax2.plot(rho, ne, color='blue', marker='o')
    ax2.fill_between(rho, ne-ne_err, ne+ne_err, color='blue', alpha=0.2)
    ax2.set_ylabel(r' $n_e$ $[10^{19} m^{-3}]$', color='blue')
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(19,19))
    # ax2.grid()
    

def plot_Te_time(shot, rho_list, tmin=0, tmax=-1):
    dthomson = get_thomson_data(shot=shot)
    
    t=dthomson['t']
    rho=dthomson['rho']
    ne=dthomson['ne']
    ne_err=dthomson['ne_err']
    Te=dthomson['Te']
    Te_err=dthomson['Te_err']       
    
    nt = len(t)
    nx = len(rho)
    
    if rho_list[0]=='all':
        rho_list=rho.values
    
    Te_to_plot = np.zeros((nt, len(rho_list)))
    Te_err_to_plot = np.zeros((nt, len(rho_list)))
    
    clist = get_cmap_list(len(rho_list), cmap_name='viridis')
    
    fig, ax = plot_1d([], [], grid=True, xlabel='time [s]', ylabel=r'$T_e$ [$eV$]')
    
    for i, rho_loc in enumerate(rho_list):
        rho_loc_ind = get_closest_ind(rho, rho_loc)
        Te_to_plot[:,i] = Te[:,rho_loc_ind]
        Te_err_to_plot[:,i] = Te_err[:,rho_loc_ind]
        
        
        ax.plot(t, Te[:, rho_loc_ind], color=clist[i], marker='+', label=r'$\rho = {:.3f}$'.format(rho[rho_loc_ind]))
        ax.fill_between(t, Te[:,rho_loc_ind]-Te_err[:, rho_loc_ind], Te[:,rho_loc_ind]+Te_err[:,rho_loc_ind], color=clist[i], alpha=0.2)

    if len(rho_list)>10:
        print('legend too large')
    else:
        plt.legend(fontsize='14')
    
    plt.title(r'#{}'.format(shot))
    plt.tight_layout()


def plot_NBI(shot):
    
    dNBI = get_NBI_data(shot)
    
    fig, ax = plot_1d([], [], grid=True, xlabel='time [s]', ylabel='NBI')
    
    tNBI = dNBI['t']
    if dNBI['tagNB1']==True:
        NB1 = dNBI['NB1']
        ax.plot(tNBI, NB1*1e3, color='red', marker='+', label='NB1')
    elif dNBI['tagNB2']==True:
        NB2 = dNBI['NB2']
        ax.plot(tNBI, NB2*1e3, color='blue', marker='x', label='NB2')
    else:
        print('No NBI data')
    
    plt.legend(fontsize='14')
    plt.title(r'#{}'.format(shot))
    plt.tight_layout()
    

def plot_ECRH(shot):
    
    dECRH = get_ECRH_data(shot)
    
    if dECRH['tagECRH'] is True:
        ecrh = dECRH['ECRH']
        ecrh_tot = ecrh.values[:,11]
        ecrh_time = ecrh.dim_0.values
        fig, ax = plot_1d(ecrh_tot, ecrh_time, color='black', marker='+', grid=True, 
                          xlabel='time [s]', ylabel='ECRH [kW]')
        plt.title(r'#{}'.format(shot))
        plt.tight_layout()
        
    else:
        print('No ECRH data')
        
def plot_Ti_time(shot, rho_list, remove_zeros=True):
    
    dcxrs = get_cxrs_data(shot)
    if dcxrs['tagCXRSti'] is False:
        print('No CXRS Ti data to plot')
        return
    else:
        t=dcxrs['t']
        rho=dcxrs['rho']
        ti=dcxrs['ti']
        ti_err=dcxrs['ti_err']
        
        nt=len(t)
        nx=len(rho)
        
   
        if rho_list[0]=='all':
            rho_list=rho.values
        
        ti_to_plot = np.zeros((nt, len(rho_list)))
        ti_err_to_plot = np.zeros((nt, len(rho_list)))
        
        clist = get_cmap_list(len(rho_list), cmap_name='viridis')
    
        fig, ax = plot_1d([], [], grid=True, xlabel='time [s]', ylabel=r'$T_i$ [$eV$]')
        for i, rho_loc in enumerate(rho_list):
            
            #print(Ti.shape)
            #nonzero_ind = np.nonzero(Ti)
            
            rho_loc_ind = get_closest_ind(rho, rho_loc)
            
            ti_to_plot[:,i] = ti[:,rho_loc_ind]
            ti_err_to_plot[:,i] = ti_err[:,rho_loc_ind]
        
        
            #ax.scatter(t, ti[:, rho_loc_ind], color=clist[i], marker='+', label=r'$\rho = {}$'.format(rho_loc))
            #ax.fill_between(t, ti[:,rho_loc_ind]-ti_err[:, rho_loc_ind], ti[:,rho_loc_ind]+ti_err[:,rho_loc_ind], color=clist[i], alpha=0.2)
            ax.errorbar(t, ti[:,rho_loc_ind], ti_err[:,rho_loc_ind], color=clist[i],linestyle='', marker='+', label=r'$\rho = {:.3f}$'.format(rho[rho_loc_ind]))


        if len(rho_list)>10:
            print('legend too large')
        else:
            plt.legend(fontsize='14')
        
        if remove_zeros:
            plt.ylim(np.mean(ti_to_plot), np.max(ti_to_plot)+np.max(ti_err_to_plot))
        
        plt.title(r'#{}'.format(shot))
        plt.tight_layout()


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


def deprecated_get_ne_B_from_raytracing(shot,time):

    # plasma = retrieve_plasma_data(shot=shot, time=time, verbose=False, outpath=None)
    
    Btor_tot = plasma.BtotStruct.B.phi
    R_B = plasma.BtotStruct.coord.R  # 1d array
    nr, nz = np.shape(Btor_tot)
    B0 = abs(np.mean(Btor_tot[nr//2, :]))

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
