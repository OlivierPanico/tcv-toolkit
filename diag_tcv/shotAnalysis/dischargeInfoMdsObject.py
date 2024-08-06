#%%
#=== general info ===#
#Author: Olivier Panico
#contact: olivier.panico@free.fr

#Goal: provides quick way of knowing the general information of
#       the studied shot. 


#=== imports ===#
#General imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#Local imports

#DBS
# from DBS.beamtracing.src.tcv.io import retrieve_plasma_data #for density and mag eq
# from DBS import definitions as defs
from DBS.beamtracing.src.equilibrium import Equilibrium2d


#dataAnalysis
from dataAnalysis.utils.plot_utils import prep_multiple_subplots, plot_1d, plot_2d, get_cmap_list
from dataAnalysis.utils.utils import get_closest_ind

#TCV
import tcv
import MDSplus as mds

### ========== ###
### PARAMETERS ###
### ========== ###
# cmap=plt.cm.viridis #'viridis'
# cmap=plt.cm.jet #'jet'
# cmap = plt.cm.plasma

from matplotlib.colors import LinearSegmentedColormap
cdict = {'red':   [[0.0,  0.0, 0.0],
                   [0.5,  1.0, 1.0],
                   [1.0,  1.0, 1.0]],
         'green': [[0.0,  0.0, 0.0],
                   [0.25, 0.0, 0.0],
                   [0.75, 1.0, 1.0],
                   [1.0,  1.0, 1.0]],
         'blue':  [[0.0,  0.0, 0.0],
                   [0.5,  0.0, 0.0],
                   [1.0,  1.0, 1.0]]}

# cmap = LinearSegmentedColormap('custom', segmentdata=cdict, N=256)

color_list = ['xkcd:blue', 'xkcd:red', 'xkcd:forest green', 'xkcd:teal', 'xkcd:orange', 'xkcd:magenta', 'brown', 'pink', 'grey', 'black']
from matplotlib.colors import ListedColormap
cmap = ListedColormap(color_list)

marker_list = ['+', 'x', 's', 'o', 'v', '^', '<', '>', 'p', 'P', '*', 'h', 'H', 'X', 'd', '|', '_']


### INTERACTIVE FIGURES ###
# from DBS.io.utils import run_line_magic
# run_line_magic('matplotlib', 'widget')
# run_line_magic('matplotlib', 'inline')

### PHYSICAL CONSTANTS ###
k_b = 1.380649E-23   #Boltzmann Constant
epsilon_0 = 8.854E-12 #Void permittivity
m_proton = 1.6726E-27    #Mass proton
m_e = 9.1094E-31    #Mass electron
e = 1.6022E-19        #Charge electron

### ============== ###
### Specifying gas ###
### ============== ###
def isHydrogen(shot):
    hydrogen={}
    hydrogen['80931'] = True
    hydrogen['80939'] = True
    hydrogen['80940'] = True
    hydrogen['80942'] = True
    hydrogen['80946'] = True
    hydrogen['80947'] = True
    hydrogen['80949'] = True
    hydrogen['80951'] = True
    
    if str(shot) not in hydrogen:
        print(' --- hydrogen not filled: putting False as default --- ')
        return False
    return hydrogen[str(shot)]




### ================= ###
### Utility functions ###
### ================= ###
def handle_nan(data, ax0_array, value=None):
    '''
    Use it if full lines of the data are filled with Nans
    Handles nan if you know whether they are located in the first or second axis
    '''
    datadf = pd.DataFrame(data)
    mask = datadf.isna() #mask True if value is Nan
    
    #mask_ax1: returns True when value is NOT nan and False when value is Nan
    #I don't know why we have to put axis=1 
    #then you can take: ax0_array[mask_ax0] and you will have the good data
    mask_ax0  = ~datadf.isna().any(axis=1)
    
    if value is None:
        print('Dropping Nans values')
        #Cleaning the array: 
        datadf_clean = datadf.dropna()
        data_np = datadf_clean.to_numpy()
        return data_np, ax0_array[mask_ax0]
    else:
        print('Replacing Nans with {}'.format(value))
        datadf_replace = datadf.fillna(value, inplace=True)
        data_np = datadf_clean.to_numpy()
        return data_np, ax0_array
   
### General plotting functions

# def plot_(array, t, rho, time_list, array_err = None, ax=None, rhomin=0, rhomax=1, tavg=0.2):
def plot_time(array, t, rho, rho_list, array_err=None, ax=None, tmin=0, tmax=-1, rhoavg=0.2, legend=True, **kwargs):
    '''
    TO PLOT IF THE ARRAY IS (RHO, TIME)
    Otherwise call the same function with np.transpose(array) & np.transpose(array_err)
    '''

    
    clist = get_cmap_list(len(rho_list), cmap)
    ind_tmin = get_closest_ind(t, tmin)
    if tmax==-1:
        ind_tmax=-1
    else:
        ind_tmax = get_closest_ind(t, tmax)+1       
             
    t = t[ind_tmin:ind_tmax]
    
    if ax is None:
        _, ax = plot_1d([], [], grid=True)
    
    for i,rholoc in enumerate(rho_list):
        rhomin = rholoc - rhoavg/2
        rhomax = rholoc + rhoavg/2
        ind_rhomin = get_closest_ind(rho, rhomin)
        ind_rhomax = get_closest_ind(rho, rhomax)
        

        if rhoavg != 0:
            array_loc = np.mean(array[ind_rhomin:ind_rhomax, ind_tmin:ind_tmax], axis=0)
        else:
            assert(ind_rhomin==ind_rhomax)
            array_loc = array[ind_rhomin, ind_tmin:ind_tmax]
        
        if 'label' in kwargs:
            label=kwargs['label']
        else:
            label=r'$\rho$={:.2f}'.format(rho[int((ind_rhomin+ind_rhomax)/2)])
        if 'color' in kwargs:
            color=kwargs['color']
        else:
            color=clist[i]
        ax.plot(t, array_loc, color=color, marker= marker_list[i],markersize=3, label=label)
        
        if array_err is not None:
            if rhoavg != 0: 
                array_err_loc = np.mean(array_err[ind_rhomin:ind_rhomax, ind_tmin:ind_tmax], axis=0)
            else:
                assert(ind_rhomin==ind_rhomax)
                array_err_loc = array_err[ind_rhomin, ind_tmin:ind_tmax]
            
            ax.fill_between(t, array_loc-array_err_loc, array_loc+array_err_loc, color=color, alpha=0.2)
    if legend:
        ax.legend(fontsize='8')
    ax.set_xlabel('T [s]')
    # ax.ticklabel_format(axis='y', style='sci', scilimits=(19,19))

    return ax  



def plot_prof(array, t, rho, time_list, array_err = None, ax=None, rhomin=0, rhomax=1, tavg=0.1):
    '''
    TO PLOT IF THE ARRAY IS (RHO, TIME)
    Otherwise call the same function with np.transpose(array) & np.transpose(array_err)
    '''
    
    clist = get_cmap_list(len(time_list), cmap)
    ind_rhomin = get_closest_ind(rho, rhomin)
    ind_rhomax = get_closest_ind(rho, rhomax)+1            
    rho = rho[ind_rhomin:ind_rhomax]
    
    if ax is None:
        _, ax = plot_1d([], [], grid=True)
    
    for i,timeloc in enumerate(time_list):
        tinit = timeloc - tavg/2
        tfin = timeloc + tavg/2
        ind_tinit = get_closest_ind(t, tinit)
        ind_tfin = get_closest_ind(t, tfin)

        if tavg != 0:
            array_loc = np.mean(array[ind_rhomin:ind_rhomax, ind_tinit:ind_tfin], axis=1)
        else:
            assert(ind_tinit==ind_tfin)
            array_loc = array[ind_rhomin:ind_rhomax, ind_tinit]
            
        ax.plot(rho, array_loc, color=clist[i], marker='+', label='t={:.2f}'.format(t[int((ind_tinit+ind_tfin)/2)]))
        
        if array_err is not None:
            if tavg != 0: 
                array_err_loc = np.mean(array_err[ind_rhomin:ind_rhomax, ind_tinit:ind_tfin], axis=1)
            else:
                assert(ind_tinit==ind_tfin)
                array_err_loc = array_err[ind_rhomin:ind_rhomax, ind_tinit]
            
            ax.fill_between(rho, array_loc-array_err_loc, array_loc+array_err_loc, color=clist[i], alpha=0.2)

    ax.legend(fontsize='12')
    ax.set_xlabel(r'$\rho$')
    # ax.ticklabel_format(axis='y', style='sci', scilimits=(19,19))

    return ax 
    
### ================= ###
### get data from tcv ###
### ================= ###
class TCVShot():
    # NOTE: Not all data is available from any Lac#, connect to standard LAC to be sure
    # data = MDSplus.Data.execute(query)
    
    def __init__(self, shot, verbose=False):
        self.shot = shot
        self.verbose = verbose

        print('')
        print('Opening shot ', shot)
        try:
            self.tree = mds.Tree('tcv_shot', self.shot)
            self.tag = True
        except:
            print('Impossible to open shot {}'.format(shot))
            self.tag = False
            pass
        
        self.isH = isHydrogen(shot)
        if self.isH==True:
            self.m_i = m_proton
        else:
            self.m_i = 2*m_proton
        
        self.cmap = cmap
        
        self.tag_th_fit=False
        
        #mag eq
        self.tag_eq_mag = False
        self.tag_btor=False
        self.tag_q=False

    def get_thomson_fit(self):
        print('\n Loading Thomson FIT')
        #electron density shape(rho,t)
        try:
            self.th_ne = self.tree.getNode('\RESULTS::thomson.profiles.auto:ne').data()
            self.th_ne_err = self.tree.getNode('\RESULTS::thomson.profiles.auto:ne:error_bar').data()
            self.tag_th_ne = True
        except:
            print('No electron density from Thomson')
            self.tag_th_ne = False
        #electron temperature shape(rho,t)
        try:
            self.th_te = self.tree.getNode('\RESULTS::thomson.profiles.auto:te').data()
            self.th_te_err = self.tree.getNode('\RESULTS::thomson.profiles.auto:te:error_bar').data()
            self.tag_th_te = True
        except:
            print('No electron temperature from Thomson')
            self.tag_th_te = False
        try:
            #time   
            self.th_time = self.tree.getNode('\RESULTS::thomson.profiles.auto:time').data()
            #rho 
            self.th_rho = self.tree.getNode('\RESULTS::thomson.profiles.auto:rho').data()
            self.tag_th_time = True
        except:
            print('No time / rho from Thomson')
            self.tag_th_time = False
            
        #tag thomson True indicates the user already tried to load thomson data
        self.tag_th_fit = True
        
        
    def get_nbi(self):
        print('\n Loading NBI')
        #NB1 launcher
        try:
            self.nb1 = self.tree.getNode('\RESULTS::NB1:POWR_TCV').data()
            self.tag_nb1 = True
            
            self.nb1_ref = self.tree.getNode('\ATLAS::NB1.DATA.MODEL:POWR_NEUTRAL').data()
            self.nb1_ref_time = self.tree.getNode('\ATLAS::NB1.DATA.MODEL:POWR_NEUTRAL').dim_of().data()
            
        except:    
            print('No NB1 data')
            self.tag_nb1 = False
        #NB2 launcher
        try:
            self.nb2 = self.tree.getNode('\RESULTS::NB2:POWR_TCV').data()
            self.tag_nb2 = True
            
            self.nb2_ref = self.tree.getNode('\ATLAS::NB2.DATA.MODEL:POWR_NEUTRAL').data()
            self.nb2_ref_time = self.tree.getNode('\ATLAS::NB2.DATA.MODEL:POWR_NEUTRAL').dim_of().data()
            
        except:
            print('No NB2 data')
            self.tag_nb2 = False
        #NBI time
        try:
            self.nbi_time = self.tree.getNode('\ATLAS::NBH.DATA.MAIN_ADC:DATA').dim_of().data()
            self.tag_nbi = True
        except:
            print('No NBI time')
            self.tag_nbi = False
        
        
    def get_ecrh(self, set_nan_to_zero = True):
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
        print('\n Loading ECRH')

        #loading real values of ECH ?
        try:
            self.ecrh = self.tree.getNode('\RESULTS::TORAY.INPUT:P_GYRO').data()
            self.ecrh_tot = self.ecrh[11,:]
            self.ecrh_time = self.tree.getNode('\RESULTS::TORAY.INPUT:P_GYRO').dim_of().data()
            self.tag_ecrh = True
        except:
            print('No ECRH data')
            self.tag_ecrh = False

        #loading value of power cluster =>  those that we program (in kW)
        try:
            self.ecrh_cluster_a = self.tree.getNode('\pcs::power_cluster_A').data()/1e3
            self.ecrh_cluster_b = self.tree.getNode('\pcs::power_cluster_B').data()/1e3
            self.ecrh_cluster_c = self.tree.getNode('\pcs::power_cluster_C').data()/1e3
            self.ecrh_cluster_a_time = self.tree.getNode('\pcs::power_cluster_A').dim_of().data()
            self.ecrh_cluster_b_time = self.tree.getNode('\pcs::power_cluster_B').dim_of().data()
            self.ecrh_cluster_c_time = self.tree.getNode('\pcs::power_cluster_C').dim_of().data()
            self.tag_ecrh_cluster = True
        except:
            print('No ECRH cluster data')
            self.tag_ecrh_cluster = False

    def get_cxrs_fit(self, set_nan_to_zero=True):
        '''
        data from cxrs: ion temperature ti / toroidal velocity vtor
        shape: (time, rho)
        '''
        
        print('\n Loading CXRS')
    
        # Ion temperature
        try:
            self.cxrs_ti = self.tree.getNode('\RESULTS::CXRS.PROFFIT:TI').data()
            self.cxrs_ti_err = self.tree.getNode('\RESULTS::CXRS.PROFFIT:TI:ERR').data()
            self.tag_cxrs_ti = True    
        except:
            print('No CXRS Ti data')
            self.tag_cxrs_ti = False
            
        # Toroidal velocity
        try:
            self.cxrs_vtor = self.tree.getNode('\RESULTS::CXRS.PROFFIT:VI_TOR').data()
            self.cxrs_vtor_err = self.tree.getNode('\RESULTS::CXRS.PROFFIT:VI_TOR:ERR').data()
            #rho vtor and time vtor = to ti time and rho ?
            self.cxrs_vtor_rho = self.tree.getNode('\RESULTS::CXRS.PROFFIT:VI_TOR').getDimensionAt(0).data()
            self.cxrs_vtor_time = self.tree.getNode('\RESULTS::CXRS.PROFFIT:VI_TOR').getDimensionAt(1).data()
            self.tag_cxrs_vtor = True
        except:
            print('No CXRS vtor data')
            self.tag_cxrs_vtor = False
        try:
            self.cxrs_time = self.tree.getNode('\RESULTS::CXRS.PROFFIT:TI').getDimensionAt(1).data()
            self.cxrs_rho = self.tree.getNode('\RESULTS::CXRS.PROFFIT:TI').getDimensionAt(0).data()
            self.tag_cxrs = True
            
            #Check if all the rho measurements are equal
            are_lines_equal = np.all(self.cxrs_rho[1:] == self.cxrs_rho[:-1], axis=0).all()
            if are_lines_equal:
                self.cxrs_rho = self.cxrs_rho[0,:] #if all lines are equal take the first one
            else:
                print('Careful: cxrs_rho change during the shot')
        except: 
            print('No CXRS time / rho')
            self.tag_cxrs = False
        
        self.tag_cxrs = True
        
        if set_nan_to_zero:
            #Careful: this assume ti and ti_err are Nans at the same time
            if self.tag_cxrs_ti:
                self.cxrs_ti, _ = handle_nan(self.cxrs_ti, self.cxrs_time, value=None)
                self.cxrs_ti_err, self.cxrs_time = handle_nan(self.cxrs_ti_err, self.cxrs_time, value=None)
            if self.tag_cxrs_vtor:
                self.cxrs_vtor, self.cxrs_vtor_time = handle_nan(self.cxrs_vtor, self.cxrs_vtor_time, value=None)
            


    def get_cxrs_raw(self):
        '''
        TO BE MODIFIED 
        THIS VERSION HAS BEEN TAKEN FROM MATTEO 
        '''
        print('\n Loading CRXS raw data')
        
        try:
            self.VtorNode_raw = self.tree.getNode(r'\results::cxrs:vi_tor')
            self.VpolNode_raw = self.tree.getNode(r'\results::cxrs:vi_pol')
            self.tag_cxrs_raw = True
        except:
            self.tag_cxrs_raw = False
            return

        self.cxrs_rho_vtor_raw = self.VtorNode_raw.getDimensionAt(0).data()
        self.cxrs_time_vtor_raw = self.VtorNode_raw.getDimensionAt(1).data()
        self.cxrs_vtor_raw = self.VtorNode_raw.data()
        self.cxrs_vtor_err_raw = self.tree.getNode('\RESULTS::CXRS:VI_TOR:ERR').data()
        self.cxrs_ni_raw = self.tree.getNode('\RESULTS::CXRS.NI').data()
        self.cxrs_ti_raw = self.tree.getNode('\RESULTS::CXRS.TI').data()
        self.cxrs_ni_err_raw = self.tree.getNode('\RESULTS::CXRS:NI:ERR').data()
        self.cxrs_ti_err_raw = self.tree.getNode('\RESULTS::CXRS:TI:ERR').data()
        self.cxrs_time_vpol_raw = self.VpolNode_raw.getDimensionAt(1).data()
        self.cxrs_vpol_err_raw = self.tree.getNode('\RESULTS::CXRS:VI_POL:ERR').data()
        self.cxrs_vpol_raw = self.VpolNode_raw.data()
        self.cxrs_rho_vpol_raw = self.VpolNode_raw.getDimensionAt(0).data()

    def get_Ip(self):
        '''
        For now the magnetic equilibrium / data are not loaded for pythonn
        '''
        
        print('\n Loading Ip')
        
        try:
            self.ip = self.tree.getNode('\MAGNETICS::IPLASMA').data()
            self.ip_time = self.tree.getNode('\MAGNETICS::IPLASMA').getDimensionAt(0).data()
            self.tag_ip = True
        except:
            print('No Ip data')
            self.tag_ip = False
        
        try:
            self.ip_ref = self.tree.getNode('\PCS::DRAW_REFS:REF_013').data()
            self.ip_ref_time = self.tree.getNode('\PCS::DRAW_REFS:REF_013').getDimensionAt(0).data()
            self.tag_ip_ref = True
        except:
            print('No Ip ref data')
            self.tag_ip_ref = False


    def get_FIR(self):
        
        print('\n Loading FIR data')
        
        try:
            self.fir_int_ne = self.tree.getNode('\RESULTS::FIR:LIN_INT_DENS').data()
            self.fir_time = self.tree.getNode('\RESULTS::FIR:LIN_INT_DENS').getDimensionAt(0).data()
            self.tag_fir = True
        except:
            print('No FIR data')
            self.tag_fir = False
            
        try: 
            self.ref_int_ne = self.tree.getNode('\PCS::DRAW_REFS:REF_021').data()
            self.ref_int_ne_time = self.tree.getNode('\PCS::DRAW_REFS:REF_021').getDimensionAt(0).data()
        except:
            print('No reference int ne data')
            self.tag_ref_int_ne = False


    def get_h_factor(self, nb_scaling_law = 4):
        '''
        nb_scaling law chooses which scaling law to use for the H factor
        0=RLW
        1=IAEA-TCV
        2=ITER98L
        3=ITER89P
        4=ITER98
        default = 4
        '''
        _list_model = ['RLW', 'IAEA-TCV', 'ITER98L', 'ITER89P', 'ITER98']
        
        print('\n Loading H factor data ; model = {}'.format(_list_model[nb_scaling_law]))
        
        try:
            self.h_factor = self.tree.getNode('\RESULTS::PROFFIT:AVG_TIME:H_SCAL').data()[:,nb_scaling_law]
            self.h_factor_time = self.tree.getNode('\RESULTS::PROFFIT:AVG_TIME:H_SCAL').getDimensionAt(1).data()
            self.h_factor_model = _list_model[nb_scaling_law]
            self.tag_h_factor = True
        except:
            print('No H factor data')
            self.tag_h_factor = False
            self.h_factor_model = ''


    def get_mag_eq_object(self, time):
        print('\n Loading magnetic equilibrium object at t={}'.format(time))
        
        try:
            equil_object = Equilibrium2d('tcv', self.shot, time)
            mag_eq = equil_object.from_shot('tcv', self.shot, time)
            self.mag_eq = mag_eq
            
            psiN = ((self.mag_eq['psi'] - self.mag_eq['psiaxe']) / (self.mag_eq['psibnd'] - self.mag_eq['psiaxe']))
            self.mag_eq['rho_psi'] = np.sqrt(psiN)
            
            self.tag_mag_eq = True
        except:
            print('No magnetic equilibrium object')
            self.mag_eq = None
            self.tag_mag_eq = False


    def get_btor(self, time):
        
        '''
        btor is a function of (r,z) but usually does not depend much on z
        '''
        print('\n Loading Btor at t={}'.format(time))
        
        try:
            self.get_mag_eq_object(time)
            self.btor = self.mag_eq['Btor']
            self.mag_eq_rgrid = self.mag_eq['rgrid']
            self.mag_eq_zgrid = self.mag_eq['zgrid']
            self.tag_btor = True
        except:
            print('No Btor data')
            self.tag_btor = False
            self.btor = None

    def get_rho_s_r_z(self, time_window, r, z, rho):
        '''
        sound larmor radius from thomson temperature and mag eq
        time_window is a list: [tinit:tfin]
        SHOULD BE POSSIBLE TO ESTIMATE RHO FROM R AND Z
        '''
        
        if self.tag_btor==False or self.tag_th_fit==False:
            try:
                self.get_thomson_fit()
                self.get_btor(np.mean(time_window))
            except:
                print(" can't load larmor radius: check btor or Thomson data")
                self.tag_rho_s = False
                self.rho_s = None         

        th_te_ind_init = get_closest_ind(self.th_time, time_window[0])
        th_te_ind_fin = get_closest_ind(self.th_time, time_window[1])
        
        th_te_prof_loc = np.mean(self.th_te[:, th_te_ind_init:th_te_ind_fin], axis=1)
        
        rho_te_loc = get_closest_ind(self.th_rho, rho)
        
        th_te_loc = th_te_prof_loc[rho_te_loc]
        
        cs = np.sqrt(e*th_te_loc/(self.m_i))  #sound speed in m/s
        print('ok sound speed')
        
        r_ind = get_closest_ind(self.mag_eq_rgrid, r)
        z_ind = get_closest_ind(self.mag_eq_zgrid, z)
        
        _btor_loc = abs(self.btor[r_ind, z_ind])

        # _btor_loc = np.mean(self.btor, axis=1) #Avg btor over z to have radial estimate
        print('ok btor ; shape:', _btor_loc.shape)
    
        
        omega_cs = (e*_btor_loc)/(self.m_i) #frequency in s**-1
        print('ok cycl freq: ', omega_cs)
        # return cs, omega_cs
        # self.rho_s = cs/omega_cs #sound larmor radius in m
        # print('ok rho_s')
        # self.tag_rho_s = True
        # print('ok tag')
        
        return cs/omega_cs
        
        
        
    def get_r_and_R(self, twindow=None):
        ''' 
        Careful => long due to matlab call 
        Taken from S. Rienacker and modified to fit in this class
        '''
        
        print('\n Loading magnetic position from tcv_coordinate_conversion at T = {}'.format(twindow))
    
        from matlabtools import run_matlab_code
        
        #We take rho from thomson ?
        if self.tag_th_fit is False:
            self.get_thomson_fit()

        #Getting R from matlab equilibrium
        _time = np.array(twindow).mean()
        _time_str = f'{_time:.2f}'
        _rho_str = f"linspace({self.th_rho[0]},{self.th_rho[-1]},{self.th_rho.size})"
        _theta_str = f"[0,90,180,270]"
        ml_code = f"[result.coord,result.time] = tcv_coordinate_conversion({self.shot},{_time_str},{_rho_str},{_theta_str},'p','LIUQE.M');"
        result = run_matlab_code(ml_code, logdir=None)
    
        R = np.mean(result.coord.R, axis=1)
        
        #Small r from mds
        r_axis = self.tree.tdiExecute(f'tcv_eq("r_axis", "liuqe.m")')
        z_axis = self.tree.tdiExecute(f'tcv_eq("z_axis", "liuqe.m")')

        _time = r_axis.getDimensionAt(0).data()
        it = (_time > twindow[0]) & (_time < twindow[1])
        r_axis = r_axis.data()[it]
        z_axis = z_axis.data()[it]

        r_axis = np.nanmean(r_axis)
        z_axis = np.nanmean(z_axis)

        d = np.sqrt((result.coord.R - r_axis)**2 + (result.coord.z - z_axis)**2) #poloidal small radius
        r = np.mean(d, axis=1) #poloidal avg 
        
        self.r_polavg = r
        self.R_polavg = R
        
        self.tag_r_and_R = True
        self.r_and_R_twindow = twindow
        
        return r, R


    def get_q(self):
        
        print('\n Loading safety factor')

        #safety factor
        key=('Q_PSI')
        self.q_psi = self.tree.tdiExecute(f'tcv_eq("{key}", "liuqe.m")').data()
        self.q_time = self.tree.tdiExecute(f'tcv_eq("{key}", "liuqe.m")').getDimensionAt(1).data()
        self.q_rho = self.tree.tdiExecute(f'tcv_eq("{key}", "liuqe.m")').getDimensionAt(0).data()

        self.tag_q = True



    def get_collisionality_prof(self, time_window, ax=None, rhomin=0, rhomax=1):
        
        if self.tag_th_fit is False:
            self.get_thomson_fit()
        
        self.get_r_and_R(time_window)
        # if time_window != self.r_and_R_twindow:
        #     self.get_r_and_R(time_window)
        
        if self.tag_q is False:
            self.get_q()
        
        
        #Loading Thomson data
        t = self.th_time
        rho = self.th_rho
        
        ind_tmin = get_closest_ind(t, time_window[0])
        ind_tmax = get_closest_ind(t, time_window[1])
        
        ind_rhomin = get_closest_ind(rho, rhomin)
        ind_rhomax = get_closest_ind(rho, rhomax)+1
        rho = rho[ind_rhomin:ind_rhomax]
        
        te = np.mean(self.th_te[ind_rhomin:ind_rhomax, ind_tmin:ind_tmax], axis=1)
        ne = np.mean(self.th_ne[ind_rhomin:ind_rhomax, ind_tmin:ind_tmax], axis=1)/10**19
        
        
        #Loading magnetic equilibrium
        ind_q_tmin = get_closest_ind(self.q_time, time_window[0])
        ind_q_tmax = get_closest_ind(self.q_time, time_window[1])
        
        q_psi = np.mean(self.q_psi[ind_q_tmin:ind_q_tmax, ind_rhomin:ind_rhomax], axis=0)
        
        r = self.r_polavg[ind_rhomin:ind_rhomax]
        R = self.R_polavg[ind_rhomin:ind_rhomax]
        inverse_aspect_ratio = r/R
              	
               
        #Loading Zeff
        zeff = np.mean(self.tree.getNode(r'\results::ibs:z_eff').data())
        self.zeff = zeff
        
        ### VARIABLES TCV ###
        # R = 0.88
        # a = 0.25
        # Z = 1
        
        # vth = np.sqrt(2*te/m_e)
        vth = 1.32620*10**7*np.sqrt(te/1e3)
        
        
        #cs = np.sqrt(e*te/m_i)  #sound speed in m/s 
        # print("cs : " + "{:e}".format(cs) + " in m/s")
        # omega_cs = (e*B)/(m_i) #frequency in s**-1
        # print("omega_cs : " + "{:e}".format(omega_cs) + " in 1/s")
        # rho_s = cs/omega_cs #sound larmor radius in m
        
        #inverse_aspect_ratio = rho*a/(R)
        
        
        lnLambda = 15 - 0.5*np.log(10*ne) + np.log(1E-3 * te) #Coulomb logarithm from X.Garbet magic booklet 
        self.lnLambda = lnLambda
        
        #Non-normalized (to divide by omega_cs or transit time ?)
        self.nu_ei = ((4*np.sqrt(2*np.pi))/3) * ((e**2/(4*np.pi*epsilon_0))**2) * np.mean(lnLambda) * ((ne*10**19)*zeff/(np.sqrt(m_e)*((e*te)**(3/2)))) #electron-ion collision frequency (from X.Garbet magic booklet)
        
        
        
        # print(((4*np.sqrt(2*np.pi))/3) * ((e**2/(4*np.pi*epsilon_0))**2)*np.mean(lnLambda)/(np.sqrt(m_e)))

        self.nu_star = self.nu_ei * (q_psi*R/(vth*inverse_aspect_ratio**(3/2)))
        self.inverse_aspect_ratio = inverse_aspect_ratio
        self.nu_star_rho = rho
        
        self.nu_star_q = q_psi
        self.nu_star_r = r
        self.nu_star_R = R
        
        if ax is None:
            fig,ax= plot_1d([], [], grid=True, xlabel=r'$\rho$', ylabel=r'$\nu_*$')
        ax.set_title('#{}'.format(self.shot))
        ax.plot(rho, self.nu_star,marker='+', label=r' T = {} s'.format(time_window))
        
        # ax.plot(rho, ne/10 , label='ne')
        # # ax.plot(rho, te/1000, label='te')
        # ax.plot(rho, q_psi/10, label='q/10')
        # ax.plot(rho, inverse_aspect_ratio, label='epsilon')
        
        # ax.plot(rho, self.inverse_aspect_ratio**(3/2), label=r'$\epsilon^{3/2}$', color='black')
        ax.legend()
             
        return ax
               
               
               
    def plot_thomson_prof(self, time_list, rhomin=0, rhomax=1, tavg=0.2):
        
        clist = get_cmap_list(len(time_list), self.cmap)
        
        if self.tag_th_fit==True:

            t = self.th_time
            rho = self.th_rho
            
            ind_rhomin = get_closest_ind(rho, rhomin)
            ind_rhomax = get_closest_ind(rho, rhomax)+1            
            rho = rho[ind_rhomin:ind_rhomax]

            fig, ax = prep_multiple_subplots(2, 1, figsize=None, axgrid=[0,1], constrained_layout=True, sharex=True, sharey=False)
            ax[0].set_title('#{} ; {:.1f} s avg'.format(self.shot, tavg))
            
            if self.tag_th_ne:
                ne = self.th_ne
                ne_err = self.th_ne_err
                
                for i,timeloc in enumerate(time_list):
                    tinit = timeloc - tavg/2
                    tfin = timeloc + tavg/2
                    ind_tinit_thomson = get_closest_ind(t, tinit)
                    ind_tfin_thomson = get_closest_ind(t, tfin)

                    ne_loc = np.mean(ne[ind_rhomin:ind_rhomax, ind_tinit_thomson:ind_tfin_thomson], axis=1)
                    ne_err_loc = np.mean(ne_err[ind_rhomin:ind_rhomax, ind_tinit_thomson:ind_tfin_thomson], axis=1)

                    ax[0].plot(rho, ne_loc, color=clist[i], marker='+', label='t={:.2f}'.format(timeloc))
                    ax[0].fill_between(rho, ne_loc-ne_err_loc, ne_loc+ne_err_loc, color=clist[i], alpha=0.2)
                    ax[0].set_ylabel(r' $n_e$ $[10^{19} m^{-3}]$')
                    ax[0].ticklabel_format(axis='y', style='sci', scilimits=(19,19))
                    ax[0].legend(fontsize='12')
            else:
                ax[0].text(0.5, 0.5, 'No ne data', horizontalalignment='center',
                        verticalalignment='center', transform=ax.transAxes)
                ax[0].set_ylabel(r'$n_e$ $[10^{19} m^{-3}]$')
            
            if self.tag_th_te:
                te = self.th_te
                te_err = self.th_te_err
                
                for i,timeloc in enumerate(time_list):
                    tinit = timeloc - tavg/2
                    tfin = timeloc + tavg/2
                    ind_tinit_thomson = get_closest_ind(t, tinit)
                    ind_tfin_thomson = get_closest_ind(t, tfin)

                    te_loc = np.array(np.mean(te[ind_rhomin:ind_rhomax, ind_tinit_thomson:ind_tfin_thomson], axis=1))
                    te_err_loc = np.array( np.mean(te_err[ind_rhomin:ind_rhomax, ind_tinit_thomson:ind_tfin_thomson], axis=1))
                    ax[1].plot(rho, te_loc, color=clist[i], marker='+')
                    ax[1].fill_between(rho, te_loc-te_err_loc, te_loc+te_err_loc, color=clist[i], alpha=0.2)
                    ax[1].set_xlabel(r'$\rho$')
                    ax[1].set_ylabel(r'$T_e$ $[eV]$')
                    ax[1].ticklabel_format(axis='y', style='sci', scilimits=(2,2))
            else:
                ax[1].text(0.5, 0.5, 'No Te data', horizontalalignment='center',
                        verticalalignment='center', transform=ax.transAxes)
                ax[1].set_xlabel(r'$\rho$')
                ax[1].set_ylabel(r'$T_e$ $[eV]$', color='red')
            
        else:
            self.get_thomson_fit()
            self.plot_thomson_prof(time, rhomin, rhomax)
            
        
            
          
    def plot_thomson_time(self, rho_list, tmin=0, tmax=-1, rho_avg=0.1):
        
        clist = get_cmap_list(len(rho_list), self.cmap)
        
        if self.tag_th_fit==True:

            t = self.th_time
            rho = self.th_rho
            nt = len(t)
            nx = len(rho)
            if rho_list[0]=='all':
                rho_list=rho

            ind_tmin = get_closest_ind(t, tmin)
            
            if tmax==-1:
                ind_tmax=-1
            else:
                ind_tmax = get_closest_ind(t, tmax)+1            
            
            t = t[ind_tmin:ind_tmax]


            fig, ax = prep_multiple_subplots(2, 1, figsize=None, axgrid=[0,1], constrained_layout=True, sharex=True, sharey=False)
            ax[0].set_title('#{} ; {:.1f} rho avg'.format(self.shot, rho_avg))
            
            if self.tag_th_ne:
                ne = self.th_ne
                ne_err = self.th_ne_err
                
                for i, rholoc in enumerate(rho_list):
                    rleft = rholoc - rho_avg/2
                    rright = rholoc + rho_avg/2
                    ind_rleft = get_closest_ind(rho, rleft)
                    ind_rright = get_closest_ind(rho, rright)

                    ne_loc = np.mean(ne[ind_rleft:ind_rright, ind_tmin:ind_tmax], axis=0)
                    ne_err_loc = np.mean(ne_err[ind_rleft:ind_rright, ind_tmin:ind_tmax], axis=0)

                    ax[0].plot(t, ne_loc, color=clist[i], marker='+', label=r'$\rho={:.2f}$'.format(rholoc))
                    ax[0].fill_between(t, ne_loc-ne_err_loc, ne_loc+ne_err_loc, color=clist[i], alpha=0.2)
                    ax[0].set_ylabel(r' $n_e$ $[10^{19} m^{-3}]$')
                    ax[0].ticklabel_format(axis='y', style='sci', scilimits=(19,19))
                    ax[0].legend(fontsize='12')
            else:
                ax[0].text(0.5, 0.5, 'No ne data', horizontalalignment='center',
                        verticalalignment='center', transform=ax.transAxes)
                ax[0].set_ylabel(r'$n_e$ $[10^{19} m^{-3}]$')
            
            if self.tag_th_te:
                te = self.th_te
                te_err = self.th_te_err
                
                for i,rholoc in enumerate(rho_list):
                    rleft = rholoc - rho_avg/2
                    rright = rholoc + rho_avg/2
                    ind_rleft = get_closest_ind(rho, rleft)
                    ind_rright = get_closest_ind(rho, rright)

                    te_loc = np.array(np.mean(te[ind_rleft:ind_rright, ind_tmin:ind_tmax], axis=0))
                    te_err_loc = np.array(np.mean(te_err[ind_rleft:ind_rright, ind_tmin:ind_tmax], axis=0))
                    ax[1].plot(t, te_loc, color=clist[i], marker='+')
                    ax[1].fill_between(t, te_loc-te_err_loc, te_loc+te_err_loc, color=clist[i], alpha=0.2)
                    ax[1].set_xlabel(r'$t$')
                    ax[1].set_ylabel(r'$T_e$ $[eV]$')
                    ax[1].ticklabel_format(axis='y', style='sci', scilimits=(2,2))
            else:
                ax[1].text(0.5, 0.5, 'No Te data', horizontalalignment='center',
                        verticalalignment='center', transform=ax.transAxes)
                ax[1].set_xlabel(r'$t$')
                ax[1].set_ylabel(r'$T_e$ $[eV]$', color='red')
            
        else:
            self.get_thomson_fit()
            plot_thomson_time(shot, rho_list, tmin, tmax, rho_avg)

            
    
    def plot_heating(self, ax=None, color=None, **kwargs):
        self.get_nbi()
        self.get_ecrh()
        
        if ax is None:
            fig, ax = plot_1d([], [], grid=True, **kwargs)
            ax.set_title('#{}'.format(self.shot))
        
        if self.tag_nbi:
            if self.tag_nb1:
                if color is None:
                    ax.plot(self.nbi_time, self.nb1*1e3, color='red', label='NB1', **kwargs)
                    ax.plot(self.nb1_ref_time, self.nb1_ref*1e3, color='xkcd:dark red', label='NB1 ref')
                else:
                    ax.plot(self.nbi_time, self.nb1*1e3,linestyle='--', marker='o', linewidth=0.2, markersize=5, markevery=500, color=color, label='NB1', **kwargs)
                
        if self.tag_nb2:
            if color is None:
                ax.plot(self.nbi_time, self.nb2*1e3, color='green', label='NB2', **kwargs)
                ax.plot(self.nb2_ref_time, self.nb2_ref*1e3, color='xkcd:dark green', label='NB2 ref')
            else:
                ax.plot(self.nbi_time, self.nb2*1e3, linestyle='dotted', color=color, label='NB2', **kwargs)
        
        if self.tag_nb1 and self.tag_nb2:
            if color is None:
                ax.plot(self.nbi_time, self.nb1*1e3 + self.nb2*1e3, color='blue', label='NBtot',**kwargs)
            else:
                ax.plot(self.nbi_time, self.nb1*1e3 + self.nb2*1e3, linestyle='-.', color=color, label='NBtot', **kwargs)
                
        if self.tag_ecrh:
            if color is None:
                ax.plot(self.ecrh_cluster_a_time, self.ecrh_cluster_a, color='xkcd:dark orange', label='ECRH ref a')
                ax.plot(self.ecrh_cluster_b_time, self.ecrh_cluster_b, color='xkcd:dirty orange', label='ECRH ref b')
                ax.plot(self.ecrh_time, self.ecrh_tot, color='xkcd:dark yellow', label='ECRH tot', **kwargs)
            else:
                ax.plot(self.ecrh_time, self.ecrh_tot, marker='s',linewidth=0.2, markersize=5, markevery=500, color=color, label='ECRH tot', **kwargs)
                

        ax.set_ylabel('Power [kW]')
        ax.set_xlabel('T [s]')
        ax.legend(fontsize='12')
        if self.tag_ecrh and self.tag_nbi:
            ax.set_xlim(max(self.ecrh_time[0], self.nbi_time[0]), min(self.ecrh_time[-1], self.nbi_time[-1]))
        
        return ax

    def plot_summary(self):
        
        
        self.get_thomson_fit()
        self.get_cxrs_fit(set_nan_to_zero=False)
        self.get_Ip()
        self.get_FIR()
        self.get_h_factor()
        
        fig ,axs = prep_multiple_subplots(3,3, figsize=(12,8), axgrid=[0,1,2,3,4,5,6,7,8], sharex=False)
        fig.suptitle('#{}'.format(self.shot))
        
        #Ip
        if self.tag_ip:
            axs[0,0].plot(self.ip_time, self.ip*1e-6, color='xkcd:green', label='Ip')
        if self.tag_ip_ref:
            axs[0,0].plot(self.ip_ref_time, self.ip_ref*1e-6, color='black', label='Ip ref')
        axs[0,0].set_ylabel('Ip [MA]')
        axs[0,0].legend(fontsize='10')
        axs[0,0].set_xlabel('T [s]')
        #heating
        self.plot_heating(ax=axs[1,0])
        if self.tag_th_time:
            axs[1,0].set_xlim(self.th_time[0], self.th_time[-1])
        #H factor
        if self.tag_h_factor:
            axs[2,0].plot(self.h_factor_time, self.h_factor, color='xkcd:dark green', label='H factor')
        axs[2,0].set_ylabel(r'$H$ {}'.format(self.h_factor_model))
        axs[2,0].set_xlabel('T [s]')
        #ne
        if self.tag_th_ne:
            plot_time(self.th_ne/10**19, self.th_time, self.th_rho, [0.6,0.8,0.95], array_err=self.th_ne_err/10**19, ax=axs[0,1], tmin=0, tmax=-1, rhoavg=0.1)
        axs[0,1].set_ylabel(r'$n_e$ $[10^{19}]$')
        #FIR ne
        if self.tag_fir:
            axs[1,1].plot(self.fir_time, self.fir_int_ne/10**19, color='xkcd:dark green', label='FIR ne')
            axs[1,1].plot(self.ref_int_ne_time, self.ref_int_ne, color='black', label='FIR ne ref')
            axs[1,1].set_ylim(0, np.max(self.fir_int_ne/10**19)*1.1)
        axs[1,1].set_ylabel(r'FIR $n_e$ $[10^{19}]$')
        if self.tag_th_time:
            axs[1,1].set_xlim(self.th_time[0], self.th_time[-1])
        axs[1,1].set_xlabel('T [s]')
        axs[1,1].legend()
        #vtor 
        if self.tag_cxrs_vtor:
            plot_time(np.transpose(self.cxrs_vtor), self.cxrs_vtor_time, self.cxrs_rho, [0.6,0.8,0.95], array_err=np.transpose(self.cxrs_vtor_err), ax=axs[2,1], tmin=0, tmax=-1, rhoavg=0, legend=False)
        axs[2,1].set_ylabel(r'$v_{tor}$ $[km/s]$')
        #Te
        if self.tag_th_te:
            plot_time(self.th_te/1e3, self.th_time, self.th_rho, [0.6,0.8,0.95], array_err=self.th_te_err/1e3, ax=axs[0,2], tmin=0, tmax=-1, rhoavg=0.1, legend=False)
        axs[0,2].set_ylabel(r'$T_e$ $[keV]$')
        #Ti
        if self.tag_cxrs_ti:
            plot_time(np.transpose(self.cxrs_ti/1e3), self.cxrs_time, self.cxrs_rho, [0.6,0.8,0.95], array_err=np.transpose(self.cxrs_ti_err/1e3), ax=axs[1,2], tmin=0, tmax=-1, rhoavg=0, legend=False)
        axs[1,2].set_ylabel(r'$T_i$ $[keV]$')
        #Te/Ti
        if self.tag_th_te:
            plot_time(self.th_te, self.th_time, self.th_rho, [0.9], array_err=self.th_te_err, ax=axs[2,2], tmin=0, tmax=-1, rhoavg=0.1, label=r'$T_e$ $\rho=0.9$', color='blue')
        if self.tag_cxrs_ti:
            plot_time(np.transpose(self.cxrs_ti), self.cxrs_time, self.cxrs_rho, [0.9], array_err=np.transpose(self.cxrs_ti_err), ax=axs[2,2], tmin=0, tmax=-1, rhoavg=0, label=r'$T_i$, $\rho=0.9$', color='red')
            axs[2,2].set_xlim(self.cxrs_time[0], self.cxrs_time[-1])
        axs[2,2].set_ylabel(r'$T_e$, $T_i$ $[eV]$')
        
        
        axs[0,1].legend(fontsize='10')
    
        # plt.tight_layout()




def plot_profiles_comparison(shot_list, time_list, plot_err=None, rhomin=None, rhomax=None, tavg=None):
    
    list_obj_tcv = [TCVShot(shot) for shot in shot_list]
    
    for i in range(len(list_obj_tcv)):
        list_obj_tcv[i].get_thomson_fit()
        list_obj_tcv[i].get_cxrs_fit(set_nan_to_zero=True)
        
    if len(shot_list) != len(time_list):
        print('\n --- choosing same time for all shots --- ')
        th_time_ind_list = [get_closest_ind(list_obj_tcv[i].th_time, time_list[0]) for i in range(len(list_obj_tcv))]
    
    th_time_ind_list = [get_closest_ind(list_obj_tcv[i].th_time, time_list[i]) for i in range(len(list_obj_tcv))]
    

    # color_list = ['blue', 'red', 'green', 'black', 'orange', 'purple', 'brown', 'pink', 'grey', 'cyan']
    # marker_list = 


    fig ,axs = prep_multiple_subplots(2,2, figsize=(8,5), axgrid=[0,1,2,3], sharex=True)
    fig.suptitle('#{} ; tavg = {}'.format(shot_list, tavg))
    
    # fig, ax_te = plot_1d([], [], grid=True, xlabel=r'$\rho$', ylabel=r'$T_e$ $[eV]$')
    # fig, ax_ne = plot_1d([], [], grid=True, xlabel=r'$\rho$', ylabel=r'$n_e$ $[10^{19} m^{-3}]$')
    # fig, ax_ti = plot_1d([], [], grid=True, xlabel=r'$\rho$', ylabel=r'$T_i$ $[eV]$')
    for i in range(len(list_obj_tcv)):
        
        # ---rho range------------------------------------- #
        th_rho_all = list_obj_tcv[i].th_rho
        if rhomin is not None:
            rho_min_ind = get_closest_ind(th_rho_all, rhomin)
        else:
            rho_min_ind = None
        if rhomax is not None:
            rho_max_ind = get_closest_ind(th_rho_all, rhomax)+1
        else:
            rho_max_ind = None
        th_rho = th_rho_all[rho_min_ind:rho_max_ind]
        # ------------------------------------------------ #
        
        # ---time average--------------------------------- #
        if tavg is not None:
            tmin_ind = get_closest_ind(list_obj_tcv[i].th_time, list_obj_tcv[i].th_time[th_time_ind_list[i]]-tavg/2)
            tmax_ind = get_closest_ind(list_obj_tcv[i].th_time, list_obj_tcv[i].th_time[th_time_ind_list[i]]+tavg/2)
            # print(tmin_ind, tmax_ind)
            th_te = np.mean(list_obj_tcv[i].th_te[rho_min_ind:rho_max_ind, tmin_ind:tmax_ind], axis=1)
            th_ne = np.mean(list_obj_tcv[i].th_ne[rho_min_ind:rho_max_ind, tmin_ind:tmax_ind], axis=1)
            
            if plot_err:
                th_te_err = np.mean(list_obj_tcv[i].th_te_err[rho_min_ind:rho_max_ind, tmin_ind:tmax_ind], axis=1)
                th_ne_err = np.mean(list_obj_tcv[i].th_ne_err[rho_min_ind:rho_max_ind, tmin_ind:tmax_ind], axis=1)    
            else:
                th_te_err = np.zeros((np.shape(th_te)))
                th_ne_err = np.zeros((np.shape(th_ne)))
        # ------------------------------------------------ #
            
        else:
            th_te = list_obj_tcv[i].th_te[rho_min_ind:rho_max_ind,th_time_ind_list[i]]
            th_ne = list_obj_tcv[i].th_ne[rho_min_ind:rho_max_ind,th_time_ind_list[i]]
            if plot_err:
                th_te_err = list_obj_tcv[i].th_te_err[rho_min_ind:rho_max_ind,th_time_ind_list[i]]
                th_ne_err = list_obj_tcv[i].th_ne_err[rho_min_ind:rho_max_ind,th_time_ind_list[i]]
            else:
                th_te_err = np.zeros((np.shape(th_te)))
                th_ne_err = np.zeros((np.shape(th_ne)))
        
        # ---plotting thomson data------------------------ #                   
        # list_obj_tcv[i].plot_heating(ax=axs[0,0], color=color_list[i])
        axs[0,1].plot(th_rho, th_te/1e3 ,linewidth=2, markersize=4, fillstyle='none',
                      label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].th_time[th_time_ind_list[i]]),  color=color_list[i], marker='')
        axs[0,1].fill_between(th_rho, (th_te-th_te_err)/1e3, (th_te+th_te_err)/1e3, color=color_list[i], alpha=0.2)    
    
        axs[0,0].plot(th_rho, th_ne/10**19,linewidth=2, markersize=4, fillstyle='none',
                      label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].th_time[th_time_ind_list[i]]), color=color_list[i], marker='')
        axs[0,0].fill_between(th_rho, (th_ne-th_ne_err)/10**19, (th_ne+th_ne_err)/10**19, color=color_list[i], alpha=0.2)
        # ------------------------------------------------ #

    axs[0,0].legend(fontsize="8")
    # axs[0,1].legend()
    # axs[1,1].legend()
    axs[0,1].set_ylabel(r'$T_e$ $[keV]$')
    axs[0,0].set_ylabel(r'$n_e$ $[10^{19} m^{-3}]$')
    axs[1,0].set_ylabel(r'$v_{tor}$ $[km/s]$')
    axs[1,1].set_ylabel(r'$T_i$ $[keV]$')
    axs[1,0].set_xlabel(r'$\rho$')
    axs[1,1].set_xlabel(r'$\rho$')
    
    
    #Data CXRS 
    try:
        if len(shot_list) != len(time_list):
            print('\n --- choosing same time for all shots --- ')
            cxrs_time_ind_list = [get_closest_ind(list_obj_tcv[i].cxrs_time, time_list[0]) for i in range(len(list_obj_tcv))]
             
        cxrs_time_ind_list = [get_closest_ind(list_obj_tcv[i].cxrs_time, time_list[i]) for i in range(len(list_obj_tcv))]
        cxrs_vtor_time_ind_list = [get_closest_ind(list_obj_tcv[i].cxrs_vtor_time, time_list[i]) for i in range(len(list_obj_tcv))]
        
        for i in range(len(list_obj_tcv)):
            
            cxrs_rho_all = list_obj_tcv[i].cxrs_rho
            if rhomin is not None:
                rho_min_ind = get_closest_ind(cxrs_rho_all, rhomin)
            else:
                rho_min_ind = None
            if rhomax is not None:
                rho_max_ind = get_closest_ind(cxrs_rho_all, rhomax)+1
            else:
                rho_max_ind = None
            
            cxrs_rho = cxrs_rho_all[rho_min_ind:rho_max_ind]
            
            if tavg is not None:
                tmin_ind = get_closest_ind(list_obj_tcv[i].cxrs_vtor_time, list_obj_tcv[i].cxrs_vtor_time[cxrs_vtor_time_ind_list[i]]-tavg/2)
                tmax_ind = get_closest_ind(list_obj_tcv[i].cxrs_vtor_time, list_obj_tcv[i].cxrs_vtor_time[cxrs_vtor_time_ind_list[i]]+tavg/2)
                # print(tmin_ind, tmax_ind)*
                
                cxrs_vtor = np.mean(list_obj_tcv[i].cxrs_vtor[tmin_ind:tmax_ind, rho_min_ind:rho_max_ind], axis=0)
                
                tmin_ind = get_closest_ind(list_obj_tcv[i].cxrs_time, list_obj_tcv[i].cxrs_time[cxrs_time_ind_list[i]]-tavg/2)
                tmax_ind = get_closest_ind(list_obj_tcv[i].cxrs_time, list_obj_tcv[i].cxrs_time[cxrs_time_ind_list[i]]+tavg/2)
               
                cxrs_ti = np.mean(list_obj_tcv[i].cxrs_ti[tmin_ind:tmax_ind, rho_min_ind:rho_max_ind], axis=0)
                
                if plot_err:
                    cxrs_vtor_err = np.mean(list_obj_tcv[i].cxrs_vtor_err[tmin_ind:tmax_ind, rho_min_ind:rho_max_ind], axis=0)
                    cxrs_ti_err = np.mean(list_obj_tcv[i].cxrs_ti_err[tmin_ind:tmax_ind, rho_min_ind:rho_max_ind], axis=0)
                else:
                    cxrs_vtor_err = np.zeros((np.shape(cxrs_vtor)))
                    cxrs_ti_err = np.zeros((np.shape(cxrs_ti)))
                
            else:
                cxrs_vtor = list_obj_tcv[i].cxrs_vtor[cxrs_vtor_time_ind_list[i], rho_min_ind:rho_max_ind]
                cxrs_ti = list_obj_tcv[i].cxrs_ti[cxrs_time_ind_list[i], rho_min_ind:rho_max_ind]

                if plot_err:
                    cxrs_vtor_err = list_obj_tcv[i].cxrs_vtor_err[cxrs_vtor_time_ind_list[i], rho_min_ind:rho_max_ind]
                    cxrs_ti_err = list_obj_tcv[i].cxrs_ti_err[cxrs_time_ind_list[i], rho_min_ind:rho_max_ind]
                else:
                    cxrs_vtor_err = np.zeros((np.shape(cxrs_vtor)))
                    cxrs_ti_err = np.zeros((np.shape(cxrs_ti)))
         
            axs[1,0].plot(cxrs_rho, cxrs_vtor , markersize=4, fillstyle='none',
                            label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].cxrs_vtor_time[cxrs_vtor_time_ind_list[i]]), color=color_list[i], marker='')
            axs[1,0].fill_between(cxrs_rho, cxrs_vtor-cxrs_vtor_err, cxrs_vtor+cxrs_vtor_err, color=color_list[i], alpha=0.2)
            
            axs[1,1].plot(cxrs_rho, cxrs_ti/1e3,markersize=4, fillstyle='none',
                        label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].cxrs_time[cxrs_time_ind_list[i]]), color=color_list[i], marker='')
            axs[1,1].fill_between(cxrs_rho, (cxrs_ti-cxrs_ti_err)/1e3, (cxrs_ti+cxrs_ti_err)/1e3, color=color_list[i], alpha=0.2)
            
      
        axs[0,0].legend(fontsize="8")

    except:
        print('No CXRS data')
    

    
    
    #     ax_te.plot(list_obj_tcv[i].th_rho, list_obj_tcv[i].th_te[:,th_time_ind_list[i]],marker='+', label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].th_time[th_time_ind_list[i]]))
    #     ax_ne.plot(list_obj_tcv[i].th_rho, list_obj_tcv[i].th_ne[:,th_time_ind_list[i]]/10**19,marker='+', label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].th_time[th_time_ind_list[i]]))
    #     ax_ti.plot(list_obj_tcv[i].cxrs_rho, list_obj_tcv[i].cxrs_ti[cxrs_time_ind_list[i],:],marker='+', label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].cxrs_time[cxrs_time_ind_list[i]]))
    # ax_te.legend()
    # ax_ne.legend()
    # ax_ti.legend()
    return axs




def plot_profiles_comparison_single_column(shot_list, time_list, plot_err=None, rhomin=None, rhomax=None, tavg=None):
    
    list_obj_tcv = [TCVShot(shot) for shot in shot_list]
    
    for i in range(len(list_obj_tcv)):
        list_obj_tcv[i].get_thomson_fit()
        list_obj_tcv[i].get_cxrs_fit(set_nan_to_zero=True)
        
    if len(shot_list) != len(time_list):
        print('\n --- choosing same time for all shots --- ')
        th_time_ind_list = [get_closest_ind(list_obj_tcv[i].th_time, time_list[0]) for i in range(len(list_obj_tcv))]
    
    th_time_ind_list = [get_closest_ind(list_obj_tcv[i].th_time, time_list[i]) for i in range(len(list_obj_tcv))]
    

    # color_list = ['blue', 'red', 'green', 'black', 'orange', 'purple', 'brown', 'pink', 'grey', 'cyan']
    # marker_list = 


    fig ,axs = prep_multiple_subplots(4,1, axgrid=[0,1,2,3], sharex=True)
    fig.suptitle('#{} ; tavg = {}'.format(shot_list, tavg))
    
    # fig, ax_te = plot_1d([], [], grid=True, xlabel=r'$\rho$', ylabel=r'$T_e$ $[eV]$')
    # fig, ax_ne = plot_1d([], [], grid=True, xlabel=r'$\rho$', ylabel=r'$n_e$ $[10^{19} m^{-3}]$')
    # fig, ax_ti = plot_1d([], [], grid=True, xlabel=r'$\rho$', ylabel=r'$T_i$ $[eV]$')
    for i in range(len(list_obj_tcv)):
        
        # ---rho range------------------------------------- #
        th_rho_all = list_obj_tcv[i].th_rho
        if rhomin is not None:
            rho_min_ind = get_closest_ind(th_rho_all, rhomin)
        else:
            rho_min_ind = None
        if rhomax is not None:
            rho_max_ind = get_closest_ind(th_rho_all, rhomax)+1
        else:
            rho_max_ind = None
        th_rho = th_rho_all[rho_min_ind:rho_max_ind]
        # ------------------------------------------------ #
        
        # ---time average--------------------------------- #
        if tavg is not None:
            tmin_ind = get_closest_ind(list_obj_tcv[i].th_time, list_obj_tcv[i].th_time[th_time_ind_list[i]]-tavg/2)
            tmax_ind = get_closest_ind(list_obj_tcv[i].th_time, list_obj_tcv[i].th_time[th_time_ind_list[i]]+tavg/2)
            # print(tmin_ind, tmax_ind)
            th_te = np.mean(list_obj_tcv[i].th_te[rho_min_ind:rho_max_ind, tmin_ind:tmax_ind], axis=1)
            th_ne = np.mean(list_obj_tcv[i].th_ne[rho_min_ind:rho_max_ind, tmin_ind:tmax_ind], axis=1)
            
            if plot_err:
                th_te_err = np.mean(list_obj_tcv[i].th_te_err[rho_min_ind:rho_max_ind, tmin_ind:tmax_ind], axis=1)
                th_ne_err = np.mean(list_obj_tcv[i].th_ne_err[rho_min_ind:rho_max_ind, tmin_ind:tmax_ind], axis=1)    
            else:
                th_te_err = np.zeros((np.shape(th_te)))
                th_ne_err = np.zeros((np.shape(th_ne)))
        # ------------------------------------------------ #
            
        else:
            th_te = list_obj_tcv[i].th_te[rho_min_ind:rho_max_ind,th_time_ind_list[i]]
            th_ne = list_obj_tcv[i].th_ne[rho_min_ind:rho_max_ind,th_time_ind_list[i]]
            if plot_err:
                th_te_err = list_obj_tcv[i].th_te_err[rho_min_ind:rho_max_ind,th_time_ind_list[i]]
                th_ne_err = list_obj_tcv[i].th_ne_err[rho_min_ind:rho_max_ind,th_time_ind_list[i]]
            else:
                th_te_err = np.zeros((np.shape(th_te)))
                th_ne_err = np.zeros((np.shape(th_ne)))
        
        # ---plotting thomson data------------------------ #                   
        # list_obj_tcv[i].plot_heating(ax=axs[0,0], color=color_list[i])
        axs[1].plot(th_rho, th_te/1e3 ,linewidth=2, markersize=4, fillstyle='none',
                      label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].th_time[th_time_ind_list[i]]),  color=color_list[i], marker='')
        axs[1].fill_between(th_rho, (th_te-th_te_err)/1e3, (th_te+th_te_err)/1e3, color=color_list[i], alpha=0.2)    
    
        axs[0].plot(th_rho, th_ne/10**19,linewidth=2, markersize=4, fillstyle='none',
                      label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].th_time[th_time_ind_list[i]]), color=color_list[i], marker='')
        axs[0].fill_between(th_rho, (th_ne-th_ne_err)/10**19, (th_ne+th_ne_err)/10**19, color=color_list[i], alpha=0.2)
        # ------------------------------------------------ #

    axs[0].legend(fontsize="8")
    # axs[0,1].legend()
    # axs[1,1].legend()
    axs[1].set_ylabel(r'$T_e$ $[keV]$')
    axs[0].set_ylabel(r'$n_e$ $[10^{19} m^{-3}]$')
    axs[3].set_ylabel(r'$v_{tor}$ $[km/s]$')
    axs[2].set_ylabel(r'$T_i$ $[keV]$')
    axs[3].set_xlabel(r'$\rho$')
      
    
    #Data CXRS 
    try:
        if len(shot_list) != len(time_list):
            print('\n --- choosing same time for all shots --- ')
            cxrs_time_ind_list = [get_closest_ind(list_obj_tcv[i].cxrs_time, time_list[0]) for i in range(len(list_obj_tcv))]
             
        cxrs_time_ind_list = [get_closest_ind(list_obj_tcv[i].cxrs_time, time_list[i]) for i in range(len(list_obj_tcv))]
        cxrs_vtor_time_ind_list = [get_closest_ind(list_obj_tcv[i].cxrs_vtor_time, time_list[i]) for i in range(len(list_obj_tcv))]
        
        for i in range(len(list_obj_tcv)):
            
            cxrs_rho_all = list_obj_tcv[i].cxrs_rho
            if rhomin is not None:
                rho_min_ind = get_closest_ind(cxrs_rho_all, rhomin)
            else:
                rho_min_ind = None
            if rhomax is not None:
                rho_max_ind = get_closest_ind(cxrs_rho_all, rhomax)+1
            else:
                rho_max_ind = None
            
            cxrs_rho = cxrs_rho_all[rho_min_ind:rho_max_ind]
            
            if tavg is not None:
                tmin_ind = get_closest_ind(list_obj_tcv[i].cxrs_vtor_time, list_obj_tcv[i].cxrs_vtor_time[cxrs_vtor_time_ind_list[i]]-tavg/2)
                tmax_ind = get_closest_ind(list_obj_tcv[i].cxrs_vtor_time, list_obj_tcv[i].cxrs_vtor_time[cxrs_vtor_time_ind_list[i]]+tavg/2)
                # print(tmin_ind, tmax_ind)*
                
                cxrs_vtor = np.mean(list_obj_tcv[i].cxrs_vtor[tmin_ind:tmax_ind, rho_min_ind:rho_max_ind], axis=0)
                
                tmin_ind = get_closest_ind(list_obj_tcv[i].cxrs_time, list_obj_tcv[i].cxrs_time[cxrs_time_ind_list[i]]-tavg/2)
                tmax_ind = get_closest_ind(list_obj_tcv[i].cxrs_time, list_obj_tcv[i].cxrs_time[cxrs_time_ind_list[i]]+tavg/2)
               
                cxrs_ti = np.mean(list_obj_tcv[i].cxrs_ti[tmin_ind:tmax_ind, rho_min_ind:rho_max_ind], axis=0)
                
                if plot_err:
                    cxrs_vtor_err = np.mean(list_obj_tcv[i].cxrs_vtor_err[tmin_ind:tmax_ind, rho_min_ind:rho_max_ind], axis=0)
                    cxrs_ti_err = np.mean(list_obj_tcv[i].cxrs_ti_err[tmin_ind:tmax_ind, rho_min_ind:rho_max_ind], axis=0)
                else:
                    cxrs_vtor_err = np.zeros((np.shape(cxrs_vtor)))
                    cxrs_ti_err = np.zeros((np.shape(cxrs_ti)))
                
            else:
                cxrs_vtor = list_obj_tcv[i].cxrs_vtor[cxrs_vtor_time_ind_list[i], rho_min_ind:rho_max_ind]
                cxrs_ti = list_obj_tcv[i].cxrs_ti[cxrs_time_ind_list[i], rho_min_ind:rho_max_ind]

                if plot_err:
                    cxrs_vtor_err = list_obj_tcv[i].cxrs_vtor_err[cxrs_vtor_time_ind_list[i], rho_min_ind:rho_max_ind]
                    cxrs_ti_err = list_obj_tcv[i].cxrs_ti_err[cxrs_time_ind_list[i], rho_min_ind:rho_max_ind]
                else:
                    cxrs_vtor_err = np.zeros((np.shape(cxrs_vtor)))
                    cxrs_ti_err = np.zeros((np.shape(cxrs_ti)))
         
            axs[3].plot(cxrs_rho, cxrs_vtor , markersize=4, fillstyle='none',
                            label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].cxrs_vtor_time[cxrs_vtor_time_ind_list[i]]), color=color_list[i], marker='')
            axs[3].fill_between(cxrs_rho, cxrs_vtor-cxrs_vtor_err, cxrs_vtor+cxrs_vtor_err, color=color_list[i], alpha=0.2)
            
            axs[2].plot(cxrs_rho, cxrs_ti/1e3,markersize=4, fillstyle='none',
                        label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].cxrs_time[cxrs_time_ind_list[i]]), color=color_list[i], marker='')
            axs[2].fill_between(cxrs_rho, (cxrs_ti-cxrs_ti_err)/1e3, (cxrs_ti+cxrs_ti_err)/1e3, color=color_list[i], alpha=0.2)
            
      
        axs[0].legend(fontsize="8")

    except:
        print('No CXRS data')
    

    
    
    #     ax_te.plot(list_obj_tcv[i].th_rho, list_obj_tcv[i].th_te[:,th_time_ind_list[i]],marker='+', label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].th_time[th_time_ind_list[i]]))
    #     ax_ne.plot(list_obj_tcv[i].th_rho, list_obj_tcv[i].th_ne[:,th_time_ind_list[i]]/10**19,marker='+', label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].th_time[th_time_ind_list[i]]))
    #     ax_ti.plot(list_obj_tcv[i].cxrs_rho, list_obj_tcv[i].cxrs_ti[cxrs_time_ind_list[i],:],marker='+', label='#{}, T={:.3f}'.format(shot_list[i], list_obj_tcv[i].cxrs_time[cxrs_time_ind_list[i]]))
    # ax_te.legend()
    # ax_ne.legend()
    # ax_ti.legend()
    return axs






            
#%% Function to store

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


