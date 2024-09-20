#%%
#=== general info ===#
#Author: Olivier Panico
#contact: olivier.panico@free.fr

#Goal: performs correlation analysis for tcv data

#=== imports ===#
#General imports
import numpy as np
import matplotlib.pyplot as plt

#Local imports (better to avoid)
from diag_tcv.dbsAnalysis.correlationAnalysis.modeCorrelation import isModeCorrelation
from diag_tcv.dbsAnalysis.correlationAnalysis.getRawData import get_raw_signals
from diag_tcv.dbsAnalysis.correlationAnalysis.correlationFunctions import full_coherence_analysis
import diag_tcv.utils.myMplStyle 

#DBS
from DBS.io.read import get_DBS_params #Channels params (dt, freq, ...)
from DBS import definitions as defs #Location of data
from DBS.io.interface import DataInterface #General class for DBS data
from DBS.processing.sigprocessing import remove_signal_edges, get_normalized_complex_signal
from DBS.beamtracing.interface import Beam3dInterface #raytracing
from DBS.beamtracing.DBSbeam import _DBSbeam
#dataAnalysis
from dataAnalysis.utils.plot_utils import plot_1d, plot_2d, my_legend, my_text, prep_multiple_subplots
from dataAnalysis.spectral.spectralAnalysis import custom_coherence, custom_time_coherence
from dataAnalysis.utils.utils import find_plateaus, my_linearRegression
#TCV
from diag_tcv.shotAnalysis.dischargeInfoMdsObject import TCVShot
from diag_tcv.dbsAnalysis.velocityProfile.velocityFunctions import get_velocity_prof_object

### ========== ###
### PARAMETERS ###
### ========== ###
# cmap='viridis'
# cmap='jet'
cmap='Dark2'
markers = ['o', '^', 's', 'D', 'x', 'v', '<', '>', 'p', 'h']
### PHYSICAL CONSTANTS ###


### ================= ###
### Utility functions ###
### ================= ###


### ============== ###
### Main functions ###
### ============== ###
class CorrelationAnalysis(TCVShot):
    '''
    Class is inherited from TCVShot so that the whole shot analysis is available
    '''    
    
    ### ============================ ###
    ### PARAMS OF DBS FOR WHOLE SHOT ###
    ### ============================ ###
    def __init__(self, shot, xmode=1, numDemod=False, machine='tcv', verbose=False, plot=False):
        super().__init__(shot=shot, verbose=verbose) #init with master class TCVShot
        
        self.machine = machine
        self.tmp_dir = defs.DATA_TMP_DIR / f'{self.machine}'

        # self.gas = isHydrogen(self.shot) moved to dischargeObject
        
        self.modeCorrelation = isModeCorrelation(self.shot)
        self.numDemod = numDemod

        
        self.xmode = xmode    
        self.plot = plot
    
        self._load_sweep_params()
        
        self.processedData = dict()
        for i in range(self.nbsweep):
            self.processedData['sweep'+str(i+1)] = dict()

        ### Identificiation of the plateaus on reference frequencies ###
        plateaus_indices = find_plateaus(self.params_ref.F)
        nb_plateaus = len(plateaus_indices)
        print(' --> {} plateaus '.format(nb_plateaus))
        self.plateaus_indices = plateaus_indices
        self.nb_plateaus = nb_plateaus


    def _load_sweep_params(self):
        if self.numDemod:
            self.params_ref = DataInterface(self.shot, isweep=False, channelval=3, machine='tcv').params
            self.params_hop = DataInterface(self.shot, isweep=False, channelval=4, machine='tcv').params
        else:
            self.params_ref = DataInterface(self.shot, isweep=False, channelval=1, machine='tcv').params
            self.params_hop = DataInterface(self.shot, isweep=False, channelval=2, machine='tcv').params

        print(' \n === #{} dbs params === '.format(self.shot))
        
        #Nb sweep and associated times
        sweep_tinit = self.params_ref.TDIFDOP + self.params_ref.t0seq
        nbsweep = len(sweep_tinit)
        self.nbsweep=nbsweep
        period = self.params_ref.Period
        print(' --> {} sweeps '.format(nbsweep))
        print('     Associated time windows : ')
        for i in range(nbsweep):
            print('     sweep {} : [{:.2f}, {:.2f}] s'.format(i+1, sweep_tinit[i], sweep_tinit[i]+period*1e-3))

        #Nbpts per freq
        assert(self.params_ref.dtAcq==self.params_hop.dtAcq)
        nbpts_per_freq_ref = int(self.params_ref.dtStep*1e-3/self.params_ref.dtAcq)
        nbpts_per_freq_hop = int(self.params_hop.dtStep*1e-3/self.params_hop.dtAcq)
        if nbpts_per_freq_ref==nbpts_per_freq_hop:
            self.same_nbpts_per_freq = True
        else:
            self.same_nbpts_per_freq = False
        print(' --> same nbpts per freq step: {}'.format(self.same_nbpts_per_freq))   
        
        #Mode correlation: ref and hop have artificially inverted doppler shift
        #If mode correlation => then correlate ref and np.conjugate(hop)
        print(' --> mode correlation: {}'.format(self.modeCorrelation))
        
        print(' --> params stored in dbs_ref.params / dbs_hop.params')

    ### ========================== ###
    ### ACCESS TO DATA FOR A SWEEP ###
    ### ========================== ###
    def load_from_sweep(self, isweep=1):  
        '''
        Load DataInterface from specific sweep
        '''
        if self.numDemod:
            data_ref = DataInterface(shot=self.shot, isweep=isweep, channelval=3, machine='tcv', verbose=self.verbose)
            data_hop = DataInterface(shot=self.shot, isweep=isweep, channelval=4, machine='tcv', verbose=self.verbose)
        else: 
            data_ref = DataInterface(shot=self.shot, isweep=isweep, channelval=1, machine='tcv', verbose=self.verbose)
            data_hop = DataInterface(shot=self.shot, isweep=isweep, channelval=2, machine='tcv', verbose=self.verbose)

        self.processedData['sweep'+str(isweep)]['data_ref'] = data_ref
        self.processedData['sweep'+str(isweep)]['data_hop'] = data_hop


    def load_from_time(self, time):
        '''
        Load DataInterface from time => will choose the right sweep
        '''
        if self.numDemod:
            data_ref = DataInterface.from_time(shot=self.shot, time=time, channelval=3, machine='tcv', verbose=self.verbose)
            data_hop = DataInterface.from_time(shot=self.shot, time=time, channelval=4, machine='tcv', verbose=self.verbose)
        else:
            data_ref = DataInterface.from_time(shot=self.shot, time=time, channelval=1, machine='tcv', verbose=self.verbose)
            data_hop = DataInterface.from_time(shot=self.shot, time=time, channelval=2, machine='tcv', verbose=self.verbose)
        
        isweep = data_ref.isweep
        self.processedData['sweep'+str(isweep)]['data_ref'] = data_ref
        self.processedData['sweep'+str(isweep)]['data_hop'] = data_hop
        
        
    def plot_freq_pattern_sweep(self, isweep=1, title=None, returnfig=False):
        '''
        Plot the pattern used for the given sweep
        '''
        
        t0F_structured_ref = self.params_ref.t0F_structured[isweep-1,:]
        freq_ref = self.params_ref.F
        
        t0F_structured_hop = self.params_hop.t0F_structured[isweep-1,:]
        freq_hop = self.params_hop.F
        
        fig, ax = plot_1d([],[], grid=True)
        if title is None:
            ax.set_title(' #{} sweep {}'.format(self.shot, isweep))
        else:
            ax.set_title(title)
            
        if self.same_nbpts_per_freq==True:
            ax.plot(t0F_structured_ref, freq_ref, color='blue', marker='o',linestyle='', label='reference channel')
            ax.plot(t0F_structured_hop, freq_hop, color='red', marker='^',linestyle='', label='hopping channel')
        else:
            ax.scatter(t0F_structured_ref, freq_ref, color='blue', marker='o',linestyle='', label='reference channel')
            ax.plot(t0F_structured_hop, freq_hop, color='red', marker='^',linestyle='', label='hopping channel')
            arrow_length = np.diff(data_ref.params.t0F)
            
            for i in range(len(freq_ref)):
                ax.arrow(t0F_structured_ref[i], freq_ref[i], arrow_length[i], 0,zorder=5, linewidth=1.5, head_width=0, head_length=0, ec='blue')
        
        my_legend(ax)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Frequency [GHz]')
        
        if returnfig:
            return fig, ax


        
    def plot_heating_with_sweeps(self, list_isweep, **kwargs):
        
        if 'markevery' not in kwargs:
            kwargs['markevery'] = 20
        
        ax = self.plot_heating(**kwargs)
        
        #Nb sweep and associated times
        sweep_tinit = self.params_ref.TDIFDOP + self.params_ref.t0seq
        nbsweep = len(sweep_tinit)
        self.nbsweep=nbsweep
        period = self.params_ref.Period
        bbox=dict(facecolor = 'white', edgecolor='black', pad=10.0)
        for i in range(nbsweep):
            if i+1 in list_isweep:
                ax.axvline(sweep_tinit[i], color='black')
                ax.axvline(sweep_tinit[i]+period*1e-3, color='black')
                # my_text(ax, sweep_tinit[i], 0.2, '{}'.format(i+1), color='black')
                ax.text(sweep_tinit[i]+0.1,
                        100,
                        '{}'.format(i+1), 
                        horizontalalignment='center',
                        verticalalignment='center',
                        bbox=bbox,
                        color='black')
           

    ### ========================== ###
    ### ACCESS TO DATA FOR A SWEEP ###
    ### ========================== ###
    def _get_raw_data_isweep(self, isweep, ifreq_list=None, time_window=None):
        '''
        load full isweep data for correlation
        if time_window is not False => then load the time window
        '''
        if self.verbose:
            print('Getting raw data for sweep {}'.format(isweep))
        
        if ifreq_list is None:
            freq_ref = self.params_ref.F
            freq_hop = self.params_hop.F
            nbfreq = len(freq_hop)
            ifreq_list=np.linspace(0,nbfreq-1, nbfreq, dtype=int)
            if self.verbose:
                print('ifreq list for isweep {} : {}'.format(isweep, ifreq_list))
            
        # self.isweep = isweep            
        # self.rawSigDic = get_raw_signals(self.shot, isweep, ifreq_list, numDemod=self.numDemod)
        
        ### Identificiation of the plateaus on reference frequencies ###
        plateaus_indices = find_plateaus(self.params_ref.F)
        nb_plateaus = len(plateaus_indices)
        
        self.processedData['sweep'+str(isweep)]['ifreq_list'] = ifreq_list
        self.processedData['sweep'+str(isweep)]['nb_plateaus'] = nb_plateaus
        self.processedData['sweep'+str(isweep)]['plateaus_indices'] = plateaus_indices
        
        self.processedData['sweep'+str(isweep)] = get_raw_signals(self.shot, isweep, ifreq_list,self.processedData['sweep'+str(isweep)],  numDemod=self.numDemod)
    

    def get_normalized_data_isweep(self, isweep, dtsart=400.e-6, dtend=100.e-6, ret=False):
        '''
        From raw data gives the normalized complex
        '''
        if self.verbose:
            print('Getting normalized data for sweep {}'.format(isweep))
            
        if 'data_ref' not in self.processedData['sweep'+str(isweep)]:
            self._get_raw_data_isweep(isweep, ifreq_list=None)
        
                
        ###LOADING RAW DATA###
        t = self.processedData['sweep'+str(isweep)]['t']
        freq_list_hop = self.processedData['sweep'+str(isweep)]['freq_list_hop']
        I_ref_list = self.processedData['sweep'+str(isweep)]['I_list_ref']
        Q_ref_list = self.processedData['sweep'+str(isweep)]['Q_list_ref']
        I_hop_list = self.processedData['sweep'+str(isweep)]['I_list_hop']
        Q_hop_list = self.processedData['sweep'+str(isweep)]['Q_list_hop']
        
        ###REMOVING EDGE OF FIRST FREQUENCY###
        t_reduced_0, _ = remove_signal_edges(t[0,:], I_ref_list[0,:], dtstart=dtsart, dtend=dtend, axis=-1)
       
        ###FILLING NEW ARRAYS### 
        z_list_ref = np.zeros((len(freq_list_hop), len(t_reduced_0)), dtype='complex')
        z_list_hop = np.zeros((len(freq_list_hop), len(t_reduced_0)), dtype='complex')
        t_reduced_list_ref = np.zeros((len(freq_list_hop), len(t_reduced_0)))
        t_reduced_list_hop = np.zeros((len(freq_list_hop), len(t_reduced_0)))
                                 
        for i in range(len(freq_list_hop)):
            t_reduced_ref_loc,I_ref_loc = remove_signal_edges(t[i,:], I_ref_list[i,:], dtstart=dtsart, dtend=dtend, axis=-1)
            _,Q_ref_loc = remove_signal_edges(t[i,:], Q_ref_list[i,:], dtstart=dtsart, dtend=dtend, axis=-1)
    
            z_ref_loc = get_normalized_complex_signal(I_ref_loc, Q_ref_loc, self.processedData['sweep'+str(isweep)]['params_ref'].phase_cor)
            
            t_reduced_hop_loc,I_hop_loc = remove_signal_edges(t[i,:], I_hop_list[i,:], dtstart=dtsart, dtend=dtend, axis=-1)
            _,Q_hop_loc = remove_signal_edges(t[i,:], Q_hop_list[i,:], dtstart=dtsart, dtend=dtend, axis=-1)
    
            z_hop_loc = get_normalized_complex_signal(I_hop_loc, Q_hop_loc, self.processedData['sweep'+str(isweep)]['params_hop'].phase_cor)
        
            z_list_ref[i,:] = z_ref_loc
            z_list_hop[i,:] = z_hop_loc
            
            t_reduced_list_ref[i,:] = t_reduced_ref_loc
            t_reduced_list_hop[i,:] = t_reduced_hop_loc 
        
        ### The data are stored for the whole sweep ###
        self.processedData['sweep'+str(isweep)]['dt'] = t[0,1] - t[0,0]
        self.processedData['sweep'+str(isweep)]['t_reduced_list_ref'] = t_reduced_list_ref
        self.processedData['sweep'+str(isweep)]['t_reduced_list_hop'] = t_reduced_list_hop
        self.processedData['sweep'+str(isweep)]['z_list_ref'] = z_list_ref
        self.processedData['sweep'+str(isweep)]['z_list_hop'] = z_list_hop
        
        if ret:
            return z_list_ref, z_list_hop, t_reduced_list_ref, t_reduced_list_hop
        
    def get_raytracing_isweep(self, isweep, retdata=False):
        '''
        Load the raytracing for the whole sweep
        '''
        

        # if self.numDemod:
        #     output_ref, beam3d_interface_ref = _DBSbeam('tcv',shot=self.shot, isweep=isweep, xmode=self.xmode, channelval=3, ifreqs='all', verbose=self.verbose, plot=self.plot, load_if_existing=True, nrays_radial=1, nrays_azimuthal=1)
        #     output_hop, beam3d_interface_hop = _DBSbeam('tcv',shot=self.shot, isweep=isweep, xmode=self.xmode, channelval=4, ifreqs='all', verbose=self.verbose, plot=self.plot, load_if_existing=True, nrays_radial=1, nrays_azimuthal=1)
        # else:
        output_ref, beam3d_interface_ref = _DBSbeam('tcv',shot=self.shot, isweep=isweep, xmode=self.xmode, channelval=1, ifreqs='all', verbose=self.verbose, plot=self.plot, load_if_existing=True, nrays_radial=1, nrays_azimuthal=1)
        output_hop, beam3d_interface_hop = _DBSbeam('tcv',shot=self.shot, isweep=isweep, xmode=self.xmode, channelval=2, ifreqs='all', verbose=self.verbose, plot=self.plot, load_if_existing=True, nrays_radial=1, nrays_azimuthal=1)

    

        # outp = interface.fetch_result()

        self.processedData['sweep'+str(isweep)]['output_ref'] = output_ref
        self.processedData['sweep'+str(isweep)]['output_hop'] = output_hop
        self.processedData['sweep'+str(isweep)]['beam3d_interface_ref'] = beam3d_interface_ref
        self.processedData['sweep'+str(isweep)]['beam3d_interface_hop'] = beam3d_interface_hop
        self.processedData['sweep'+str(isweep)]['rho_list_ref'] = output_ref.rho
        self.processedData['sweep'+str(isweep)]['rho_list_hop'] = output_hop.rho
        
        if retdata:
            return output_ref.rho, output_hop.rho
    
    
    def get_delta(self, isweep, retdata=False):
        '''
        Get the delta between the two beams
        '''
        
        if 'rho_list_ref' not in self.processedData['sweep'+str(isweep)]:
                self.get_raytracing_isweep(isweep)
                
        nbfreq = len(self.processedData['sweep'+str(isweep)]['rho_list_ref'])
            
        r_ref = np.zeros((nbfreq))
        z_ref = np.zeros((nbfreq))
        r_hop = np.zeros((nbfreq))
        z_hop = np.zeros((nbfreq))
        delta = np.zeros((nbfreq))
        
        for i in range(nbfreq):
            if self.verbose:
                print('get delta, ifreq: ', i)
                print(self.processedData['sweep'+str(isweep)]['output_ref'])
            r_ref[i] = self.processedData['sweep'+str(isweep)]['output_ref'].dif[i].x
            z_ref[i] = self.processedData['sweep'+str(isweep)]['output_ref'].dif[i].z

            r_hop[i] =self.processedData['sweep'+str(isweep)]['output_hop'].dif[i].x
            z_hop[i] = self.processedData['sweep'+str(isweep)]['output_hop'].dif[i].z

            if r_ref[i]-r_hop[i]>0:
                delta[i] = -np.sqrt((r_ref[i]-r_hop[i])**2 + (z_ref[i]-z_hop[i])**2)
            else:
                delta[i] = np.sqrt((r_ref[i]-r_hop[i])**2 + (z_ref[i]-z_hop[i])**2)
    
        self.processedData['sweep'+str(isweep)]['delta'] = delta
        
        if retdata:
            return delta
    
    
    def get_rho_s_plateau(self, isweep, retdata=False):
        '''
        Get the normalization rho_s for each plateau in the sweep
        '''
        if 'delta' not in self.processedData['sweep'+str(isweep)]:
            self.get_delta(isweep)
        
        nbfreq = len(self.processedData['sweep'+str(isweep)]['delta'])
        r_ref = np.zeros((nbfreq))
        z_ref = np.zeros((nbfreq))
        
        sweep_tinit  = self.params_ref.TDIFDOP + self.params_ref.t0seq[isweep-1]
        length_sweep = self.params_ref.Period*1e-3
        
        for i in range(nbfreq):
            r_ref[i] = self.processedData['sweep'+str(isweep)]['output_ref'].dif[i].x
            z_ref[i] = self.processedData['sweep'+str(isweep)]['output_ref'].dif[i].z

        
        rho_s_plateaus = np.zeros((self.nb_plateaus))
        for i in range(self.nb_plateaus):
            r_ref_loc = r_ref[self.plateaus_indices[i][0]:self.plateaus_indices[i][1]+1]
            z_ref_loc = z_ref[self.plateaus_indices[i][0]:self.plateaus_indices[i][1]+1]
            rho_s_loc = self.get_rho_s_r_z(time_window=[sweep_tinit, sweep_tinit+length_sweep], r=np.mean(r_ref_loc),
                                           z=np.mean(z_ref_loc), rho=np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][self.plateaus_indices[i][0]:self.plateaus_indices[i][1]+1]))
            #put the rho_s in cm
            rho_s_loc = rho_s_loc*100
            rho_s_plateaus[i] = rho_s_loc
            
        self.processedData['sweep'+str(isweep)]['rho_s_plateaus'] = rho_s_plateaus
        
        if retdata:
            return rho_s_plateaus
    
    
    def get_correlation_isweep(self, isweep, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode='amp', plot=False, retdata=False):
        """
        isweep : int, sweep number
        """
        
        if self.verbose:
            print('Performing correlation for sweep {}'.format(isweep))
        
        if 'z_list_hop' not in self.processedData['sweep'+str(isweep)]:
            print(' --- getting normalized data for sweep {}'.format(isweep) + ' --- ')
            self.get_normalized_data_isweep(isweep=isweep, dtsart=400.e-6, dtend=100.e-6)
    

        ifreq_list = np.arange(0, len(self.processedData['sweep'+str(isweep)]['freq_list_hop']))

        # arrays to fill
        corr_list = np.zeros((len(ifreq_list), nperseg), dtype=complex)
        maxcorr_list = np.zeros((len(ifreq_list)))
        
        spectral_coh_list = np.zeros((len(ifreq_list), nperseg), dtype=complex)
        raw_maxspectralcoh_list = np.zeros((len(ifreq_list)))
        fit_maxspectralcoh_list = np.zeros((len(ifreq_list)))
        err_fit_maxspectralcoh_list = np.zeros((len(ifreq_list)))
        
        scipy_corr_list = np.zeros((len(ifreq_list), nperseg), dtype=complex)
        scipy_maxcorr_list = np.zeros((len(ifreq_list)))
        
        dt=self.processedData['sweep'+str(isweep)]['dt']
        
        if plot:
            fig, axfullcoherence=plot_1d([], [], grid=True)
        else:
            axfullcoherence = None
            
        for i, ifreq in enumerate(ifreq_list):
            zref_loc = self.processedData['sweep'+str(isweep)]['z_list_ref'][ifreq]
            zhop_loc = self.processedData['sweep'+str(isweep)]['z_list_hop'][ifreq]

            # for each zref and zhop we compute the full coherence analysis
            dictFullCohAnalysis = full_coherence_analysis(zref_loc, zhop_loc, dt, nperseg=nperseg, noverlap=noverlap, window=window,
                                                        remove_mean=remove_mean, mode=mode, plot=plot,ax=axfullcoherence, verbose=self.verbose)
            
            tcorr_loc = dictFullCohAnalysis['tcorr_spec']
            corr_loc = dictFullCohAnalysis['corr']
            fcsd_spectral = dictFullCohAnalysis['fcsd']
            spectral_coh_loc = dictFullCohAnalysis['spectral_coh']
            tcorr_scipy = dictFullCohAnalysis['tcorr_scipy']
            corr_scipy = dictFullCohAnalysis['corr_scipy']
            
            maxcorr = dictFullCohAnalysis['max_corr']
            raw_maxspectralcoh = dictFullCohAnalysis['max_raw_spectral_coh']
            scipy_maxcorr = dictFullCohAnalysis['max_corr_scipy']
            fit_maxspectralcoh = dictFullCohAnalysis['max_fit_spectral_coh']
            err_fit_maxspectralcoh = dictFullCohAnalysis['err_max_fit_spectral_coh']

            #We fill with full functions and with the max values
            corr_list[i,:] = corr_loc
            spectral_coh_list[i,:] = spectral_coh_loc
            scipy_corr_list[i,:] = corr_scipy
            
            maxcorr_list[i] = maxcorr
            scipy_maxcorr_list[i] = scipy_maxcorr
            raw_maxspectralcoh_list[i] = raw_maxspectralcoh
            fit_maxspectralcoh_list[i] = fit_maxspectralcoh
            err_fit_maxspectralcoh_list[i] = err_fit_maxspectralcoh
            
        corr_list = np.array(np.real(corr_list))
        maxcorr_list = np.array(maxcorr_list)
        spectral_coh_list = np.array((spectral_coh_list))
        raw_maxspectralcoh_list = np.array(raw_maxspectralcoh_list)
        fit_maxspectralcoh_list = np.array(fit_maxspectralcoh_list)
        err_fit_maxspectralcoh_list = np.array(err_fit_maxspectralcoh_list)
        scipy_corr_list = np.array((scipy_corr_list))
        scipy_maxcorr_list = np.array(scipy_maxcorr_list)
        
        #We save the results in the processedData dictionary   
        self.processedData['sweep'+str(isweep)]['tcorr'] = tcorr_loc
        self.processedData['sweep'+str(isweep)]['corr_list'] = corr_list
        self.processedData['sweep'+str(isweep)]['maxcorr_list'] = maxcorr_list
        
        self.processedData['sweep'+str(isweep)]['fcsd_spectral'] = fcsd_spectral
        self.processedData['sweep'+str(isweep)]['spectral_coh_list'] = spectral_coh_list
        self.processedData['sweep'+str(isweep)]['raw_maxspectralcoh_list'] = raw_maxspectralcoh_list
        self.processedData['sweep'+str(isweep)]['fit_maxspectralcoh_list'] = fit_maxspectralcoh_list
        self.processedData['sweep'+str(isweep)]['err_fit_maxspectralcoh_list'] = err_fit_maxspectralcoh_list
        
        self.processedData['sweep'+str(isweep)]['scipy_corr_list'] = scipy_corr_list
        self.processedData['sweep'+str(isweep)]['scipy_maxcorr_list'] = scipy_maxcorr_list
        
        if retdata:
            return tcorr_loc, corr_list, maxcorr_list, fcsd_spectral, spectral_coh_list, raw_maxspectralcoh_list, fit_maxspectralcoh_list, err_fit_maxspectralcoh_list, scipy_corr_list, scipy_maxcorr_list

    def get_velocity_isweep(self, isweep, retdata=False):
        '''
            Get the velocity for a given sweep and save it in processedData
        '''
        velocity_prof = get_velocity_prof_object(self.shot, isweep_list = [isweep])
        self.processedData['sweep'+str(isweep)]['velocity_prof'] = velocity_prof
        if retdata:
            return velocity_prof



    ### ======================================================================================== ###
    
    def prepare_coherence_for_plot(self, isweep, add_pt_zero=True, x_norm='rho_s', mode='amp',  load_if_existing=False, retdata=False):
        """
        Prepare the coherence for a given sweep 
        
        x_norm : str, 
            - 'rho_s': delta given in rho_s units
            - 'rho_hop' : position of hop channel
            - 'cm': delta in cm
        
        mode: str,
            - 'amp' : amplitude
            - 'phase' : phase
            - 'full' : full signal    

        """
        
        nb_plateaus = self.nb_plateaus
        plateau_list = np.arange(nb_plateaus)

        ### if data not loaded for the sweep we load them
        if load_if_existing is False or 'maxcoh_list' not in self.processedData['sweep'+str(isweep)]:
            print(' --- getting coherence data for sweep {}'.format(isweep) + ' --- ')
            self.get_correlation_isweep(isweep, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False, retdata=False)
        ### delta is computed depending on the chosen normalization
        if load_if_existing is False or 'delta' not in self.processedData['sweep'+str(isweep)]:
            self.get_delta(isweep)
        if load_if_existing is False or 'rho_s_plateaus' not in self.processedData['sweep'+str(isweep)]:
            self.get_rho_s_plateau(isweep)
        

        rho_list_hop = self.processedData['sweep'+str(isweep)]['rho_list_hop']
        
        if x_norm == 'rho_hop':
            delta = rho_list_hop
        else:
            delta = self.processedData['sweep'+str(isweep)]['delta']
        
        rho_s_plateaus = self.processedData['sweep'+str(isweep)]['rho_s_plateaus']

        self.processedData['sweep'+str(isweep)]['prepData'] = dict()
        self.processedData['sweep'+str(isweep)]['prepData']['nb_plateaus'] = nb_plateaus
        self.processedData['sweep'+str(isweep)]['prepData']['plateau_list'] = plateau_list
        
        for i, plat in enumerate(plateau_list):
            
            #Step 0: choose the plateau
            plat_indices = self.plateaus_indices[plat]
            
            #Step 1: load delta for a given plateau
            delta_loc = delta[plat_indices[0]:plat_indices[1]+1]
            rho_loc = np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][plat_indices[0]:plat_indices[1]+1])
            rho_s_loc = rho_s_plateaus[plat]
            
            rho_hop_loc = rho_list_hop[plat_indices[0]:plat_indices[1]+1]
            rho_hop_err_plus = np.max(rho_hop_loc) - rho_loc
            rho_hop_err_minus = rho_loc - np.min(rho_hop_loc)
            #Step 2: if normalization by rho_s
            if x_norm=='rho_s':
                delta_loc = delta_loc/rho_s_loc
            
            
            #Step 3: load arrays for the given plateau
            maxcorr_list         = self.processedData['sweep'+str(isweep)]['maxcorr_list'][plat_indices[0]:plat_indices[1]+1]
            raw_maxspectralcoh_list = self.processedData['sweep'+str(isweep)]['raw_maxspectralcoh_list'][plat_indices[0]:plat_indices[1]+1]
            fit_maxspectralcoh_list = self.processedData['sweep'+str(isweep)]['fit_maxspectralcoh_list'][plat_indices[0]:plat_indices[1]+1]
            err_fit_maxspectralcoh_list = self.processedData['sweep'+str(isweep)]['err_fit_maxspectralcoh_list'][plat_indices[0]:plat_indices[1]+1]
            
            #Step 4: add the point zero
            if add_pt_zero:
                if x_norm == 'rho_hop':
                    delta_loc = np.concatenate(([rho_loc], delta_loc))
                    maxcorr_list = np.concatenate(([1], maxcorr_list))
                    raw_maxspectralcoh_list = np.concatenate(([1], raw_maxspectralcoh_list))
                    fit_maxspectralcoh_list = np.concatenate(([1], fit_maxspectralcoh_list))
                    err_fit_maxspectralcoh_list = np.concatenate(([0], err_fit_maxspectralcoh_list))
                else:
                    delta_loc = np.concatenate(([0], delta_loc))
                    maxcorr_list = np.concatenate(([1], maxcorr_list))
                    raw_maxspectralcoh_list = np.concatenate(([1], raw_maxspectralcoh_list))
                    fit_maxspectralcoh_list = np.concatenate(([1], fit_maxspectralcoh_list))
                    err_fit_maxspectralcoh_list = np.concatenate(([0], err_fit_maxspectralcoh_list))
                    
                #sorting the arrays along the delta components
                inds = delta_loc.argsort()
                delta_loc = delta_loc[inds]
                maxcorr_list = maxcorr_list[inds]
                raw_maxspectralcoh_list = raw_maxspectralcoh_list[inds]
                fit_maxspectralcoh_list = fit_maxspectralcoh_list[inds]
                err_fit_maxspectralcoh_list = err_fit_maxspectralcoh_list[inds]
            
            self.processedData['sweep'+str(isweep)]['prepData']['plateau'+str(plat)] = {'rho_loc':rho_loc, 'rho_s_loc':rho_s_loc, 'delta':delta_loc,
                                                                                        'rho_hop_err_plus':rho_hop_err_plus, 'rho_hop_err_minus':rho_hop_err_minus,
                                                                                        'maxcorr':maxcorr_list, 'raw_maxspectralcoh':raw_maxspectralcoh_list,
                                                                                        'fit_maxspectralcoh':fit_maxspectralcoh_list, 'err_fit_maxspectralcoh_list':err_fit_maxspectralcoh_list}
            

        if retdata:
            return self.processedData['sweep'+str(isweep)]['prepData']
    

    
    def plot_coherence_delta_isweep(self, isweep, plateau_list=None, ax=None, plot_fit_spec = True, plot_raw_spec=False, plot_corr=False, ylog=True, caption=True, add_pt_zero=True, x_norm='rho_s', mode='amp', load_if_existing=False, retdata=False, **kwargs):
        
        
        prepData=self.prepare_coherence_for_plot(isweep, add_pt_zero=add_pt_zero, x_norm=x_norm, mode=mode, load_if_existing=load_if_existing, retdata=True)
       
        nb_plateaus = prepData['nb_plateaus']
        plateau_list = prepData['plateau_list']
        if ax is None:
            create_ax = True
        else:
            create_ax = False
        default_colors = plt.cm.Dark2(np.linspace(0,1,len(plateau_list)+2))
        for i, plat in enumerate(plateau_list):

                rho_loc            = prepData['plateau'+str(plat)]['rho_loc']
                rho_s_loc          = prepData['plateau'+str(plat)]['rho_s_loc']
                delta_loc          = prepData['plateau'+str(plat)]['delta']
                maxcorr            = prepData['plateau'+str(plat)]['maxcorr']
                raw_maxspectralcoh = prepData['plateau'+str(plat)]['raw_maxspectralcoh']
                fit_maxspectralcoh = prepData['plateau'+str(plat)]['fit_maxspectralcoh']
                err_fit_maxspectralcoh = prepData['plateau'+str(plat)]['err_fit_maxspectralcoh_list']
                if create_ax:
                    fig, ax = plot_1d([], [], grid=True)
                    
                default_kwargs_fit_spec = {'label':'fit spec','marker':'o', 'color':default_colors[i]}
                default_kwargs_raw_spec = {'label':'raw spec','marker':'^', 'color':default_colors[i]}
                default_kwargs_corr = {'label':'corr','marker':'s', 'color':default_colors[i]}
                default_kwargs_fit_spec.update(kwargs)
                default_kwargs_raw_spec.update(kwargs)
                default_kwargs_corr.update(kwargs)
                
                if plot_fit_spec:
                    ax.errorbar(delta_loc, fit_maxspectralcoh,err_fit_maxspectralcoh, **default_kwargs_fit_spec)
                
                if plot_raw_spec:
                    ax.plot(delta_loc, raw_maxspectralcoh, linestyle='-.',  **default_kwargs_raw_spec)
                    
                if plot_corr:
                    ax.plot(delta_loc, maxcorr, fillstyle='none', linestyle='--', **default_kwargs_corr)
                    
                if caption:
                    ax.set_title(r' #{} sweep {} $\rho_\psi$ = {:.2f}'.format(self.shot, isweep, rho_loc))
                    my_legend(ax, loc='upper left')
                    my_text(ax, 0.8,0.2, r'$\rho_s$ = {:.1f} mm'.format(rho_s_loc*10))
                    ax.set_ylabel('correlation')
                    ax.set_ylim(0.1,1.1)
                    
                    if x_norm == 'rho_hop':
                        ax.set_xlabel(r'$\rho_{hop}$')
                        ax.axhline(1, color='black', linestyle='--')
                        ax.axvline(rho_loc, color='black', linestyle='--')
                    elif x_norm == 'delta':
                        ax.set_xlabel(r'$\Delta$ [cm]')
                        ax.axhline(1, color='black', linestyle='--')
                        ax.axvline(0, color='black', linestyle='--')
                    elif x_norm == 'rho_s':
                        ax.set_xlabel(r'$\Delta/\rho_s$')
                        ax.axhline(1, color='black', linestyle='--')
                        ax.axvline(0, color='black', linestyle='--')
                
                if ylog:
                    ax.set_yscale('log')

    
      
    
    def plot_velocity_isweep(self, isweep, plateau_list=None, ax=None, several_axes=False, x_norm='rho_hop', load_if_existing=True,caption=True, retdata=False, **kwargs):
        ''' plot the velocity for a given sweep and plateaus
            different options for x scale: 
                - rho hopping 
                - delta (cm)
                - delta (rho_s)
        '''

        if plateau_list is None:
            plateau_list = np.arange(self.nb_plateaus)

    

        default_colors = plt.cm.Dark2(np.linspace(0,1,len(plateau_list)+2))
            
        ### if data not loaded for the sweep we load them
        if load_if_existing is False or 'velocity_prof' not in self.processedData['sweep'+str(isweep)]:
            print(' --- getting velocity data for sweep {}'.format(isweep) + ' --- ')
            self.get_velocity_isweep(isweep)
        ### delta is computed depending on the chosen normalization
        if load_if_existing is False or 'delta' not in self.processedData['sweep'+str(isweep)]:
            self.get_delta(isweep)
        if load_if_existing is False or 'rho_s_plateaus' not in self.processedData['sweep'+str(isweep)]:
            self.get_rho_s_plateau(isweep)
        
        if x_norm == 'rho_hop':
            delta = self.processedData['sweep'+str(isweep)]['rho_list_hop']
        else:
            delta = self.processedData['sweep'+str(isweep)]['delta']
        
        rho_s_plateaus = self.processedData['sweep'+str(isweep)]['rho_s_plateaus']
        
        for i, plat in enumerate(plateau_list):
            
            #Step 0: choose the plateau
            plat_indices = self.plateaus_indices[plat]
            
            #Step 1: load delta for a given plateau
            delta_loc = delta[plat_indices[0]:plat_indices[1]+1]
            rho_loc = np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][plat_indices[0]:plat_indices[1]+1])
            rho_s_loc = rho_s_plateaus[plat]
            #Step 2: if normalization by rho_s
            if x_norm=='rho_s':
                delta_loc = delta_loc/rho_s_loc

        
            #Step 3: load arrays for the given plateau
            v_perp_loc = self.processedData['sweep'+str(isweep)]['velocity_prof'].v_perp.values[plat_indices[0]:plat_indices[1]+1]
            dv_perp_loc = self.processedData['sweep'+str(isweep)]['velocity_prof'].dv_perp[plat_indices[0]:plat_indices[1]+1]
            dv_perp_up_loc = self.processedData['sweep'+str(isweep)]['velocity_prof'].dv_perp_up[plat_indices[0]:plat_indices[1]+1]
            dv_perp_low_loc = self.processedData['sweep'+str(isweep)]['velocity_prof'].dv_perp_low[plat_indices[0]:plat_indices[1]+1]
            
            err = np.zeros((2, len(v_perp_loc)))
            err[1] = dv_perp_up_loc
            err[0] = dv_perp_low_loc
           
                
            #Step 5: plot the arrays
            if i==0 and ax is None:
                fig, ax = plot_1d([], [], grid=True)
            elif i>0 and several_axes:
                fig, ax = plot_1d([], [], grid=True)
                
            
            default_kwargs = {'marker':markers[i], 'color':default_colors[i]}
            default_kwargs.update(kwargs)
            
            ax.errorbar(delta_loc, v_perp_loc, err, **default_kwargs)
            
            if caption:
                ax.set_title(r' #{} sweep {} $\rho_\psi$ = {:.2f}'.format(self.shot, isweep, rho_loc))
                ax.set_ylabel(r'$v_\perp$ [m/s]')
                if x_norm == 'rho_hop':
                    ax.set_xlabel(r'$\rho_{hop}$')
                    ax.axvline(rho_loc, color='black', linestyle='--')
                elif x_norm == 'delta':
                    ax.set_xlabel(r'$\Delta$ [cm]')
                    ax.axvline(0, color='black', linestyle='--')
                elif x_norm == 'rho_s':
                    ax.set_xlabel(r'$\Delta/\rho_s$')
                    ax.axvline(0, color='black', linestyle='--')

            plt.tight_layout()
    
    

    def plot_correlation_velocity(self, isweep_list, plateau_list=None, mode='amp', x_norm='rho_hop', load_if_existing=True, retdata=False):
        '''
        Will create double figures: up the correlation and down the velocity
        '''

        if plateau_list is None:
            plateau_list = np.arange(self.nb_plateaus)
            
        colors = plt.cm.Dark2(np.linspace(0,1,len(isweep_list)+2))
        
        for i in range(len(plateau_list)):
            fig, axs = prep_multiple_subplots(2,1, axgrid=[0,1],figsize=(10,8), sharex=True)
            rho_list = []
            rho_s_list = []
            for j in range(len(isweep_list)):

                isweep_loc = isweep_list[j]
                
                if 'delta' not in self.processedData['sweep'+str(isweep_loc)]:
                    self.get_delta(isweep_loc)
                if 'rho_s_plateaus' not in self.processedData['sweep'+str(isweep_loc)]:
                    self.get_rho_s_plateau(isweep_loc)
                
                rho_loc = np.mean(self.processedData['sweep'+str(isweep_loc)]['rho_list_ref'][self.plateaus_indices[plateau_list[i]][0]:self.plateaus_indices[plateau_list[i]][1]+1])
                rho_s_loc = self.processedData['sweep'+str(isweep_loc)]['rho_s_plateaus'][plateau_list[i]]
                
                self.plot_coherence_delta_isweep(isweep_loc, plateau_list=[plateau_list[i]], ax = axs[0], plot_fit_spec = True, plot_raw_spec=False, plot_corr=False,
                                            ylog=True, add_pt_zero=True, x_norm=x_norm, mode=mode, load_if_existing=load_if_existing, caption=False,
                                            color=colors[j], marker=markers[j], linestyle='', label=r'isweep = {} ; $\rho$ = {:.2f} ; $\rho_s$ = {:.1f} mm'.format(isweep_loc, rho_loc, rho_s_loc*10))
                self.plot_velocity_isweep(isweep_loc, plateau_list=[plateau_list[i]], ax=axs[1], x_norm=x_norm,
                                          load_if_existing=load_if_existing, caption=False, color=colors[j], marker=markers[j], linestyle='')    
        
        
                rho_list.append(rho_loc)
                rho_s_list.append(rho_s_loc)
                                        
            axs[0].set_title(r'#{}'.format(self.shot))
        
            axs[0].set_ylabel('correlation')
            axs[0].set_ylim(0.1,1.1)
            axs[1].set_ylabel(r'$v_\perp$ [m/s]')
            if x_norm == 'rho_hop':
                axs[1].set_xlabel(r'$\rho_{hop}$')
            elif x_norm == 'delta':
                axs[1].set_xlabel(r'$\Delta$ [cm]')
            elif x_norm == 'rho_s':
                axs[1].set_xlabel(r'$\Delta/\rho_s$')
        
        
            # text_str = '\n'.join([r'$\rho_s=$' + f'{val*10:.1f}' for val in rho_s_list])
            # my_text(axs[0], 0.8,0.2, text_str)
            
            my_legend(axs[0], fontsize='12')
            
            plt.tight_layout()


   
    ### ================================== ###
    ### diagnostics using specific freqs ###
    ### ================================== ###        
    def get_correlation_isweep_plateau(self, isweep, ifreq_list=None, plateau=None, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode='amp', plot=False, ret=False):
        """
        isweep : int, sweep number
        ifreq_list : list of int, list of frequencies to be used 
        plateau: int, which frequency plateau to use
        """
        
        if self.verbose:
            print('Performing correlation for sweep {}'.format(isweep))
        
        if 'z_list_hop' not in self.processedData['sweep'+str(isweep)]:
            print(' --- getting normalized data for sweep {}'.format(isweep) + ' --- ')
            self.get_normalized_data_isweep(isweep=isweep, dtsart=400.e-6, dtend=100.e-6)
    
        #If no direct list of frequencies is given
        if ifreq_list is None:
            #We look if a plateau has been given
            if plateau is not None:
                ifreq_list = np.arange(self.processedData['sweep'+str(isweep)]['plateaus_indices'][plateau][0], self.processedData['sweep'+str(isweep)]['plateaus_indices'][plateau][1]+1)
                if self.verbose:
                    print(' --- using plateau {} : {} --- '.format(plateau, ifreq_list))
            #if both are none, then we use the full list
            else:
                ifreq_list = np.arange(0, len(self.processedData['sweep'+str(isweep)]['freq_list_hop']))
        else:
            if self.verbose:
                print(' --- using ifreq_list : {} --- '.format(ifreq_list))
            ifreq_list = np.array(ifreq_list)

        coh_list = np.zeros((len(ifreq_list), nperseg), dtype=complex)
        spectral_coh_list = np.zeros((len(ifreq_list), nperseg), dtype=complex)
        maxcoh_list = np.zeros((len(ifreq_list)), dtype=complex)
        maxspectralcoh_list = np.zeros((len(ifreq_list)), dtype=complex)
        
        dt=self.processedData['sweep'+str(isweep)]['dt']
        
        if plot:
            fig, axfullcoherence=plot_1d([], [], grid=True)
        else:
            axfullcoherence = None
            
        for i, ifreq in enumerate(ifreq_list):
            zref_loc = self.processedData['sweep'+str(isweep)]['z_list_ref'][ifreq]
            zhop_loc = self.processedData['sweep'+str(isweep)]['z_list_hop'][ifreq]

            tcorr_loc, coh_loc, fcsd_spectral, spectral_coh_loc = full_coherence_analysis(zref_loc, zhop_loc, dt, nperseg=nperseg, noverlap=noverlap, window=window,
                                                                                          remove_mean=remove_mean, mode=mode, plot=plot,ax=axfullcoherence, verbose=self.verbose)
            maxcoh = np.max(abs(coh_loc))
            maxspectralcoh = np.max(spectral_coh_loc)

            coh_list[i,:] = coh_loc
            spectral_coh_list[i,:] = spectral_coh_loc
            maxcoh_list[i] = maxcoh
            maxspectralcoh_list[i] = maxspectralcoh
            
            coh_list = np.array(np.real(coh_list))
            spectral_coh_list = np.array((spectral_coh_list))
            
            
        self.processedData['sweep'+str(isweep)]['tcorr'] = tcorr_loc
        self.processedData['sweep'+str(isweep)]['fcsd_spectral'] = fcsd_spectral
        self.processedData['sweep'+str(isweep)]['coh_list'] = coh_list
        self.processedData['sweep'+str(isweep)]['spectral_coh_list'] = spectral_coh_list
        self.processedData['sweep'+str(isweep)]['maxcoh_list'] = maxcoh_list
        self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'] = maxspectralcoh_list
        
        if ret:
            return tcorr_loc, coh_list, fcsd_spectral, spectral_coh_list, maxcoh_list, maxspectralcoh_list


    def plot_coherence_delta_isweep2(self, isweep, plateau_list=None, ax=None, plot_spectral=True, plot_pearson=False, ylog=True, caption=True, add_pt_zero=True, x_norm='rho_s', mode='amp', load_if_existing=False, retdata=False, **kwargs):
        """
        Plot the coherence for a given sweep 
        Inside this sweep you can choose to plot only a plateau 
        
        x_norm : str, 
            - 'rho_s': delta given in rho_s units
            - 'rho_hop' : position of hop channel
            - 'cm': delta in cm
        """
        
        nb_plateaus = self.nb_plateaus
        if plateau_list is None:
            plateau_list = np.arange(nb_plateaus)


        if ax is None:
            create_ax = True
        else:
            create_ax = False
            
        default_colors = plt.cm.Dark2(np.linspace(0,1,len(plateau_list)+2))
            
        ### if data not loaded for the sweep we load them
        if load_if_existing is False or 'maxcoh_list' not in self.processedData['sweep'+str(isweep)]:
            print(' --- getting coherence data for sweep {}'.format(isweep) + ' --- ')
            self.get_correlation_isweep(isweep, nperseg=1024, noverlap=512, window=None, remove_mean=True,mode=mode, plot=False, retdata=False)
        ### delta is computed depending on the chosen normalization
        if load_if_existing is False or 'delta' not in self.processedData['sweep'+str(isweep)]:
            self.get_delta(isweep)
        if load_if_existing is False or 'rho_s_plateaus' not in self.processedData['sweep'+str(isweep)]:
            self.get_rho_s_plateau(isweep)
        
        if x_norm == 'rho_hop':
            delta = self.processedData['sweep'+str(isweep)]['rho_list_hop']
        else:
            delta = self.processedData['sweep'+str(isweep)]['delta']
        
        rho_s_plateaus = self.processedData['sweep'+str(isweep)]['rho_s_plateaus']
        
        for i, plat in enumerate(plateau_list):
            
            #Step 0: choose the plateau
            plat_indices = self.plateaus_indices[plat]
            
            #Step 1: load delta for a given plateau
            delta_loc = delta[plat_indices[0]:plat_indices[1]+1]
            rho_loc = np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][plat_indices[0]:plat_indices[1]+1])
            rho_s_loc = rho_s_plateaus[plat]
            #Step 2: if normalization by rho_s
            if x_norm=='rho_s':
                delta_loc = delta_loc/rho_s_loc
           
           
            #Step 3: load arrays for the given plateau
            maxcoh_list         = self.processedData['sweep'+str(isweep)]['maxcoh_list'][plat_indices[0]:plat_indices[1]+1]
            maxspectralcoh_list = self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'][plat_indices[0]:plat_indices[1]+1]
            
            #Step 4: add the point zero
            if add_pt_zero:
                if x_norm == 'rho_hop':
                    delta_loc = np.concatenate(([rho_loc], delta_loc))
                    maxcoh_list = np.concatenate(([1], maxcoh_list))
                    maxspectralcoh_list = np.concatenate(([1], maxspectralcoh_list))
                else:
                    delta_loc = np.concatenate(([0], delta_loc))
                    maxcoh_list = np.concatenate(([1], maxcoh_list))
                    maxspectralcoh_list = np.concatenate(([1], maxspectralcoh_list))
                
                #sorting the arrays along the delta components
                inds = delta_loc.argsort()
                delta_loc = delta_loc[inds]
                maxcoh_list = maxcoh_list[inds]
                maxspectralcoh_list = maxspectralcoh_list[inds]
                
            #Step 5: plot the arrays
            if create_ax:
                fig, ax = plot_1d([], [], grid=True)
                
            default_kwargs_spec = {'label':'spectral','marker':'o', 'color':default_colors[i]}
            default_kwargs_pear = {'label':'pearson','marker':'o', 'color':default_colors[i]}
            default_kwargs_spec.update(kwargs)
            default_kwargs_pear.update(kwargs)
                
            if plot_spectral:
                ax.plot(delta_loc, maxspectralcoh_list,  **default_kwargs_spec)
            if plot_pearson:
                ax.plot(delta_loc, maxcoh_list, fillstyle='none', linestyle='--', **default_kwargs_pear)
           
            if caption:
                ax.set_title(r' #{} sweep {} $\rho_\psi$ = {:.2f}'.format(self.shot, isweep, rho_loc))
                my_legend(ax, loc='upper left')
                my_text(ax, 0.8,0.2, r'$\rho_s$ = {:.1f} mm'.format(rho_s_loc*10))
                ax.set_ylabel('correlation')
                ax.set_ylim(0.1,1.1)
                
                if x_norm == 'rho_hop':
                    ax.set_xlabel(r'$\rho_{hop}$')
                    ax.axhline(1, color='black', linestyle='--')
                    ax.axvline(rho_loc, color='black', linestyle='--')
                elif x_norm == 'delta':
                    ax.set_xlabel(r'$\Delta$ [cm]')
                    ax.axhline(1, color='black', linestyle='--')
                    ax.axvline(0, color='black', linestyle='--')
                elif x_norm == 'rho_s':
                    ax.set_xlabel(r'$\Delta/\rho_s$')
                    ax.axhline(1, color='black', linestyle='--')
                    ax.axvline(0, color='black', linestyle='--')
                
            if ylog:
                ax.set_yscale('log')


 
    def deprecated_plot_coherence_delta(self, isweep_list, ax_low_f=None, ax_high_f=None,  plot_spectral=True, plot_pearson=True,plot_rho=False, ylog=False, add_pt_zero=True, norm_rho_s=True, mode='full', retdata=False):
        
        print(' Careful this function is deprecated, use plot_coherence_delta_isweep instead')
        
        if ax_low_f is None:
            fig, ax_low_f = plot_1d([],[], grid=True)
        if ax_high_f is None: 
            fig, ax_high_f = plot_1d([],[], grid=True)
        
        sweep_tinit = self.params_ref.TDIFDOP + self.params_ref.t0seq
        period = self.params_ref.Period*1e-3
        rho_s_list_low_f = []
        rho_s_list_high_f = []   
        
        for i, isweep in enumerate(isweep_list):
            # self.get_correlation_isweep(isweep, mode=mode)  
            if 'maxcoh_list' not in self.processedData['sweep'+str(isweep)]:
                self.get_correlation_isweep(isweep, mode=mode)
            if 'rho_list_ref' not in self.processedData['sweep'+str(isweep)]:
                self.get_raytracing_isweep(isweep)
            
            
            r_ref = np.zeros((len(self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'])))
            z_ref = np.zeros((len(self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'])))
            r_hop = np.zeros((len(self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'])))
            z_hop = np.zeros((len(self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'])))
            delta = np.zeros((len(self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'])))
            
            for i in range(len(self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'])):
                r_ref[i] = self.processedData['sweep'+str(isweep)]['output_ref'].dif[i].x
                z_ref[i] = self.processedData['sweep'+str(isweep)]['output_ref'].dif[i].z
                
                r_hop[i] =self.processedData['sweep'+str(isweep)]['output_hop'].dif[i].x
                z_hop[i] = self.processedData['sweep'+str(isweep)]['output_hop'].dif[i].z

                if r_ref[i]-r_hop[i]>0:
                    delta[i] = -np.sqrt((r_ref[i]-r_hop[i])**2 + (z_ref[i]-z_hop[i])**2)
                else:
                    delta[i] = np.sqrt((r_ref[i]-r_hop[i])**2 + (z_ref[i]-z_hop[i])**2)
            
            
            if norm_rho_s:
                sweep_tinit_loc = sweep_tinit[isweep-1]
                period_loc = period
                
                rho_s_loc_low_f = self.get_rho_s_r_z([sweep_tinit_loc, sweep_tinit_loc+period_loc], np.mean(r_ref[:20]), np.mean(z_ref[:20]), np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20]))
                rho_s_loc_high_f = self.get_rho_s_r_z([sweep_tinit_loc, sweep_tinit_loc+period_loc], np.mean(r_ref[20:40]), np.mean(z_ref[20:40]), np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40]))
                
                #put the rho_s in cm
                rho_s_loc_low_f = rho_s_loc_low_f*100
                rho_s_loc_high_f = rho_s_loc_high_f*100
                
                rho_s_list_low_f.append(rho_s_loc_low_f)
                rho_s_list_high_f.append(rho_s_loc_high_f)
                
                delta[:20] = delta[:20]/rho_s_loc_low_f
                delta[20:40] = delta[20:40]/rho_s_loc_high_f
            
            
            if add_pt_zero:
                delta_low_f = np.concatenate(([0], delta[:20]))
                delta_high_f = np.concatenate(([0], delta[20:40]))
            
                val_low_f_spec = np.concatenate(([1], self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'][:20]))
                val_high_f_spec = np.concatenate(([1], self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'][20:40]))
                
                val_low_f_pear = np.concatenate(([1], self.processedData['sweep'+str(isweep)]['maxcoh_list'][:20]))
                val_high_f_pear = np.concatenate(([1], self.processedData['sweep'+str(isweep)]['maxcoh_list'][20:40]))
            
                #sorting the arrays along the delta components
                inds_low_f = delta_low_f.argsort()
                print('inds_low_f' , inds_low_f)
                inds_high_f = delta_high_f.argsort()
                delta_low_f = delta_low_f[inds_low_f]
                delta_high_f = delta_high_f[inds_high_f]
                val_low_f_spec = val_low_f_spec[inds_low_f]
                val_high_f_spec = val_high_f_spec[inds_high_f]
                val_low_f_pear = val_low_f_pear[inds_low_f]
                val_high_f_pear = val_high_f_pear[inds_high_f]
                
                
            
            else:
                delta_low_f = delta[:20]
                delta_high_f = delta[20:40]
            
                val_low_f_spec = self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'][:20]
                val_high_f_spec = self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'][20:40]
                
                val_low_f_pear = self.processedData['sweep'+str(isweep)]['maxcoh_list'][:20]
                val_high_f_pear = self.processedData['sweep'+str(isweep)]['maxcoh_list'][20:40]
                
            if plot_spectral:
                ax_low_f.plot(delta_low_f,
                              val_low_f_spec,
                              label=r'spec ; isweep {} ; $\rho$ = {:.2f}'.format(isweep, np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20])),
                              markersize=5,
                              linewidth=1)
                ax_high_f.plot(delta_high_f,
                               val_high_f_spec,
                               label=r'spec ; isweep {} ; $\rho$ = {:.2f}'.format(isweep, np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40])),
                                markersize=5,
                                linewidth=1)
            if plot_pearson:
                
                ax_low_f.plot(delta_low_f,
                              val_low_f_pear,
                                label=r'pear ; isweep {} ; $\rho$ = {:.2f}'.format(isweep, np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20])),
                                 markersize=5,
                              linewidth=1)
                ax_high_f.plot(delta_high_f,
                               val_high_f_pear,
                               label=r'pear ; isweep {} ; $\rho$ = {:.2f}'.format(isweep, np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40])),
                                markersize=5,
                              linewidth=1)

        
        ax_low_f.set_title(r'#{}'.format(self.shot))
        ax_high_f.set_title(r'#{}'.format(self.shot))
          
        my_legend(ax_low_f, fontsize="12")
        my_legend(ax_high_f, fontsize="12")
        
        if norm_rho_s:
            ax_low_f.set_xlabel(r'$\Delta/\rho_s$')
            ax_high_f.set_xlabel(r'$\Delta/\rho_s$')
        else:
            ax_low_f.set_xlabel(r'$\Delta$ $[cm]$')
            ax_high_f.set_xlabel(r'$\Delta$ $[cm]$')
            
        ax_high_f.set_ylabel(r'correlation')
        ax_low_f.set_ylabel(r'correlation')
        
        
        # ax_low_f.set_xlim(-8, 3)
        # ax_high_f.set_xlim(-8, 3)

        # ax_low_f.set_ylim(0.1, 1.1)
        # ax_high_f.set_ylim(0.1, 1.1)
        
        
        
        if ylog:
            ax_low_f.set_yscale('log')
            ax_high_f.set_yscale('log')

        

        
        if norm_rho_s and plot_rho:
            plot_1d(np.array(rho_s_list_low_f)*10, label='low f', grid=True)
            plt.plot(np.array(rho_s_list_high_f)*10, label='high f')
            plt.ylabel(r'$\rho_s$ [mm]')
            plt.xlabel('isweep')
            plt.legend()
        
        if retdata:
            data=dict()
            data['delta_low_f'] = delta_low_f
            data['val_low_f_spec'] = val_low_f_spec
            data['val_high_f_spec'] = val_high_f_spec
            data['val_low_f_pear'] = val_low_f_pear
            data['val_high_f_pear'] = val_high_f_pear
            data['rho_s_low_f'] = np.array(rho_s_list_low_f)
            data['rho_low_f'] = np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20])
            
            data['delta_high_f'] = delta_high_f
            data['val_low_f_spec'] = val_low_f_spec
            data['val_high_f_spec'] = val_high_f_spec
            data['val_low_f_pear'] = val_low_f_pear
            data['val_high_f_pear'] = val_high_f_pear
            data['rho_s_high_f'] = np.array(rho_s_list_high_f)
            data['rho_high_f'] = np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40])
            
            return  ax_low_f, ax_high_f, delta, data
        
        return ax_low_f, ax_high_f, delta
    
    
    
    def plot_coherence(self, isweep_list, ax_low_f=None, ax_high_f=None, delta_rho=True, plot_spectral=True, plot_pearson=True, ylog=False):

        if ax_low_f is None:
            fig, ax_low_f = plot_1d([],[], grid=True)
        if ax_high_f is None: 
            fig, ax_high_f = plot_1d([],[], grid=True)
            
        for i, isweep in enumerate(isweep_list):
            if 'maxcoh_list' not in self.processedData['sweep'+str(isweep)]:
                self.get_correlation_isweep(isweep)
            if 'rho_list_ref' not in self.processedData['sweep'+str(isweep)]:
                self.get_raytracing_isweep(isweep)
            
            
            
            if plot_spectral:
                if delta_rho:
                    ax_low_f.plot(self.processedData['sweep'+str(isweep)]['rho_list_hop'][:20]-self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20],
                                  self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'][:20],
                                  label=r'spec ; isweep {} ; $\rho$ = {:.2f}'.format(isweep, np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20])))
                    ax_high_f.plot(self.processedData['sweep'+str(isweep)]['rho_list_hop'][20:40]-self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40],
                                   self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'][20:40],
                                   label=r'spec ; isweep {} ; $\rho$ = {:.2f}'.format(isweep, np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40])))
                else:
                    ax_low_f.plot(self.processedData['sweep'+str(isweep)]['rho_list_hop'][:20], self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'][:20],
                                  label='spec ; isweep {}'.format(isweep))
                    ax_high_f.plot(self.processedData['sweep'+str(isweep)]['rho_list_hop'][20:40], self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'][20:40],
                                   label='spec ; isweep {}'.format(isweep))
            if plot_pearson:
                if delta_rho:
                    ax_low_f.plot(self.processedData['sweep'+str(isweep)]['rho_list_hop'][:20]-self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20],
                                  self.processedData['sweep'+str(isweep)]['maxcoh_list'][:20],
                                  label=r'pear ; isweep {} ; $\rho$ = {:.2f}'.format(isweep, np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20])))
                    ax_high_f.plot(self.processedData['sweep'+str(isweep)]['rho_list_hop'][20:40]-self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40],
                                   self.processedData['sweep'+str(isweep)]['maxcoh_list'][20:40],
                                   label=r'pear ; isweep {} ; $\rho$ = {:.2f}'.format(isweep, np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40])))
                else:
                    ax_low_f.plot(self.processedData['sweep'+str(isweep)]['rho_list_hop'][:20], self.processedData['sweep'+str(isweep)]['maxcoh_list'][:20], label='pearson ; isweep {}'.format(isweep))
                    ax_high_f.plot(self.processedData['sweep'+str(isweep)]['rho_list_hop'][20:40], self.processedData['sweep'+str(isweep)]['maxcoh_list'][20:40], label='pearson ; isweep {}'.format(isweep))
                
            # ax_low_f.axvline(np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20]))
            # ax_high_f.axvline(np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40]))
            # ax_low_f.text(np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20]), -.05, 'ref')
            # ax_high_f.text(np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40]), -.05, 'ref')
        
        
        ax_low_f.set_title(r'#{}'.format(self.shot))
        ax_high_f.set_title(r'#{}'.format(self.shot))
          
        my_legend(ax_low_f, fontsize="12")
        my_legend(ax_high_f, fontsize="12")
        
        if delta_rho:
            ax_low_f.set_xlabel(r'$\Delta_\rho$')
            ax_high_f.set_xlabel(r'$\Delta_\rho$')
        else:
            ax_low_f.set_xlabel(r'$\rho$')
            ax_high_f.set_xlabel(r'$\rho$')
        
        ax_high_f.set_ylabel(r'correlation')
        ax_low_f.set_ylabel(r'correlation')
        if ylog:
            ax_low_f.set_yscale('log')
            ax_high_f.set_yscale('log')
    
        return ax_low_f, ax_high_f
    
    



def plot_correlation_slopes(xdata, ydata, rho_s=None, rho_loc = None, ind_turb=None, ind_aval=None, ind_turb_pos=None, exclude_aval=None, caption=True, xunit='rho_s', ax=None, ylog=True, retdata=False, **kwargs):
    '''
    Assumes that xdata is given in delta/rho_s
    
    xunit: 'rho_s' or 'cm'
    '''
    
    default_kwargs = {'label':'spectral','marker':'o', 'color':'teal', 'linestyle':''}
    default_kwargs.update(kwargs)
    
    if xunit=='cm':
        xdata = xdata*rho_s
        
    if ax is None:
        fig, ax = plot_1d([], [], grid=True)
    
    ax.plot(xdata, (ydata), **default_kwargs)
    
    if caption:
        if rho_s is not None:
            my_text(ax, 0.15, 0.87, r'$\rho_s$ = {:.2f} mm'.format(rho_s*10), fontsize=12, color='k')
        if rho_loc is not None:
            my_text(ax, 0.35, 0.87, r'$\rho$ = {:.2f}'.format(rho_loc), fontsize=12, color='k')
       

    if ind_turb is not None:
        ind_turb_min=ind_turb[0]    
        ind_turb_max=ind_turb[1]
        popt, perr, rsquared = my_linearRegression(xdata[ind_turb_min:ind_turb_max], np.log(ydata[ind_turb_min:ind_turb_max]), mode='affine')
        ax.plot(xdata[ind_turb_min:ind_turb_max], np.exp(popt[0]*xdata[ind_turb_min:ind_turb_max]+ popt[1]), 'r', marker='')
        lc = 1/popt[0]
        lcerr = perr[0]/popt[0]**2
        print('Lc = {:.2f} +/- {:.2f} ; R² = {:.2f}'.format(lc, lcerr, rsquared))
        
        if caption:
            if xunit=='rho_s':
                my_text(ax, 0.8, 0.4, r'$l_c \approx$ {:.2f} +/- {:.2f} $\rho_s$'.format(lc, lcerr), fontsize=12, color='r')
            elif xunit=='cm':
                my_text(ax, 0.8, 0.4, r'$l_c \approx$ {:.2f} +/- {:.2f} $\cm$'.format(lc, lcerr), fontsize=12, color='r')
        
    else:
        lc = None
        lcerr = None
            
        
    if ind_aval is not None:
        ind_aval_min=ind_aval[0]
        ind_aval_max=ind_aval[1]
        if exclude_aval is not None:
            ### remove the point exclude
            xdata_excluded = np.delete(xdata[ind_aval_min:ind_aval_max], exclude_aval)
            ydata_excluded = np.delete(ydata[ind_aval_min:ind_aval_max], exclude_aval)
            popt, perr, rsquared = my_linearRegression(xdata_excluded, np.log(ydata_excluded), mode='affine')
            ax.scatter(xdata[ind_aval_min:ind_aval_max][exclude_aval], ydata[ind_aval_min:ind_aval_max][exclude_aval], color='red', marker='x', s=120)
        else:     
            popt, perr, rsquared = my_linearRegression(xdata[ind_aval_min:ind_aval_max], np.log(ydata[ind_aval_min:ind_aval_max]), mode='affine')
        ax.plot(xdata[ind_aval_min:ind_aval_max], np.exp(popt[0]*xdata[ind_aval_min:ind_aval_max] + popt[1]), 'g', marker='')
        laval = 1/popt[0]
        lavalerr = perr[0]/popt[0]**2
        Caval = ydata[ind_aval_max]
        print('La = {:.2f} +/- {:.2f} ; R² = {:.2f} ; Caval = {:.2f}'.format(laval, lavalerr, rsquared, Caval))
        Cavalerr = 0.1
        if caption:
            if xunit=='rho_s':
                my_text(ax, 0.8,0.25, r'$L_a \approx$ {:.2f} +/- {:.2f} $\rho_s$'.format(laval, lavalerr), fontsize=12, color='g')
            elif xunit=='cm':
                my_text(ax, 0.8,0.25, r'$L_a \approx$ {:.2f} +/- {:.2f} cm'.format(laval, lavalerr), fontsize=12, color='g')
            
    else:
        laval = None
        lavalerr = None
        Caval = None
        Cavalerr = None

    if ind_turb_pos is not None:
        ind_turb_min=ind_turb_pos[0]
        ind_turb_max=ind_turb_pos[1]
        popt, perr, rsquared = my_linearRegression(xdata[ind_turb_min:ind_turb_max], np.log(ydata[ind_turb_min:ind_turb_max]), mode='affine')
        ax.plot(xdata[ind_turb_min:ind_turb_max], np.exp(popt[0]*xdata[ind_turb_min:ind_turb_max]+ popt[1]), 'b', marker='')
        lcplus = abs(1/popt[0])
        lcpluserr = perr[0]/popt[0]**2
        print('Lc+ = {:.2f} +/- {:.2f} ; R² = {:.2f}'.format(lcplus, lcpluserr, rsquared))
        
        if caption:
            if xunit=='rho_s':
                my_text(ax, 0.8,0.1, r'$L_c \approx$ {:.2f} +/- {:.2f} $\rho_s$'.format(lcplus, lcpluserr), fontsize=12, color='b')
            elif xunit=='cm':
                my_text(ax, 0.8,0.1, r'$L_c \approx$ {:.2f} +/- {:.2f} cm'.format(lcplus, lcpluserr), fontsize=12, color='b')
    else:
        lcplus = None
        lcpluserr = None            
            

    if ylog:
        plt.yscale('log')
    
    plt.ylim(0.07,1.1)
    plt.ylabel('correlation')
    plt.axhline(1, color='black', linestyle='--')
    plt.axvline(0, color='black', linestyle='--')
    
    if xunit=='rho_s':
        plt.xlabel(r'$\Delta$ $[\rho_s]$')
    elif xunit=='cm':
        plt.xlabel(r'$\Delta$ $[cm]$')
    
    if retdata:
        data=dict()
        data['lc'] = lc
        data['lcerr'] = lcerr
        data['laval'] = laval
        data['lavalerr'] = lavalerr
        data['lcplus'] = lcplus
        data['lcpluserr'] = lcpluserr
        data['Caval'] = Caval
        data['Cavalerr'] = Cavalerr
        return ax, data
    








# if __name__ == '__main__':
    # a=CorrelationAnalysis(82609, numDemod=True)
    # a.plot_freq_pattern_sweep()
    # tcorr_loc, coh_list, fcsd_spectral, spectral_coh_list, maxcoh_list, maxspectralcoh_list = a.get_correlation_isweep(isweep=4, plateau=3, mode='amp', ret=True)

    # rho_ref, rho_hop = a.get_raytracing_isweep(isweep=4, ret=True)
    # fig, ax, im = plot_2d(np.transpose(coh_list), rho_hop, tcorr_loc, vmin=0, vmax=1)
    # ax.set_ylim(-2e-5, 2e-5)
    # ax.set_xlabel(r'$\rho$')
    # ax.set_ylabel(r'time delay')
    # plt.tight_layout()

    # fig, ax, im = plot_2d(np.transpose(abs(spectral_coh_list)**2), rho_hop, fcsd_spectral, vmin=0, vmax=1)
    # ax.set_ylim(-2e6, 2e6)
    # ax.set_xlabel(r'$\rho$')
    # ax.set_ylabel(r'frequency')
    # plt.tight_layout()




