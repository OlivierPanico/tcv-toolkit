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
from dataAnalysis.utils.plot_utils import plot_1d, my_legend, my_text
from dataAnalysis.spectral.spectralAnalysis import custom_coherence, custom_time_coherence
 
#TCV
from diag_tcv.shotAnalysis.dischargeInfoMdsObject import TCVShot


### ========== ###
### PARAMETERS ###
### ========== ###
cmap='viridis'
cmap='jet'

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
        
        self.modeCorrelation = isModeCorrelation(self.shot)
        self.numDemod = numDemod
    
        self.xmode = xmode    
        self.plot = plot
    
        self._load_sweep_params()
        
        self.processedData = dict()
        for i in range(self.nbsweep):
            self.processedData['sweep'+str(i+1)] = dict()

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
        
        
    def plot_freq_pattern_sweep(self, isweep=1):
        '''
        Plot the pattern used for the given sweep
        '''
        
        t0F_structured_ref = self.params_ref.t0F_structured[isweep-1,:]
        freq_ref = self.params_ref.F
        
        t0F_structured_hop = self.params_hop.t0F_structured[isweep-1,:]
        freq_hop = self.params_hop.F
        
        fig, ax = plot_1d([],[], grid=True)
        ax.set_title(' #{} sweep {}'.format(self.shot, isweep))
        
        if self.same_nbpts_per_freq==True:
            ax.plot(t0F_structured_ref, freq_ref, color='blue', label='ref')
            ax.plot(t0F_structured_hop, freq_hop, color='red', label='hop')
        else:
            ax.scatter(t0F_structured_ref, freq_ref, color='blue',label='ref')
            ax.plot(t0F_structured_hop, freq_hop, color='red', label='hop')
            arrow_length = np.diff(data_ref.params.t0F)
            
            
            for i in range(len(freq_ref)):
                ax.arrow(t0F_structured_ref[i], freq_ref[i], arrow_length[i], 0,zorder=5, linewidth=1.5, head_width=0, head_length=0, ec='blue')
        
        my_legend(ax)
        ax.set_xlabel('time [s]')
        ax.set_ylabel('freq [GHz]')



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
        self.processedData['sweep'+str(isweep)]['ifreq_list'] = ifreq_list
        self.processedData['sweep'+str(isweep)] = get_raw_signals(self.shot, isweep, ifreq_list,self.processedData['sweep'+str(isweep)],  numDemod=self.numDemod)
    
    def get_normalized_data_isweep(self, isweep, ifreq_list=None, dtsart=400.e-6, dtend=100.e-6, ret=False):
        '''
        From raw data gives the normalized complex
        '''
        if self.verbose:
            print('Getting normalized data for sweep {}'.format(isweep))
            
        if 'data_ref' not in self.processedData['sweep'+str(isweep)]:
            self._get_raw_data_isweep(isweep, ifreq_list)
        
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
        
        self.processedData['sweep'+str(isweep)]['dt'] = t[0,1] - t[0,0]
        self.processedData['sweep'+str(isweep)]['t_reduced_list_ref'] = t_reduced_list_ref
        self.processedData['sweep'+str(isweep)]['t_reduced_list_hop'] = t_reduced_list_hop
        self.processedData['sweep'+str(isweep)]['z_list_ref'] = z_list_ref
        self.processedData['sweep'+str(isweep)]['z_list_hop'] = z_list_hop
        
        if ret:
            return z_list_ref, z_list_hop, t_reduced_list_ref, t_reduced_list_hop
        
    def get_correlation_isweep(self, isweep, ifreq_list=None, nperseg=1024, noverlap=512, window=None, remove_mean=True, plot=False):
        
        if self.verbose:
            print('Performing correlation for sweep {}'.format(isweep))
        
        if 'z_list_hop' not in self.processedData['sweep'+str(isweep)]:
            self.get_normalized_data_isweep(isweep=isweep, ifreq_list=ifreq_list, dtsart=400.e-6, dtend=100.e-6)
    

        coh_list = np.zeros((len(self.processedData['sweep'+str(isweep)]['ifreq_list']), nperseg))
        spectral_coh_list = np.zeros((len(self.processedData['sweep'+str(isweep)]['ifreq_list']), nperseg))
        maxcoh_list = np.zeros((len(self.processedData['sweep'+str(isweep)]['ifreq_list'])))
        maxspectralcoh_list = np.zeros((len(self.processedData['sweep'+str(isweep)]['ifreq_list'])))
        
        dt=self.processedData['sweep'+str(isweep)]['dt']
        for i, ifreq in enumerate(self.processedData['sweep'+str(isweep)]['ifreq_list']):
            zref_loc = self.processedData['sweep'+str(isweep)]['z_list_ref'][ifreq]
            zhop_loc = self.processedData['sweep'+str(isweep)]['z_list_hop'][ifreq]

            tcorr_loc, coh_loc, fcsd_spectral, spectral_coh_loc = full_coherence_analysis(zref_loc, zhop_loc, dt, nperseg=nperseg, noverlap=noverlap, window=window, remove_mean=remove_mean, plot=self.plot, verbose=self.verbose)
            maxcoh = np.max(abs(coh_loc))
            maxspectralcoh = np.max(spectral_coh_loc)

            coh_list[i,:] = coh_loc
            spectral_coh_list[i,:] = spectral_coh_loc
            maxcoh_list[i] = maxcoh
            maxspectralcoh_list[i] = maxspectralcoh
        
        self.processedData['sweep'+str(isweep)]['coh_list'] = coh_list
        self.processedData['sweep'+str(isweep)]['spectral_coh_list'] = spectral_coh_list
        self.processedData['sweep'+str(isweep)]['maxcoh_list'] = maxcoh_list
        self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'] = maxspectralcoh_list
        

    def get_raytracing_isweep(self, isweep, ifreq_list=None):
        
        ###real frequencies of hop and ref
        # freq_list_hop = self.processedData['sweep'+str(isweep)]['freq_list_hop']
        # freq_list_ref = self.processedData['sweep'+str(isweep)]['freq_list_ref']
        
        ### index of frequencies to be used
        ifreq_list = self.processedData['sweep'+str(isweep)]['ifreq_list']
        
        # t_reduced_list_ref = self.processedData['sweep'+str(isweep)]['t_reduced_list_ref']
        
        # freq_list_hop = self.corrSigDic['freq_list_hop']
        # t = self.corrSigDic['t']
        # tinit = t[0]
        # tfin = t[-1]
        # interface = Beam3dInterface.from_shot('tcv', shot=self.shot, channelval=2, modex=1, twindow=[tinit, tfin], freq_choice=freq_list_hop)
        
        if self.numDemod:
            output_ref, beam3d_interface_ref = _DBSbeam('tcv',shot=self.shot, isweep=isweep, xmode=self.xmode, channelval=3, ifreqs=ifreq_list, verbose=self.verbose, plot=self.plot, load_if_existing=True)
            output_hop, beam3d_interface_hop = _DBSbeam('tcv',shot=self.shot, isweep=isweep, xmode=self.xmode, channelval=4, ifreqs=ifreq_list, verbose=self.verbose, plot=self.plot, load_if_existing=True)
        else:
            output_ref, beam3d_interface_ref = _DBSbeam('tcv',shot=self.shot, isweep=isweep, xmode=self.xmode, channelval=1, ifreqs=ifreq_list, verbose=self.verbose, plot=self.plot, load_if_existing=True)
            output_hop, beam3d_interface_hop = _DBSbeam('tcv',shot=self.shot, isweep=isweep, xmode=self.xmode, channelval=2, ifreqs=ifreq_list, verbose=self.verbose, plot=self.plot, load_if_existing=True)

        # outp = interface.fetch_result()

        self.processedData['sweep'+str(isweep)]['output_ref'] = output_ref
        self.processedData['sweep'+str(isweep)]['output_hop'] = output_hop
        self.processedData['sweep'+str(isweep)]['beam3d_interface_ref'] = beam3d_interface_ref
        self.processedData['sweep'+str(isweep)]['beam3d_interface_hop'] = beam3d_interface_hop
        self.processedData['sweep'+str(isweep)]['rho_list_ref'] = output_ref.rho
        self.processedData['sweep'+str(isweep)]['rho_list_hop'] = output_hop.rho
        
    
            
    def plot_coherence_delta(self, isweep_list, ax_low_f=None, ax_high_f=None,  plot_spectral=True, plot_pearson=True, ylog=False):

        if ax_low_f is None:
            fig, ax_low_f = plot_1d([],[], grid=True)
        if ax_high_f is None: 
            fig, ax_high_f = plot_1d([],[], grid=True)
        
        
        sweep_tinit = self.params_ref.TDIFDOP + self.params_ref.t0seq
        period = self.params_ref.Period*1e-3
        print('shape period', np.shape(period))
              
              
        for i, isweep in enumerate(isweep_list):
            if 'maxcoh_list' not in self.processedData['sweep'+str(isweep)]:
                self.get_correlation_isweep(isweep)
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
            
            
            sweep_tinit_loc = sweep_tinit[isweep-1]
            period_loc = period
            print(sweep_tinit_loc, sweep_tinit_loc+period_loc)
            rho_s_loc_low_f = self.get_rho_s_r_z([sweep_tinit_loc, sweep_tinit_loc+period_loc], np.mean(r_ref[:20]), np.mean(z_ref[:20]), np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20]))
            rho_s_loc_high_f = self.get_rho_s_r_z([sweep_tinit_loc, sweep_tinit_loc+period_loc], np.mean(r_ref[20:40]), np.mean(z_ref[20:40]), np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40]))
            
            #put the rho_s in cm
            rho_s_loc_low_f = rho_s_loc_low_f*100
            rho_s_loc_high_f = rho_s_loc_high_f*100
            
            print(rho_s_loc_low_f, rho_s_loc_high_f)
            
            if plot_spectral:
                ax_low_f.plot(delta[:20]/rho_s_loc_low_f,
                                self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'][:20],
                                label=r'spec ; isweep {} ; $\rho$ = {:.2f}'.format(isweep, np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20])))
                ax_high_f.plot(delta[20:40]/rho_s_loc_high_f,
                                self.processedData['sweep'+str(isweep)]['maxspectralcoh_list'][20:40],
                                label=r'spec ; isweep {} ; $\rho$ = {:.2f}'.format(isweep, np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40])))
            
            if plot_pearson:
                
                ax_low_f.plot(delta[:20]/rho_s_loc_low_f,
                                self.processedData['sweep'+str(isweep)]['maxcoh_list'][:20],
                                label=r'pear ; isweep {} ; $\rho$ = {:.2f}'.format(isweep, np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][:20])))
                ax_high_f.plot(delta[20:40]/rho_s_loc_high_f,
                                self.processedData['sweep'+str(isweep)]['maxcoh_list'][20:40],
                                label=r'pear ; isweep {} ; $\rho$ = {:.2f}'.format(isweep, np.mean(self.processedData['sweep'+str(isweep)]['rho_list_ref'][20:40])))

        
        ax_low_f.set_title(r'#{}'.format(self.shot))
        ax_high_f.set_title(r'#{}'.format(self.shot))
          
        my_legend(ax_low_f, fontsize="12")
        my_legend(ax_high_f, fontsize="12")
        

        ax_low_f.set_xlabel(r'$\Delta/\rho_s$')
        ax_high_f.set_xlabel(r'$\Delta/\rho_s$')

        ax_high_f.set_ylabel(r'correlation')
        ax_low_f.set_ylabel(r'correlation')
        if ylog:
            ax_low_f.set_yscale('log')
            ax_high_f.set_yscale('log')

        
        ax_low_f.set_xlim(-8, 3)
        ax_high_f.set_xlim(-8, 3)

        ax_low_f.set_ylim(0.1, 1)
        ax_high_f.set_ylim(0.1, 1)
        
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
    
    
    
    def plot_heating_with_sweeps(self, list_isweep):
    
        ax = self.plot_heating()
        
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
           
           
           
# if __name__ == '__main__':
#     a=CorrelationAnalysis(79797)
#     isweep_list=[3,4,5,6,7,8]
#     ifreq_list=np.linspace(20,39,20)
#     fig, ax = plot_1d([],[])
#     for i,sweeploc in enumerate(isweep_list):
#         a.wrapper_coherence_analysis(sweeploc, ifreq_list)
#         ax.plot(a.corrSigDic['freq_list_hop']-a.corrSigDic['freq_list_ref'],
#                 a.corrSigDic['coh_max_list'], marker='+')
    
    # a=CorrelationAnalysis(80338)
    # isweep_list=[2,3,4,5]
    # ifreq_list=np.linspace(20,39,20)
    # fig, ax = plot_1d([],[], grid=True, xlabel='delta freq', ylabel='coherence')
    # for i,sweeploc in enumerate(isweep_list):
    #     a.wrapper_coherence_analysis(sweeploc, ifreq_list)
    
    #     ax.plot(a.corrSigDic['freq_list_hop']-a.corrSigDic['freq_list_ref'],
    #             a.corrSigDic['coh_max_list'], marker='+', label='sweep = {} T=[{:.2f},{:.2f}]'.format(sweeploc, a.corrSigDic['t'][0][0], a.corrSigDic['t'][-1][0]))
    
    # # plt.yscale('log')
    # plt.legend()
    # plt.title('#{}'.format(a.shot))
    
    # isweep_list=[2,3,4,5]
    # ifreq_list=np.linspace(1,19,19)
    # fig, ax = plot_1d([],[], grid=True, xlabel='delta freq', ylabel='coherence')
    # for i,sweeploc in enumerate(isweep_list):
    #     a.wrapper_coherence_analysis(sweeploc, ifreq_list)
    
    #     ax.plot(a.corrSigDic['freq_list_hop']-a.corrSigDic['freq_list_ref'],
    #             a.corrSigDic['coh_max_list'], marker='+', label='sweep = {} T=[{:.2f},{:.2f}]'.format(sweeploc, a.corrSigDic['t'][0][0], a.corrSigDic['t'][-1][0]))
    
    # # plt.yscale('log')
    # plt.legend()
    # plt.title('#{}'.format(a.shot))
    
# %%

# a=CorrelationAnalysis(79797)
# isweep_list=[8]
# ifreq_list=np.linspace(20,39,20)
# fig, ax = plot_1d([],[], grid=True)
# for i,sweeploc in enumerate(isweep_list):
#     a.wrapper_coherence_analysis(sweeploc, ifreq_list)
#     ax.plot(a.rho_list_hop-a.rho_list_ref,
#             a.corrSigDic['coh_max_list'], marker='+')


# a=CorrelationAnalysis(80257)
# isweep_list=[2]
# ifreq_list=np.linspace(20,39,20)
# for i,sweeploc in enumerate(isweep_list):
#     a.wrapper_coherence_analysis(sweeploc, ifreq_list, method='time')
   
# #%%
# fig, ax = plot_1d([],[], grid=True)
# ax.plot(a.rho_list_hop[20:40], a.corrSigDic['coh_max_list'], marker='+', color='blue')
# ax.axvline(np.mean(a.rho_list_ref[20:40]), color='black')
# ax.set_xlabel(r'$\Delta_\rho$')
# ax.set_ylabel('coherence')
# plt.yscale('log')
# plt.title('#{} ; sweep {}'.format(80257, 2))


# a=CorrelationAnalysis(80914)
# isweep_list=[1]
# ifreq_list=np.linspace(1,19,19)
# fig, ax = plot_1d([],[])
# for i,sweeploc in enumerate(isweep_list):
#     a.wrapper_coherence_analysis(sweeploc, ifreq_list)
#     ax.plot(a.corrSigDic['freq_list_hop']-a.corrSigDic['freq_list_ref'],
#             a.corrSigDic['coh_max_list'], marker='+')


# a=CorrelationAnalysis(80949, numDemod=True)
# a.plot_freq_pattern_sweep()
# isweep_list=[4]
# ifreq_list=np.linspace(1,19,19)
# fig, ax = plot_1d([],[], grid=True)
# for i,sweeploc in enumerate(isweep_list):
#     a.wrapper_coherence_analysis(sweeploc, ifreq_list)
#     ax.plot(a.corrSigDic['freq_list_hop']-a.corrSigDic['freq_list_ref'],
#             a.corrSigDic['coh_max_list'], marker='+', color='blue')

# plt.xlabel('delta freq [GHz]')
# plt.ylabel('coherence')
    
    
# a=CorrelationAnalysis(81002, numDemod=True)
# a.plot_freq_pattern_sweep()
# isweep_list=[5]
# ifreq_list=np.linspace(20,39,20)
# fig, ax = plot_1d([],[], grid=True)
# for i,sweeploc in enumerate(isweep_list):
#     a.wrapper_coherence_analysis(sweeploc, ifreq_list)
#     ax.plot(a.corrSigDic['freq_list_hop']-a.corrSigDic['freq_list_ref'],
#             a.corrSigDic['coh_max_list'], marker='+')

# plt.xlabel('delta freq [GHz]')
# plt.ylabel('coherence')
    
    
# %% DEPRECATED

    # def get_normalized_data_list(self):
    #     '''
    #     From raw data gives the normalized complex data without the transitory region at beginning of freq
    #     '''
        
    #     t = self.rawSigDic['t']
    #     freq_list_hop = self.rawSigDic['freq_list_hop']
    #     I_ref_list = self.rawSigDic['I_list_ref']
    #     Q_ref_list = self.rawSigDic['Q_list_ref']
    #     I_hop_list = self.rawSigDic['I_list_hop']
    #     Q_hop_list = self.rawSigDic['Q_list_hop']
        
    #     t_reduced, _ = remove_signal_edges(t[0,:], I_ref_list[0,:], dtstart=400.e-6, dtend=100.e-6, axis=-1)
       
    #     z_list_ref = np.zeros((len(freq_list_hop), len(t_reduced)), dtype='complex')
    #     z_list_hop = np.zeros((len(freq_list_hop), len(t_reduced)), dtype='complex')
        
    #     t_reduced_list_ref = np.zeros((len(freq_list_hop), len(t_reduced)))
    #     t_reduced_list_hop = np.zeros((len(freq_list_hop), len(t_reduced)))
        
    #     for i in range(len(freq_list_hop)):
    #         t_reduced_ref_loc,I_ref_loc = remove_signal_edges(t[i,:], I_ref_list[i,:], dtstart=400.e-6, dtend=100.e-6, axis=-1)
    #         _,Q_ref_loc = remove_signal_edges(t[i,:], Q_ref_list[i,:], dtstart=400.e-6, dtend=100.e-6, axis=-1)
    
    #         z_ref_loc = get_normalized_complex_signal(I_ref_loc, Q_ref_loc, self.rawSigDic['params_ref'].phase_cor)
            
    #         t_reduced_hop_loc,I_hop_loc = remove_signal_edges(t[i,:], I_hop_list[i,:], dtstart=400.e-6, dtend=100.e-6, axis=-1)
    #         _,Q_hop_loc = remove_signal_edges(t[i,:], Q_hop_list[i,:], dtstart=400.e-6, dtend=100.e-6, axis=-1)
    
    #         z_hop_loc = get_normalized_complex_signal(I_hop_loc, Q_hop_loc, self.rawSigDic['params_hop'].phase_cor)
        
    #         z_list_ref[i,:] = z_ref_loc
    #         z_list_hop[i,:] = z_hop_loc
            
    #         t_reduced_list_ref[i,:] = t_reduced_ref_loc
    #         t_reduced_list_hop[i,:] = t_reduced_hop_loc 
        
    #     corrSigDic = dict()
    #     corrSigDic['t'] = t
    #     corrSigDic['t_reduced_list_ref'] = t_reduced_list_ref
    #     corrSigDic['t_reduced_list_hop'] = t_reduced_list_hop
    #     corrSigDic['freq_list_ref'] = self.rawSigDic['freq_list_ref']
    #     corrSigDic['freq_list_hop'] = freq_list_hop
    #     corrSigDic['z_list_ref'] = z_list_ref
    #     corrSigDic['z_list_hop'] = z_list_hop
        
    #     self.corrSigDic = corrSigDic
        
    # def get_coherence_list(self, method='spectral'):
    #     '''
    #     return the correlation spectra and max for the whole list of freq
    #     method: 
    #     spectral => using normalized cpsd
    #     time => using pearson correlation coefficients
    #     '''
    #     freq_list_hop = self.corrSigDic['freq_list_hop']

    #     #params
    #     nperseg=1024
    #     noverlap=512
        
    #     coh_list = np.zeros((len(freq_list_hop), nperseg))
    #     coh_max_list = np.zeros((len(freq_list_hop)))
        
    #     for i in range(len(freq_list_hop)):
    #         z_ref = self.corrSigDic['z_list_ref'][i]
    #         z_hop = self.corrSigDic['z_list_hop'][i]
            
    #         dt = self.corrSigDic['t'][i,1] - self.corrSigDic['t'][i,0]
            
    #         if self.modeCorrelation==True:
    #             z_ref = np.conjugate(z_ref)
            
    #         if method=='time':
    #             f, coh = custom_time_coherence(z_hop, z_ref, nperseg=1024, noverlap=512)
    #         else:
    #             f,coh = custom_coherence((z_hop), (z_ref), nperseg=1024,noverlap=512, dt=dt, window='hanning',remove_mean=True)
                
                
    #         coh_list[i,:] = coh            
    #         coh_max_list[i] = np.max(coh)
            
    #     self.corrSigDic['coh_list'] = coh_list
    #     self.corrSigDic['coh_max_list'] = coh_max_list
