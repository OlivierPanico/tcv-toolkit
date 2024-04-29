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
from diag_tcv.dbsAnalysis.correlationAnalysis.channelsCorrelation import get_raw_signals

#DBS
from DBS.io.read import get_DBS_params #Channels params (dt, freq, ...)
from DBS import definitions as defs #Location of data
from DBS.io.interface import DataInterface #General class for DBS data
from DBS.processing.sigprocessing import remove_signal_edges, get_normalized_complex_signal
from DBS.beamtracing.interface import Beam3dInterface #raytracing
from DBS.beamtracing.DBSbeam import _DBSbeam
#dataAnalysis
from dataAnalysis.utils.plot_utils import plot_1d
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
    
    def __init__(self, shot, xmode=1, numDemod=False, machine='tcv', verbose=False):
        super().__init__(shot=shot, verbose=verbose) #init with master class TCVShot
        
        self.machine = machine
        self.tmp_dir = defs.DATA_TMP_DIR / f'{self.machine}'
        
        self.modeCorrelation = isModeCorrelation(self.shot)
        self.numDemod = numDemod
    
        self.xmode = xmode    
    
        self._load_sweep_params()

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


    def load_from_sweep(self, isweep=1):  
        '''
        Load DataInterface from specific sweep
        '''
        if self.numDemod:
            self.data_ref = DataInterface(shot=self.shot, isweep=isweep, channelval=3, machine='tcv', verbose=self.verbose)
            self.data_hop = DataInterface(shot=self.shot, isweep=isweep, channelval=4, machine='tcv', verbose=self.verbose)
        else: 
            self.data_ref = DataInterface(shot=self.shot, isweep=isweep, channelval=1, machine='tcv', verbose=self.verbose)
            self.data_hop = DataInterface(shot=self.shot, isweep=isweep, channelval=2, machine='tcv', verbose=self.verbose)

    def load_from_time(self, time):
        '''
        Load DataInterface from time => will choose the right sweep
        '''
        if self.numDemod:
            self.data_ref = DataInterface.from_time(shot=self.shot, time=time, channelval=3, machine='tcv', verbose=self.verbose)
            self.data_hop = DataInterface.from_time(shot=self.shot, time=time, channelval=4, machine='tcv', verbose=self.verbose)
        else:
            self.data_ref = DataInterface.from_time(shot=self.shot, time=time, channelval=1, machine='tcv', verbose=self.verbose)
            self.data_hop = DataInterface.from_time(shot=self.shot, time=time, channelval=2, machine='tcv', verbose=self.verbose)
        
    def plot_freq_pattern_sweep(self, isweep=1):
        '''
        Plot the pattern used for the given sweep
        '''
        
        t0F_structured_ref = self.params_ref.t0F_structured[isweep-1,:]
        freq_ref = self.params_ref.F
        
        t0F_structured_hop = self.params_hop.t0F_structured[isweep-1,:]
        freq_hop = self.params_hop.F
        
        fig, ax = plot_1d([],[], grid=True)
        ax.set_title(' #{} sweep {} ; freq pattern'.format(self.shot, isweep))
        
        if self.same_nbpts_per_freq==True:
            ax.plot(t0F_structured_ref, freq_ref, color='blue', marker='+', label='ref')
            ax.plot(t0F_structured_hop, freq_hop, color='red', marker='+', label='hop')
        else:
            ax.scatter(t0F_structured_ref, freq_ref, color='blue', marker='+', label='ref')
            ax.plot(t0F_structured_hop, freq_hop, color='red', marker='+', label='hop')
            arrow_length = np.diff(data_ref.params.t0F)
            
            
            for i in range(len(freq_ref)):
                ax.arrow(t0F_structured_ref[i], freq_ref[i], arrow_length[i], 0,zorder=5, linewidth=1.5, head_width=0, head_length=0, ec='blue')
        
        ax.legend()
        ax.set_xlabel('time [s]')
        ax.set_ylabel('freq [GHz]')





    def get_raw_data_list(self, isweep, ifreq_list=None, time_window=None):
        '''
        load full isweep data for correlation
        if time_window is not False => then load the time window
        '''
        if ifreq_list is None:
            freq_ref = self.dbs_ref.params.F
            freq_hop = self.dbs_hop.params.F
            
            ifreq_list=freq_hop
        
        self.isweep = isweep            
        self.rawSigDic = get_raw_signals(self.shot, isweep, ifreq_list, numDemod=self.numDemod)
        
    
    def get_normalized_data_list(self):
        '''
        From raw data gives the normalized complex data without the transitory region at beginning of freq
        '''
        
        t = self.rawSigDic['t']
        freq_list_hop = self.rawSigDic['freq_list_hop']
        I_ref_list = self.rawSigDic['I_list_ref']
        Q_ref_list = self.rawSigDic['Q_list_ref']
        I_hop_list = self.rawSigDic['I_list_hop']
        Q_hop_list = self.rawSigDic['Q_list_hop']
        
        t_reduced, _ = remove_signal_edges(t[0,:], I_ref_list[0,:], dtstart=400.e-6, dtend=100.e-6, axis=-1)
       
        z_list_ref = np.zeros((len(freq_list_hop), len(t_reduced)), dtype='complex')
        z_list_hop = np.zeros((len(freq_list_hop), len(t_reduced)), dtype='complex')
        
        t_reduced_list_ref = np.zeros((len(freq_list_hop), len(t_reduced)))
        t_reduced_list_hop = np.zeros((len(freq_list_hop), len(t_reduced)))
        
        for i in range(len(freq_list_hop)):
            t_reduced_ref_loc,I_ref_loc = remove_signal_edges(t[i,:], I_ref_list[i,:], dtstart=400.e-6, dtend=100.e-6, axis=-1)
            _,Q_ref_loc = remove_signal_edges(t[i,:], Q_ref_list[i,:], dtstart=400.e-6, dtend=100.e-6, axis=-1)
    
            z_ref_loc = get_normalized_complex_signal(I_ref_loc, Q_ref_loc, self.rawSigDic['params_ref'].phase_cor)
            
            t_reduced_hop_loc,I_hop_loc = remove_signal_edges(t[i,:], I_hop_list[i,:], dtstart=400.e-6, dtend=100.e-6, axis=-1)
            _,Q_hop_loc = remove_signal_edges(t[i,:], Q_hop_list[i,:], dtstart=400.e-6, dtend=100.e-6, axis=-1)
    
            z_hop_loc = get_normalized_complex_signal(I_hop_loc, Q_hop_loc, self.rawSigDic['params_hop'].phase_cor)
        
            z_list_ref[i,:] = z_ref_loc
            z_list_hop[i,:] = z_hop_loc
            
            t_reduced_list_ref[i,:] = t_reduced_ref_loc
            t_reduced_list_hop[i,:] = t_reduced_hop_loc 
        
        corrSigDic = dict()
        corrSigDic['t'] = t
        corrSigDic['t_reduced_list_ref'] = t_reduced_list_ref
        corrSigDic['t_reduced_list_hop'] = t_reduced_list_hop
        corrSigDic['freq_list_ref'] = self.rawSigDic['freq_list_ref']
        corrSigDic['freq_list_hop'] = freq_list_hop
        corrSigDic['z_list_ref'] = z_list_ref
        corrSigDic['z_list_hop'] = z_list_hop
        
        self.corrSigDic = corrSigDic
        
    def get_coherence_list(self, method='spectral'):
        '''
        return the correlation spectra and max for the whole list of freq
        method: 
        spectral => using normalized cpsd
        time => using pearson correlation coefficients
        '''
        freq_list_hop = self.corrSigDic['freq_list_hop']

        #params
        nperseg=1024
        noverlap=512
        
        coh_list = np.zeros((len(freq_list_hop), nperseg))
        coh_max_list = np.zeros((len(freq_list_hop)))
        
        for i in range(len(freq_list_hop)):
            z_ref = self.corrSigDic['z_list_ref'][i]
            z_hop = self.corrSigDic['z_list_hop'][i]
            
            dt = self.corrSigDic['t'][i,1] - self.corrSigDic['t'][i,0]
            
            if self.modeCorrelation==True:
                z_ref = np.conjugate(z_ref)
            
            if method=='time':
                f, coh = custom_time_coherence(z_hop, z_ref, nperseg=1024, noverlap=512)
            else:
                f,coh = custom_coherence((z_hop), (z_ref), nperseg=1024,noverlap=512, dt=dt, window='hanning',remove_mean=True)
                
                
            coh_list[i,:] = coh            
            coh_max_list[i] = np.max(coh)
            
        self.corrSigDic['coh_list'] = coh_list
        self.corrSigDic['coh_max_list'] = coh_max_list
    
    def get_raytracing_data(self):
        freq_list_hop = self.corrSigDic['freq_list_hop']
        t = self.corrSigDic['t']
        tinit = t[0]
        tfin = t[-1]
        # interface = Beam3dInterface.from_shot('tcv', shot=self.shot, channelval=2, modex=1, twindow=[tinit, tfin], freq_choice=freq_list_hop)
        
        if self.numDemod:
            output_ref, beam3d_interface_ref = _DBSbeam('tcv',shot=self.shot, isweep=self.isweep, xmode=self.xmode, channelval=3, ifreqs=freq_list_hop, verbose=True, plot=True, load_if_existing=True)
            output_hop, beam3d_interface_hop = _DBSbeam('tcv',shot=self.shot, isweep=self.isweep, xmode=self.xmode, channelval=4, ifreqs=freq_list_hop, verbose=True, plot=True, load_if_existing=True)
        else:
            output_ref, beam3d_interface_ref = _DBSbeam('tcv',shot=self.shot, isweep=self.isweep, xmode=self.xmode, channelval=1, ifreqs=freq_list_hop, verbose=True, plot=True, load_if_existing=True)
            output_hop, beam3d_interface_hop = _DBSbeam('tcv',shot=self.shot, isweep=self.isweep, xmode=self.xmode, channelval=2, ifreqs=freq_list_hop, verbose=True, plot=True, load_if_existing=True)

        # outp = interface.fetch_result()
        
        self.output_ref = output_ref
        self.output_hop = output_hop
        
        self.rho_list_ref = output_ref.rho
        self.rho_list_hop = output_hop.rho
        
    
    def wrapper_coherence_analysis(self, isweep, ifreq_list, method='spectral'):
        self.get_raw_data_list(isweep, ifreq_list)
        self.get_normalized_data_list()
        self.get_coherence_list(method=method)
        # self.get_raytracing_data()
        
 

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
    
    
# %%
