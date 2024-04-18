#%%

#=== general info ===#
#Author: Olivier Panico
#contact: olivier.panico@free.fr

#Goal: performs correlation analysis for tcv data

#=== imports ===#
#General imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.io 

#Local imports (better to avoid)

#dataAnalysis
from dataAnalysis.utils.plot_utils import plot_1d, plot_2d
 
#TCV
from diag_tcv.dischargeInfoMdsObject import TCVShot

#DataPath

# data_path = '/home/panico/NoTivoli/astra_data/ASTRA_TCV_80376_0.400_1.935.mat'
data_path = '/home/panico/NoTivoli/astra_data/ASTRA_TCV_80328_0.400_2.000.mat'

class AstraAnalysis(TCVShot):
    
    def __init__(self, shot, verbose=False):
        super().__init__(shot=shot, verbose=verbose) #init with master class TCVShot
    
    def load_data(self, path):
        self.data = scipy.io.loadmat(path)
    
    def get_fluxes(self):
        self.rho = self.data['out']['RHOPSI'][0,0]
        self.T = self.data['out']['T'][0,0][0]
        self.QIEFF = self.data['out']['QIEFF'][0,0]
        self.QEEFF = self.data['out']['QEEFF'][0,0]
        
    def plot_fluxes(self, itime):
        fig, ax = plot_1d([],[], grid=True)
        ax.plot(self.rho[:,itime], self.QEEFF[:,itime], label='QEEFF')
        ax.plot(self.rho[:,itime], self.QIEFF[:,itime], label='QIEFF')
        ax.set_title(f'Time: {self.T[itime]}')
        ax.legend()
    
    def plot_fluxes_2d(self):
        X,Y = np.meshgrid(self.rho[:,100], self.T)
        Z = np.transpose((self.QEEFF/(self.QEEFF+self.QIEFF)))
        fig, ax,_ = plot_2d(Z,X,Y, cmap=plt.cm.jet)
        ax.contour(X,Y,Z, [0.5], colors='white')
        
    
# %%
