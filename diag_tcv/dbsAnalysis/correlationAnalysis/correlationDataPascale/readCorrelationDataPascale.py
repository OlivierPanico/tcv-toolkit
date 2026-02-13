#%% This file is aimed at reading the .mat files produced by the pmusic analysis of Pascale. The data are stored under dopcor.maxCpmus for the correlation values. All the sweeps and frequencies are concatenated in the same array. The code needs to extract one sweep (chosen by user) from this array. 

import scipy.io
import numpy as np
import os

def openPascaleCorrMatFile(shot_number):
    """
    Reads correlation data from .mat files produced by the pmusic analysis of Pascale.

    Parameters:
    shot_number (int): The shot number to read data for.

    Returns:
    mat_data (dict): The loaded .mat file data.
    """
    
    # start by loading the .mat file
    base_path = '/home/panico/tcv-toolkit/diag_tcv/dbsAnalysis/correlationAnalysis/correlationDataPascale/'  # Update this path to your data
    
    #name of the file is shot_number_corV12.mat
    file_name = f"{shot_number}_corV12.mat"
    file_path = os.path.join(base_path, file_name)
    
    mat_data = scipy.io.loadmat(file_path)
    return mat_data

def extractPascaleCorrSweepData(mat_data, sweep_index):
    """
    Extracts data for a specific sweep from the loaded .mat file data.

    Parameters:
    mat_data (dict): The loaded .mat file data.
    sweep_index (int): The index of the sweep to extract.

    Returns:
    sweep_data (np.ndarray): The extracted data for the specified sweep.
    """
    
    # Access NbStep: the number of frequencies per sweep
    NbStep = mat_data['dopcor']['NbStep'][0,0][0,0]
    
    # The correlation data for the sweep is isweep*NbStep:(isweep+1)*NbStep
    maxCpmus_allsweep = mat_data['dopcor']['maxCpmus'][0,0]
    maxCpmus = maxCpmus_allsweep[0, (sweep_index-1)*NbStep:(sweep_index)*NbStep]
    
    # correlation on the full complex signal
    Cmax_allsweep = mat_data['dopcor']['Cmax'][0,0]
    Cmax = Cmax_allsweep[0, (sweep_index-1)*NbStep:(sweep_index)*NbStep]
    
    # correlation on the amplitude signal
    Campmax_allsweep = mat_data['dopcor']['Campmax'][0,0]
    Campmax = Campmax_allsweep[0, (sweep_index-1)*NbStep:(sweep_index)*NbStep]
    
    return maxCpmus, Cmax, Campmax
   

# %%
