import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys
import pathlib
import os
from os import fspath
DEFAULTPATH = fspath(pathlib.Path(__file__).parent.absolute())
 

def load_well_data():
    '''
    Loads the well data for PLM1 and PLM6. Includes water level data from Tokunaga et al. (2019).
    '''
    PLM1 = {
        'Easting': 331012.302,
        'Northing': 4309678.924,
        'Elevation': 2786.739,
        'x':None,
        'screen_depth': [6.3, 7.3]}
    PLM6 = {
        'Easting': 331087.75,
        'Northing': 4309791.197,
        'Elevation': 2759.575,
        'x':None,
        'screen_depth': [6.1, 9.15]}
    
    PLM_waterlevel_file = DEFAULTPATH + '/data/PLM_waterlevel.csv'
    return PLM1, PLM6, pd.read_csv(PLM_waterlevel_file)


def load_simulated_data(folder):
    '''
    Loads simulated (observation) data from a specified simulation folder
    '''
    if not folder.endswith('/'):
        folder = folder + '/'
        
    datafile = folder + 'simulation-obs-0.tec'
    obs = pd.read_csv(datafile, delimiter='  | ', header=None, skiprows=1)
    obs.columns = ["Time [day]", "PLM1 Liquid Head [m]", "PLM1 Liquid Pressure [Pa]", "PLM1 Liquid Saturation", "PLM6 Liquid Head [m]", "PLM6 Liquid Pressure [Pa]", "PLM6 Liquid Saturation"]
    return obs
    
    
def get_number_of_reals(mcfolder):
    return sum(os.path.isdir(i) for i in os.listdir(mcfolder))
    
    
def load_mc_data(mcfolder, N='all'):
    '''
    Loads a list of data from the realizations specified by N, located in \'mcfolder\'
    Options for N:
      \'all\': gets the total number of realizations in the folder, and loads all of them
      list of ints: loads each realization in the list
      int: loads from realization 0 to the realization specified by the int
    '''
    if not mcfolder.endswith('/'):
        mcfolder = mcfolder + '/'
        
    if N is 'all':
        N2 = range(get_number_of_reals(mcfolder))
    elif isinstance(N, list):
        N2 = N
    elif isinstance(N, int):
        N2 = range(N)
    else:
        sys.exit('Error: \'N\' is not \'all\', a list, nor an integer')
    
    data = [None] * len(N2)
    for i in range(N2): 
        folder = mcfolder + str(i) + '/simulation-obs-0.tec'
        data[i] = load_simulated_data(folder)

    return data
    
    
def plot_realization(mcfolder, real, wells=['PLM1','PLM6'], show_observed=False, linewidth=3, cPLM1=[.5,.5,.5], cPLM6=[.68,.85,.9], tshift=365, tstart=360, dt=30):
    '''
    Plots water level elevation data from one or multiple realizations.
    Inputs
      mcfolder: folder of the Monte Carlo data
      real: Realization(s) to plot 
         \'all\': plots all realizations in the mcfolder
         list of integers: plots all realizations specified by integers in the list
         range object: plots all realizations in the range
         int: plots the realization specified by the integer
       show_observed: if True, loads data from Tokunaga et al. (2019) and plots alongside the realization
       linewidth: line width for each line
       cPLM1: color of the PLM1 line
       cPLM6: color of the PLM6 line
       tshift: time (days) to shift the simulated water table data. Used to shift where Tetsu's data lines up with the simulated data.
       tstart: time (days) at which to start displaying simulated water table data. Used to discard the first few simulation years
       dt: time (days) between snapshots
    '''

    istart = int(np.ceil(tstart / dt))

    rho = 1000
    g = 9.8
    pref = 101325.
    
    PLM1, PLM6, tetsu = load_well_data()
    
    if not mcfolder.endswith('/'):
        mcfolder = mcfolder + '/'
    
    if real is 'all':
        N = range(get_number_of_reals(mcfolder))
    elif isinstance(real, list) or isinstance(real, range):
        N = real
    elif isinstance(real, int):
        N = [real]
    else:
        sys.exit('Error: Format of argument \'real\' not recognized')
    
    
    for i in N:
        sim = load_simulated_data(mcfolder+str(i))
        
        

        if 'PLM6' in wells or wells == 'PLM6':
            piez6 = (sim["PLM6 Liquid Pressure [Pa]"]-pref) / rho / g
            wte6 = piez6 + PLM6['Elevation'] - np.mean(PLM6['screen_depth'])
            wtd6 = PLM6['Elevation'] - wte6
            plt.plot(sim['Time [day]'][istart:]-tshift,wtd6[12:], color=cPLM6, label='PLM6 (real '+str(i)+')', linewidth=linewidth)
            
        if 'PLM1' in wells or wells == 'PLM1':
            piez1 = (sim["PLM1 Liquid Pressure [Pa]"]-pref) / rho / g
            wte1 = piez1 + PLM1['Elevation'] - np.mean(PLM1['screen_depth'])
            wtd1 = PLM1['Elevation'] - wte1
            plt.plot(sim['Time [day]'][istart:]-tshift,wtd1[12:], color=cPLM1, label='PLM1 (real '+str(i)+')', linewidth=linewidth)
  
    if show_observed:
        if 'PLM6' in wells or wells == 'PLM6':
            plt.plot(tetsu['day'], PLM6['Elevation'] - tetsu['PLM6 Piezometer WTE'],'-b', linewidth=linewidth,label='PLM6 (Obs)')
        if 'PLM1' in wells or wells == 'PLM1':
            plt.plot(tetsu['day'], PLM1['Elevation'] - tetsu['PLM1 Piezometer WTE'],'-k', linewidth=linewidth, label='PLM1 (Obs)')
        
    plt.gca().invert_yaxis()
    plt.ylabel('Water Depth (mbgs)',fontsize=16)
    plt.xlabel('Days since Jan. 1, 2017',fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    return



    

