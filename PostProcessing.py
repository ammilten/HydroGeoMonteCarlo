import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

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
    return sum([os.path.isdir(mcfolder+i) for i in os.listdir(mcfolder)])
    
    
def load_mc_data(mcfolder, N='all'):
    '''
    Loads a list of data from the realizations specified by N, located in \'mcfolder\'
    Options for N:
      \'all\': gets the total number of realizations in the folder, and loads all of them
      list of ints: loads each realization in the list
      int: loads from realization 0 to the realization specified by the int
      
    Outputs a list of dataframes
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
    data[0] = load_simulated_data(mcfolder + '0')
    nt = len(data[0])
    for i in N2: 
        folder = mcfolder + str(i)
        try:
            data[i] = load_simulated_data(folder)
            if len(data[i]) != nt:
                print('Did not converge: '+folder)
        except:
            print('Could not load '+folder)

    data = [data[i] for i in range(len(data)) if data[i] is not None]
    data = [data[i] for i in range(len(data)) if len(data[i]) == nt]
    
    return data
    
    
def pressure2wtd(pressure, screen_depth, rho=1000, g=9.8, pref=101325):
    if not isinstance(pressure, np.ndarray):
        pressure = pressure.values

    piez = (pressure - pref) / rho / g
    wtd = screen_depth - piez
    return wtd
    
    
def plot_simulation(folder, wells=['PLM1','PLM6'], show_observed=False, linewidth=3, cPLM1=[.5,.5,.5], cPLM6=[.68,.85,.9], tshift=365, tstart=360, dt=30):
    istart = int(np.ceil(tstart / dt))
    
    PLM1, PLM6, tetsu = load_well_data()
    
    if not folder.endswith('/'):
        folder = folder + '/'

    sim = load_simulated_data(folder)
    if 'PLM6' in wells or wells == 'PLM6':
        wtd6 = pressure2wtd(sim["PLM6 Liquid Pressure [Pa]"], np.mean(PLM6['screen_depth']))
        plt.plot(sim['Time [day]'][istart:]-tshift,wtd6[12:], color=cPLM6, label='PLM6 (model)', linewidth=linewidth)
            
    if 'PLM1' in wells or wells == 'PLM1':
        wtd1 = pressure2wtd(sim["PLM1 Liquid Pressure [Pa]"], np.mean(PLM1['screen_depth']))
        plt.plot(sim['Time [day]'][istart:]-tshift,wtd1[12:], color=cPLM1, label='PLM1 (model)', linewidth=linewidth)
  
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
            wtd6 = pressure2wtd(sim["PLM6 Liquid Pressure [Pa]"], np.mean(PLM6['screen_depth']))
            plt.plot(sim['Time [day]'][istart:]-tshift,wtd6[12:], color=cPLM6, label='PLM6 (real '+str(i)+')', linewidth=linewidth)
            
        if 'PLM1' in wells or wells == 'PLM1':
            wtd1 = pressure2wtd(sim["PLM1 Liquid Pressure [Pa]"], np.mean(PLM1['screen_depth']))
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

def WLPCA(mcfolder, N, tshift=1460):
    '''
    Fetches data in mcfolder and returns a PCA object
    '''

    dat = load_mc_data(mcfolder, N=N)
    PLM1, PLM6, tetsu = load_well_data()
    
    PLM1_pressure = np.concatenate([dat[i]['PLM1 Liquid Pressure [Pa]'].values[np.newaxis,:] for i in range(len(dat))], axis=0)
    PLM6_pressure = np.concatenate([dat[i]['PLM6 Liquid Pressure [Pa]'].values[np.newaxis,:] for i in range(len(dat))], axis=0)

    PLM1_wtd = pressure2wtd(PLM1_pressure, np.mean(PLM1['screen_depth']))
    PLM6_wtd = pressure2wtd(PLM6_pressure, np.mean(PLM6['screen_depth']))    
    
    WTD = np.concatenate((PLM1_wtd, PLM6_wtd), axis=1)
    
    interp_days = dat[0]['Time [day]']-tshift
    conditions = np.array(interp_days.values >= np.min(tetsu['day'].values)) & np.array(interp_days.values <= np.max(tetsu['day'].values))
    interp_days = interp_days[conditions]
    
    it2 = np.concatenate((interp_days, interp_days))
    WTD = WTD[:,np.where(it2)[0]]

    PLM1_tetsu = PLM6['Elevation'] - tetsu['PLM6 Piezometer WTE']
    PLM6_tetsu = PLM6['Elevation'] - tetsu['PLM6 Piezometer WTE']
    PLM1_interp = np.interp(interp_days, tetsu['day'], PLM1_tetsu)
    PLM6_interp = np.interp(interp_days, tetsu['day'], PLM6_tetsu) 
    WTD_tetsu = np.concatenate((PLM1_interp, PLM6_interp),axis=0)
    WTD_tetsu = WTD_tetsu[:,np.newaxis].T
    
    
    pca = PCA()
    pca.fit(WTD)
    scores = pca.transform(WTD)
    scores_tetsu = pca.transform(WTD_tetsu)
    
    return pca, scores, scores_tetsu


    

