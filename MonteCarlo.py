# Author: Alex Miltenberger ammilten@stanford.edu

# This file contains:
# 1. A function to generate samples from a dictionary, and save to a csv
#    - Dictionary structure: 
#       dict = {'column_title':scipy_distribution_obj}
#       dict = {'column_title':constant_float_or_int}
# 2. A wrapper to take a .mc (ascii) file and turn it into the dict above
# 3. A wrapper to take a row in the parameter table and generate a realization folder
# 4. A wrapper to run a realization folder
# 5. Utilities for importing, overwriting, command line access

import numpy as np
import pandas as pd
import scipy.stats as st
from pathlib import Path
import warnings
from Realization import Realization

def reformat(sim, params):
    params2 = {
        'efile':sim['efile'],        
        'dip':params['dip'],
        'H':params['H'],
        'xpos':params['xpos'],
        'xtra':sim['xtra'],
        'Q':sim['Q'],
        'deplist':sim['deplist'], #FracZone only
        'area':sim['area'],
        'dep':sim['dep'],
        'shdep':params['shdep'], #ShallowLayerFZ only
        'fpdep':params['fpdep'], #FracZoneFloodplain only
        'fplim':sim['fplim'] #FracZoneFloodplain only
    } 
    props = {
        'marker': [0, 1, 2, 3, 4],
        'por': [params['por_bk'], params['por_fz'], params['por_bk'], params['por_sl'], params['por_fp']],
        'permX': [params['Kh_bk'], params['Kh_fz'], params['Kh_bk'], params['Kh_sl'], params['Kh_fp']], 
        'permZ': [params['Kv_bk'], params['Kv_fz'], params['Kv_bk'], params['Kv_sl'], params['Kv_fp']]
    }
    sim2 = {
        'FINAL_TIME': '{:.1f}d0 yr'.format(sim['tend']), 
        'OBSERVATION_INTERVAL': '{:.1f}d0 day'.format(sim['obs_int']), 
        'ALPHA': '{:.2e}'.format(params['alpha']).replace('e','d'), #Brooks-Corey
        'LAMBDA': '{:.2f}d0'.format(params['lambda']), #Brooks-Corey
        'LIQUID_RESIDUAL_SATURATION': '{:.2f}d0'.format(params['satresid']), #Brooks-Corey
        'MAX_RECHARGE': '{:f}d0'.format(params['max_recharge']), #m/day
        'PRESSURE_RIVER': '{:.2f}d0'.format(sim['riv_pressure']) #Pa
    } 
    return params2, props, sim2

def run(sim, params, folder):
    meshtype = 'FracZoneFloodplain'
    params2, props, sim2 = reformat(sim, params)    
    folder2 = make_directory(folder)
    real = Realization(meshtype, params2, props, sim2, folder=folder2)
    real.realize()
    return

def make_directory(mcfolder):
    if not mcfolder.endswith('/'):
        mcfolder = mcfolder + '/'
    
    try:
        Path(mcfolder).mkdir()
    except OSError as e:
        warnings.warn('Warning: '+mcfolder+' already exists')
        
    return mcfolder

def sample(parameter_dict, N, seed=111):
    '''
    This function takes a parameter dictionary and creates a pandas dataframe, with the keys as column headers and entries as either random samples from the distribution or constant parameters
    '''
    
    tbl = pd.DataFrame()
    for key, value in parameter_dict.items():
        if isinstance(value, int) or isinstance(value, float):
            tbl[key] = np.ones(N) * value
        else:
            tbl[key] = value.rvs(N)
    return tbl
    
    
def import_simulation(fname):
    '''
    This function takes a .mc file and reads it into parameter and simulation dictionaries
    '''
    sim = None
    params = None
    return sim, params

class MonteCarlo:
    def __init__(self, mcfolder, sim=None, params=None, from_file=None):
        '''
        This constructor does...
        '''
        self.mcfolder = make_directory(mcfolder)
        
        
        if from_file is not None:
            self.sim, self.params = import_simulation(from_file)
        else:
            self.sim = sim
            self.params = params
        return
        
    def SampleParameters(self, N):
        '''
        This method does...
        '''
        self.tbl = sample(self.params, N)
        self.tbl.to_csv(self.mcfolder+'parameter_table.csv')
        return self.tbl

    def realize(self, number, parallel=False, overwrite=False):
        '''
        This method does...
        '''
        if number is 'all':
            print('Not yet implemented')
        else:
            parameters = self.tbl.loc[number,:].to_dict()
            run(self.sim, parameters, self.mcfolder+str(number))
            
        return
            
            
            
