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

import pickle as pkl
import sys
import time

def reformat(sim, params, anisotropy_ratio=False):
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
    
    por = [params['por_bk'], params['por_fz'], params['por_bk'], params['por_sl'], params['por_fp'], params['por_ts'], params['por_ss'], params['por_ws'], params['por_bd']]
    permX = [params['Kh_bk'], params['Kh_fz'], params['Kh_bk'], params['Kh_sl'], params['Kh_fp'], params['Kh_ts'], params['Kh_ss'], params['Kh_ws'], params['Kh_bd']]
    
    if anisotropy_ratio:    
        permZ = [params['Kh_bk']*params['Kr_bk'], params['Kh_fz']*params['Kr_fz'], params['Kh_bk']*params['Kr_bk'], params['Kh_sl']*params['Kr_sl'], params['Kh_fp']*params['Kr_fp'], params['Kh_ts']*params['Kr_ts'], params['Kh_ss']*params['Kr_ss'], params['Kh_ws']*params['Kr_ws'], params['Kh_bd']*params['Kr_bd']]
    else:
        permZ = [params['Kv_bk'], params['Kv_fz'], params['Kv_bk'], params['Kv_sl'], params['Kv_fp'], params['Kv_ts'], params['Kv_ss'], params['Kv_ws'], params['Kv_bd']]
    
        
    props = {
        'marker': [0, 1, 2, 3, 4, 5, 6, 7, 8],
        'por': por,
        'permX': permX, 
        'permZ': permZ
    }
    sim2 = {
        'FINAL_TIME': '{:.1f}d0 yr'.format(sim['tend']), 
        'OBSERVATION_INTERVAL': '{:.1f}d0 day'.format(sim['obs_int']), 
        'ALPHA': '{:.2e}'.format(params['alpha']).replace('e','d'), #Brooks-Corey
        'M': '{:.2f}d0'.format(params['m']), #Van Genuchten
        'LIQUID_RESIDUAL_SATURATION': '{:.2f}d0'.format(params['satresid']), #Brooks-Corey
        'MAX_RECHARGE': '{:f}d0'.format(params['max_recharge']), #m/day
        'PRESSURE_RIVER': '{:.2f}d0'.format(sim['riv_pressure']) #Pa
    } 
    return params2, props, sim2
    
def check_for_results(folder, exists):
    
    #sys.exit('Error: check_for_results has not been finished')
    return exists

def run(sim, params, folder, meshtype, overwrite=False, num=None, aniso=True):
    complete = False
    failed = False

    params2, props, sim2 = reformat(sim, params, anisotropy_ratio=aniso)    
    folder2, exists = make_directory(folder)
    if exists:
        exists = check_for_results(folder2,exists)
        
    if not exists:
        print('Simulating ' + folder2)
        try:
            real = Realization(meshtype, params2, props, sim2, folder=folder2)
            real.realize()
            complete = True
        except:
            failed = True
            print('    Failed: ' + folder2)
    elif exists and overwrite:
        print('Overwriting realization '+str(num)+'.')
        try:
            real = Realization(meshtype, params2, props, sim2, folder=folder2)
            real.realize()
            complete = True
        except:
            failed = True
            print('    Failed: ' + folder2)
            
    elif exists and not overwrite:
        print('Skipping realization '+str(num)+' (already exists).')
    else:
        sys.exit('Overwrite status encountered an error.')
        
    return complete, failed
    
def setup(sim, params, folder, meshtype, overwrite=False, num=None, aniso=True):
    complete = False
    failed = False

    params2, props, sim2 = reformat(sim, params, anisotropy_ratio=aniso)    
    folder2, exists = make_directory(folder)
    if not exists:
        print('Preparing ' + folder2)
        try:
            real = Realization(meshtype, params2, props, sim2, folder=folder2)
#            real.realize()
            complete = True
        except:
            failed = True
            print('    Failed: ' + folder2)
    elif exists and overwrite:
        print('Overwriting setup for realization '+str(num)+'.')
        try:
            real = Realization(meshtype, params2, props, sim2, folder=folder2)
#            real.realize()
            complete = True
        except:
            failed = True
            print('    Failed: ' + folder2)
            
    elif exists and not overwrite:
        print('Skipping realization '+str(num)+' (already exists).')
    else:
        sys.exit('Overwrite status encountered an error.')
        
    return complete, failed

def make_directory(mcfolder, show_warning=False):
    exists = False
    if not mcfolder.endswith('/'):
        mcfolder = mcfolder + '/'
    
    try:
        Path(mcfolder).mkdir()
    except OSError as e:
        exists = True
        if show_warning:
            warnings.warn(mcfolder+' already exists')
        
    return mcfolder, exists

def sample(parameter_dict, N, seed=111):
    '''
    This function takes a parameter dictionary and creates a pandas dataframe, with the keys as column headers and entries as either random samples from the distribution or constant parameters
    '''
#    np.random.seed(seed)
    tbl = pd.DataFrame()
    for key, value in parameter_dict.items():
        np.random.seed(seed)
        seed = seed+1 # for each key use a different seed
        if isinstance(value, int) or isinstance(value, float):
            tbl[key] = np.ones(N) * value
        else:
            tbl[key] = value.rvs(N)
    return tbl
    
    
def import_simulation(fname):
    '''
    This function takes a .mc file and reads it into parameter and simulation dictionaries
    '''
    with open(fname, 'rb') as f:
        sim = pkl.load(f)
        params = pkl.load(f)
        
    return sim, params
    
def save_simulation_setup(sim, params, fname):
    with open(fname,'wb') as f:
        pkl.dump(sim, f)
        pkl.dump(params, f)

class MonteCarlo:
    def __init__(self, mcfolder, sim=None, params=None, from_file=None, anisotropy_ratio=True, overwrite=False):
        '''
        This constructor does...
        '''
        
        self.mcfolder, exists = make_directory(mcfolder, show_warning=True)
        if exists and not overwrite:
            sys.exit('Folder already exists. Exiting. Use \"overwrite=True\" to overwrite data in this folder. You can still skip over specific realizations using \"MonteCarlo.realize(\'all\', overwrite=False)\".')
             
        self.anisotropy_ratio = anisotropy_ratio
        
        if from_file is not None:
            self.sim, self.params = import_simulation(from_file)
        else:
            self.sim = sim
            self.params = params
            save_simulation_setup(self.sim, self.params, self.mcfolder+'mc_setup.pkl')
            
        return
        
        
    def SampleParameters(self, N):
        '''
        This method does...
        '''
        self.tbl = sample(self.params, N)
        self.tbl.to_csv(self.mcfolder+'parameter_table.csv')
        return self.tbl
        
    def SetupRealization(self, number, overwrite=False, meshtype='FracZoneFloodplain'):
        '''
        
        '''
        nreals = 0
        nfails = 0
        st = time.time()
        if number is 'all':
            for i in range(len(self.tbl)):
                parameters = self.tbl.loc[i,:].to_dict()
                complete, fail = setup(self.sim, parameters, self.mcfolder+str(i), meshtype, num=i, overwrite=overwrite, aniso=self.anisotropy_ratio)
                nreals += complete
                nfails += fail
        else:
            parameters = self.tbl.loc[number,:].to_dict()
            complete, fail = setup(self.sim, parameters, self.mcfolder+str(number), meshtype, num=number, overwrite=overwrite, aniso=self.anisotropy_ratio)
            nreals += complete
            nfails += fail
            
        return

    def Realize(self, number, parallel=False, overwrite=False, meshtype='FracZoneFloodplain'):
        '''
        This method does...
        '''
        nreals = 0
        nfails = 0
        st = time.time()
        if number is 'all':
            for i in range(len(self.tbl)):
                parameters = self.tbl.loc[i,:].to_dict()
                complete, fail = run(self.sim, parameters, self.mcfolder+str(i), meshtype, num=i, overwrite=overwrite, aniso=self.anisotropy_ratio)
                nreals += complete
                nfails += fail
        else:
            parameters = self.tbl.loc[number,:].to_dict()
            complete, fail = run(self.sim, parameters, self.mcfolder+str(number), meshtype, num=number, overwrite=overwrite, aniso=self.anisotropy_ratio)
            nreals += complete
            nfails += fail
            
        et = time.time()
        simtime_mins = np.round((et-st)/60, 2)
        print('Realizations simulated: ' + str(nreals))
        print('Failed realizations: ' + str(nfails))
        print('Simulation time: ' + str(simtime_mins) + ' minutes')
        
        return
            
            
            
