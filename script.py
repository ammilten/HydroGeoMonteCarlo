import MonteCarlo as MC
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt

meshtype = 'Uniform'
mcfolder = 'data/test'
pflotran_path = '/home/ammilten/pflotran/src/pflotran/pflotran'

# These parameters are constant for all simulations and cannot be uncertain
simulation = {
    'efile':None,
    'Q':30,
    'deplist':[2.5], #FracZone only
    'area':500,
    'dep':100,
    'xtra':5,
    'fplim':[450,650],
    'tend':6, #yr
    'obs_int':30, #day
    'riv_pressure': 5000 #pa
}

# These parameters may change for each simulation
parameters = {
    'dip':st.uniform(60,80),
    'H':st.truncnorm(0,200,loc=75,scale=25),
    'xpos':st.truncnorm(0,500,loc=380,scale=30),
    'shdep':1, #ShallowLayerFZ only
    'fpdep':st.uniform(0.25, 2.75), #FracZoneFloodplain only
    'Kh_bk':st.loguniform(1e-15, 1e-12),
    'Kr_bk':st.uniform(0.1, 1.25),
    'por_bk':st.uniform(0.01, 0.07),
    'Kh_fz':st.loguniform(1e-14, 1e-12),
    'Kr_fz':st.uniform(0.1, 1.25),
    'por_fz':st.uniform(0.01, 0.15),
    'Kh_sl':5e-15, # Shallow Layer
    'Kr_sl':.5, # Shallow Layer
    'por_sl':0.03, # Shallow Layer
    'thx_ts':st.uniform(0.3, 0.7),
    'Kh_ts':st.loguniform(1e-15, 1e-12), #8.8e-13,
    'Kr_ts':st.uniform(0.1, 1.25),
    'por_ts':st.uniform(0.01, 0.25),
    'thx_ss':st.uniform(0.3, 0.7),
    'Kh_ss':st.loguniform(1e-15, 1e-12), #7.17e-13,
    'Kr_ss':st.uniform(0.1, 1.25),
    'por_ss':st.uniform(0.01, 0.25),
    'thx_ws':st.uniform(1., 4.),
    'Kh_ws':st.loguniform(1e-15, 1e-12), #1e-12,
    'Kr_ws':st.uniform(0.1, 1.25),
    'por_ws':st.uniform(0.01, 0.25),
    'Kh_bd':st.loguniform(1e-15, 1e-12), #1e-14,
    'Kr_bd':st.uniform(0.1, 1.25),
    'por_bd':st.uniform(0.01, 0.25),
    'Kh_fp':st.loguniform(5e-14, 5e-11),
    'Kr_fp':st.uniform(0.1, 1.25),
    'por_fp':st.uniform(0.01, 0.25),
    'alpha':st.uniform(5e-5, 5e-4 - 5e-5),
    'm':st.uniform(0.1, 0.9),
    'satresid':st.uniform(0.05, 0.15),
    'max_recharge':st.uniform(0.0005, 0.0045)
}

ex = MC.MonteCarlo(mcfolder, overwrite=True, pflotran_path=pflotran_path) #initialize an empty MonteCarlo simulation
ex.sim = simulation #set MonteCarlo simulation parameters
ex.params = parameters #set MonteCarlo uncertain parameters
ex.SampleParameters(N=1) #sample some parameters
ex.Realize('all', overwrite=True, meshtype=meshtype, parallel=False, nproc=8) #realize all parameters, with option to overwrite

