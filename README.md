# HydroGeoMonteCarlo

Notes:
1. PFLOTRAN must be installed
2. Will only work on my computer so far because the PFLOTRAN path is hardcoded


Example usage

```python
import MonteCarlo as MC
import scipy.stats as st
    
mcfolder = '/home/ammilten/Desktop/FZFPtest'
   
# These parameters are constant for all simulations
simulation = {
    'efile':'/home/ammilten/Documents/SCGSR/pflotran_dev/topo2.csv',
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
    'Kh_sl':5e-15,
    'Kr_sl':.5,
    'por_sl':0.03,
    'Kh_fp':st.loguniform(5e-14, 5e-11),
    'Kr_fp':st.uniform(0.1, 1.25),
    'por_fp':st.uniform(0.01, 0.25),
    'alpha':st.uniform(5e-5, 5e-4 - 5e-5),
    'm':st.uniform(0.1, 0.9),
    'satresid':st.uniform(0.05, 0.15),
    'max_recharge':st.uniform(0.0005, 0.0045)
}
    
mc_object = MC.MonteCarlo(mcfolder) # initialize object using destination folder
mc_object.sim = simulation # simulation setup
mc_object.params = parameters # uncertain parameters
mc_object.SampleParameters(N=5) # generate parameter table
mc_object.Realize('all') # perform simulations

# Data are stored at the path mc_object.mcfolder
```
    
    
