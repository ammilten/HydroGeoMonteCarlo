from GeneralMesh import Mesh
from Properties import Properties
from Regions import Regions
from Simulation import Simulation

import os 

import os

def fix_xmf(filename):
    dirname = os.path.dirname(filename)
    
    with open(filename, 'r') as file:
        filedata = file.read()
        
    filedata = filedata.replace(dirname+'/simulation.h5', 'simulation.h5')
    
    with open(filename, 'w') as file:
        file.write(filedata)
        
    return
        
def fix_all_xmfs(folder):
    for file in os.listdir(folder):
        if file.endswith(".xmf"):
            fix_xmf(folder + '/' + file)
    
    return
 

class Realization:
    def __init__(self, meshtype, meshparams, props, simparams, folder=None):

        self.folder = folder
#        print('  meshing '+folder)
        self.Mesh = Mesh(meshtype, meshparams, folder=folder)
#        print('  assigning properties to '+folder)
        self.Properties = Properties(self.Mesh, props, folder=folder)
#        print('  defining regions for '+folder)
        self.Regions = Regions(self.Mesh, folder=folder)
#        print('  preparing simulation files for '+folder)
        self.Simulation = Simulation(simparams, folder=folder)
        
    @classmethod
    def default(self):
        meshtype = 'ShallowLayerFZ'
        folder = '/home/ammilten/Desktop/blah/'
        params = {
            'efile':'/home/ammilten/Documents/SCGSR/pflotran_dev/topo2.csv',
            'dip':120,
            'H':100,
            'xpos':400,
            'xtra':5,
            'Q':30,
            'deplist':[2.5],
            'area':500,
            'dep':100,
            'shdep':30}
        props = {
            'marker': [0, 1, 2, 3],
            'por': [0.03, 0.07, 0.03, 0.10],
            'permX': [3.2e-15, 4e-14, 3.2e-15, 3.2e-12],
            'permZ': [1.6e-15, 2e-14, 1.6e-15, 1.6e-12]}
        sim = {
            'FINAL_TIME': '6.d0 yr', #85
            'OBSERVATION_INTERVAL': '30 day', #95, 107
            'ALPHA': '1.d-4', #71
            'LAMBDA': '0.3d0', #72
            'LIQUID_RESIDUAL_SATURATION': '0.1d0', #73
            'MAX_RECHARGE': '0.00175', #189
            'PRESSURE_RIVER': '5000'} #207
            
        return Realization(meshtype, params, props, sim, folder=folder)
        
        
    def realize(self, pflotran_path='/home/ammilten/pflotran/src/pflotran/pflotran', nproc=1):
        cmd = "mpirun -n " + str(nproc) + " " + pflotran_path + " -pflotranin " + self.Simulation.file
        print(cmd)
        os.system(cmd)
        
        fix_all_xmfs(self.folder)
        
        return

if __name__ == "__main__":
    meshtype = 'ShallowLayerFZ'
    folder = '/home/ammilten/Desktop/blah/'
    params = {
        'efile':'/home/ammilten/Documents/SCGSR/pflotran_dev/topo2.csv',
        'dip':120,
        'H':100,
        'xpos':400,
        'xtra':5,
        'Q':30,
        'deplist':[2.5],
        'area':500,
        'dep':100,
        'shdep':30}
    props = {
        'marker': [0, 1, 2, 3],
        'por': [0.03, 0.07, 0.03, 0.10],
        'permX': [3.2e-15, 4e-14, 3.2e-15, 3.2e-12],
        'permZ': [1.6e-15, 2e-14, 1.6e-15, 1.6e-12]}
    sim = {
        'FINAL_TIME': '6.d0 yr', #85
        'OBSERVATION_INTERVAL': '30 day', #95, 107
        'ALPHA': '1.d-4', #71
        'LAMBDA': '0.3d0', #72
        'LIQUID_RESIDUAL_SATURATION': '0.1d0', #73
        'MAX_RECHARGE': '0.00175', #189
        'PRESSURE_RIVER': '5000'} #207
    Realization.default() 
    #real = Realization(meshtype, params, props, sim)
    
