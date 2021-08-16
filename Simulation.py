import pathlib
from os import fspath
DEFAULTPATH = fspath(pathlib.Path(__file__).parent.absolute())


class Simulation:
    def __init__(self, simparams, folder=None, template=None):
        '''
        This constructor copies the template "in" file to the destination folder, and adds in the parameters in "simparams" to the empty spaces
        '''
        if template is None:
            template = DEFAULTPATH + '/template.in'
        
        self.simparams = simparams
        
        if folder is not None:
            self.file = folder+'simulation.in'
        
            lines = open(template).read().splitlines()
            lines[70] = lines[70] + simparams['ALPHA']
            lines[71] = lines[71] + simparams['LAMBDA']
            lines[72] = lines[72] + simparams['LIQUID_RESIDUAL_SATURATION']
            lines[84] = lines[84] + simparams['FINAL_TIME']
            lines[94] = lines[94] + simparams['OBSERVATION_INTERVAL']
            lines[107] = lines[107] + simparams['OBSERVATION_INTERVAL']
            lines[189] = lines[189] + simparams['MAX_RECHARGE']
            lines[207] = lines[207] + simparams['PRESSURE_RIVER']
            open(self.file,'w').write('\n'.join(lines))
        return
        

if __name__ == "__main__":
    folder = '/home/ammilten/Desktop/blah/'
    sim = {
        'FINAL_TIME': '6.d0 yr', #85
        'OBSERVATION_INTERVAL': '30 day', #95, 107
        'ALPHA': '1.d-4', #71
        'LAMBDA': '0.3d0', #72
        'LIQUID_RESIDUAL_SATURATION': '0.1d0', #73
        'MAX_RECHARGE': '0.00175', #189
        'PRESSURE_RIVER': '5000'} #207
        
    Simulation(sim, folder=folder)
        
