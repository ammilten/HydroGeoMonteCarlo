import h5py
import pandas as pd
import numpy as np

def exportRegionalDataset(mesh3D, props, fname):
    '''
    
    props is a dataframe with the following columns {marker, por, permX, permZ}
    '''

    hf = h5py.File(fname,'w')
    cellIds = np.arange(1, mesh3D.cellCount()+1)
    
    marker = np.arange(mesh3D.cellCount())
    for cell in mesh3D.cells():
        marker[cell.id()] = cell.marker()
        
    por = np.ones(mesh3D.cellCount())
    permX = np.ones(mesh3D.cellCount())
    permZ = np.ones(mesh3D.cellCount())
    
    #print(np.unique(marker))
    #print(props)
    for i in range(len(props['marker'])):
        por[marker==props['marker'][i]] = props['por'][i]
        permX[marker==props['marker'][i]] = props['permX'][i]
        permZ[marker==props['marker'][i]] = props['permZ'][i]
        
    
    hf.create_dataset('Cell Ids',data=cellIds)
    hf.create_dataset('Porosity', data=por)
    hf.create_dataset('Permeability', data=permX)
    hf.create_dataset('PermX', data=permX)
    hf.create_dataset('PermY', data=permX)
    hf.create_dataset('PermZ', data=permZ)
    hf.close()
    return
    
    
class Properties:
    def __init__(self, GM, props, folder=None):
        '''
        This class takes a 
        '''
        self.props = props
        
        if folder is not None:
            self.file = folder+'properties.h5'
            exportRegionalDataset(GM.mesh, props, self.file)
        
        
        
        
if __name__ == "__main__":
    props = {
        'marker': [0, 1, 2, 3],
        'por': [0.03, 0.07, 0.03, 0.10],
        'permX': [3.2e-15, 4e-14, 3.2e-15, 3.2e-12],
        'permZ': [1.6e-15, 2e-14, 1.6e-15, 1.6e-12]}
    props = pd.DataFrame(props)
    print(props)
        
