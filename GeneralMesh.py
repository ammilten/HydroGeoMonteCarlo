from MeshTypes import FracZoneMesh, ShallowLayerFZMesh, FracZoneWithFloodplain, Uniform, FourLayers
import os
import sys

def exportPFLOTRANmesh(mesh3D, fname='mesh3d.ugi'):
    temp_name = 'temp_mesh3d.vtk'
#     fname = 'mesh3d.ugi'
    mesh3D.exportVTK(temp_name)
    with open(fname,"w") as f:
        f.write(str(mesh3D.cellCount())+" "+str(mesh3D.nodeCount())+"\n")
        with open(temp_name,"r") as tmp:
            txt = tmp.readlines()
            node_start = 5
            node_end = node_start + mesh3D.nodeCount()
            cell_start = node_end + 1
            cell_end = cell_start + mesh3D.cellCount()
            
            for i in range(cell_start,cell_end):
                nums = txt[i][1:].strip().split('\t')
                nums_swapped = nums[2] + " " + nums[1] + " " + nums[0] + " " + nums[5] + " " + nums[4] + " " + nums[3] + "\n"
                f.write("W " + nums_swapped)
            for i in range(node_start, node_end):
                pts = txt[i].strip().split('\t')
                swapped = pts[0]+" "+pts[2]+" "+pts[1]+"\n"
                f.write(swapped)
    os.remove(temp_name)
    return



class Mesh:
    def __init__(self, meshtype, params, folder=None):
        '''
        This constructor identifies which type of mesh (scenario) to create, uses the parameters provided to create the mesh, and then saves the mesh as a .ugi in the specified folder
    
        Inputs
          meshtype: string that identifies which type of mesh to construct. Options are below.
          params: dict of parameters needed for the corresponding meshtype
          folder: string specifying which folder to place the mesh.ugi in. If None, then the mesh.ugi file is not saved
      
        Outputs
          An instance of the class containing the following attributes:
            meshtype (input)
            params (input)
            folder (input)
            mesh: copied from the mesh class of pygimli
            topo: Nx2 numpy array of x,z coordinates of topography used in the mesh
    
        meshtype, params combos:
            'FracZone'
            'ShallowLayerFZ'
            'FracZoneFloodplain
            'FourLayers'
            'Uniform'
        '''
        self.meshtype = meshtype
        self.params = params
        
        print('   mesh is of type '+meshtype)
        
        if meshtype == 'FracZone':
            M = FracZoneMesh.FZMesh(
                efile = params['efile'],
                dip = params['dip'],
                H = params['H'],
                xpos = params['xpos'],
                xtra = params['xtra'],
                area = params['area'],
                deplist = params['deplist'],
                Q = params['Q'],
                dep = params['dep'])            
        elif meshtype == 'ShallowLayerFZ':
            M = ShallowLayerFZMesh.SLFZMesh(
                efile = params['efile'],
                dip = params['dip'],
                H = params['H'],
                xpos = params['xpos'],
                xtra = params['xtra'],
                area = params['area'],
                Q = params['Q'],
                dep = params['dep'],
                shdep = params['shdep'])
        elif meshtype == 'FracZoneFloodplain':
            M = FracZoneWithFloodplain.FPFZMesh(
                efile = params['efile'],
                dip = params['dip'],
                H = params['H'],
                xpos = params['xpos'],
                xtra = params['xtra'],
                area = params['area'],
                Q = params['Q'],
                dep = params['dep'],
                fpdep = params['fpdep'],
                fplim = params['fplim'])
        elif meshtype == 'Uniform':
            print('   mesh is confirmed Uniform')
            M = Uniform.UniformMesh(
                efile = params['efile'],
                xtra = params['xtra'],
                area = params['area'],
                Q = params['Q'],
                deplist = params['deplist'],
                dep = params['dep'])
        elif meshtype == 'FourLayers':
            M = FourLayers.FourLayerMesh(
                efile = params['efile'],
                xtra = params['xtra'],
                area = params['area'],
                Q = params['Q'],
                layerthx = params['layerthx'],
                dep = params['dep'])
        else:
            print('ERROR: Incorrect meshtype ('+meshtype+')')
            sys.exit('ERROR: Incorrect meshtype. Options are -FracZone- -ShallowLayerFZ- -Uniform- -FourLayers- -FracZoneFloodplain-')
           
        self.mesh = M.mesh3D
        self.mesh2d = M.mesh
        self.topo = M.topo
        self.geom = M.geom
        #M.show_mesh()
        
        if folder is not None:
            self.file = folder + 'mesh3d.ugi'
            exportPFLOTRANmesh(self.mesh, fname=self.file)
        


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
    
    M = Mesh(meshtype, params, folder=folder)
    
    
    #print(Mesh('example/dir').folder)
