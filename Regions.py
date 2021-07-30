import numpy as np
import pandas as pd

def is_on_topo(topo, xq, zq):
    z_interp = np.interp(xq, topo[:,0], topo[:,1])
    return np.abs(z_interp-zq) < 0.01


def is_in_channel(xq, xlim):
    out = (xq >= xlim[0]) & (xq <= xlim[1])
    return out
    
    
def extract_surface_faces(obj, fname='srf.ss'):
    markers = [False] * obj.mesh.nodeCount()
    i = 0
    for node in obj.mesh.nodes():
        markers[i] = is_on_topo(obj.topo, node.x(), node.y())  
        i = i + 1

    ids = np.where(markers)[0]
    srf_cells = []
    for cell in obj.mesh.cells():
        inds = [None]*cell.nodeCount()
        node_ids = np.zeros(cell.nodeCount())
        i = 0
        for node in cell.allNodes():
            inds[i] = (node.id() in ids)
            node_ids[i] = node.id()
            i = i + 1
        if sum(inds) == 4:
            srf_cells.append(np.array(node_ids[inds]))


    srf_cells = np.flipud(np.array(srf_cells))
    df = pd.DataFrame(srf_cells.astype(np.int))
    df.insert(0,"Shp",['Q' for i in range(srf_cells.shape[0])])
    df.to_csv(fname,sep=" ",header=[str(srf_cells.shape[0]),'','','',''], index=False) 
    return
    
    
def extract_left_of_channel(obj, rivxlim=[580.5,605.5], rivdep=0.5, fname='left.ss'):
    markers = [False] * obj.mesh.nodeCount()
    i = 0
    for node in obj.mesh.nodes():
        istopo = is_on_topo(obj.topo, node.x(), node.y()) 
        markers[i] = (node.x() <= rivxlim[0]) & istopo
        i = i + 1

    ids = np.where(markers)[0]
    srf_cells = []
    srf_xmin = []
    for cell in obj.mesh.cells():
        inds = [None]*cell.nodeCount()
        node_ids = np.zeros(cell.nodeCount())
        i = 0
        minx = 10000
        for node in cell.allNodes():
            inds[i] = (node.id() in ids)
            node_ids[i] = node.id()
            if node.x() < minx:
                minx = node.x()
            i = i + 1
#         print(sum(inds))
        if sum(inds) == 4:
#             print(cell.id())
            srf_xmin.append(minx)
            srf_cells.append(np.array(node_ids[inds]))

    sorted_by_x_inds = np.argsort(srf_xmin)
    srf_cells = np.array(srf_cells)
    srf_cells = srf_cells[sorted_by_x_inds]
    df = pd.DataFrame(srf_cells.astype(np.int))
    df.insert(0,"Shp",['Q' for i in range(srf_cells.shape[0])])
    df.to_csv(fname,sep=" ",header=[str(srf_cells.shape[0]),'','','',''], index=False)
    return





def extract_right_of_channel(obj, rivxlim=[580.5,605.5], rivdep=0.5, fname='right.ss'):
    markers = [False] * obj.mesh.nodeCount()
    i = 0
    for node in obj.mesh.nodes():
        istopo = is_on_topo(obj.topo, node.x(), node.y())  
        markers[i] = (node.x() >= rivxlim[1]) & istopo
        i = i + 1

    ids = np.where(markers)[0]

    srf_cells = []
    for cell in obj.mesh.cells():
        inds = [None]*cell.nodeCount()
        node_ids = np.zeros(cell.nodeCount())
        i = 0
        for node in cell.allNodes():
            inds[i] = (node.id() in ids)
            node_ids[i] = node.id()
            i = i + 1
        if sum(inds) == 4:
            srf_cells.append(np.array(node_ids[inds]))


    srf_cells = np.array(srf_cells)
    df = pd.DataFrame(srf_cells.astype(np.int))
    df.insert(0,"Shp",['Q' for i in range(srf_cells.shape[0])])
    df.to_csv(fname,sep=" ",header=[str(srf_cells.shape[0]),'','','',''], index=False)
    return

def extract_channel_faces(obj, rivxlim=[580.5,605.5], rivdep=0.5, fname='riv.ss'):
    markers = [False] * obj.mesh.nodeCount()
    i = 0
    for node in obj.mesh.nodes():
        istopo = is_on_topo(obj.topo, node.x(), node.y())  
        markers[i] = is_in_channel(node.x(), rivxlim) & istopo
        i = i + 1

    ids = np.where(markers)[0]

    srf_cells = []
    for cell in obj.mesh.cells():
        inds = [None]*cell.nodeCount()
        node_ids = np.zeros(cell.nodeCount())
        i = 0
        for node in cell.allNodes():
            inds[i] = (node.id() in ids)
            node_ids[i] = node.id()
            i = i + 1
        if sum(inds) == 4:
            srf_cells.append(np.array(node_ids[inds]))


    srf_cells = np.array(srf_cells)
    df = pd.DataFrame(srf_cells.astype(np.int))
    df.insert(0,"Shp",['Q' for i in range(srf_cells.shape[0])])
    df.to_csv(fname,sep=" ",header=[str(srf_cells.shape[0]),'','','',''], index=False)
    return
    
class Regions:
    def __init__(self, GM, folder):
        if folder is not None:
            self.srf = folder+'srf.ss'
            self.left = folder+'left.ss'
            self.right = folder+'right.ss'
            self.riv = folder+'riv.ss'
        
            extract_surface_faces(GM, fname=self.srf)
            extract_left_of_channel(GM, fname=self.left)
            extract_right_of_channel(GM, fname=self.right)
            extract_channel_faces(GM, fname=self.riv)
            return
        else:
            print('WARNING: no folder specified for regions')
            return 'no folder specified'    

     

