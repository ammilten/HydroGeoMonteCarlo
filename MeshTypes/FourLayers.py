
# -----------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate 

import pygimli as pg
import pygimli.meshtools as mt
import pygimli.physics.ert as ert

import time
import pathlib
from os import fspath
DEFAULTPATH = fspath(pathlib.Path(__file__).parent.absolute())

# ----------Create Topo Polygon ----------------------
def createBackground(efile, xextra, botdep, shdep):
#    eloc = pd.read_csv(efile,sep='\t',header=None)
#    nelec = eloc.values.shape[0]
#    topo = eloc.values[:,(1,3)]
    
    topo = np.genfromtxt(efile ,delimiter=',',skip_header=True)
    topo[:,1] = topo[:,1] - shdep
    
    e_bot = np.floor(np.min(topo[:,1]))
    xL = np.floor(np.min(topo[:,0]))
    xR = np.ceil(np.max(topo[:,0]))

    L_ext = np.array([[xL-xextra, topo[0,1]]])
    R_ext = np.array([[xR+xextra, topo[topo.shape[0]-1,1]]])
        
    bots = np.array([[R_ext[0,0], e_bot-botdep+shdep],[L_ext[0,0],e_bot-botdep+shdep]])
    combo = np.concatenate((topo,R_ext,bots,L_ext))

    bpoly = mt.createPolygon(combo, isClosed=True, marker=8)

    return topo, bpoly

# ----------Create Shallow Layer ----------------------
def createLayer2(efile, xextra, shdep):
#    eloc = pd.read_csv(efile,sep='\t',header=None)
#    nelec = eloc.values.shape[0]
#    topo = eloc.values[:,(1,3)]
    
    topo = np.genfromtxt(efile ,delimiter=',',skip_header=True)
    
    e_bot = np.floor(np.min(topo[:,1]))
    xL = np.floor(np.min(topo[:,0]))
    xR = np.ceil(np.max(topo[:,0]))

    L_ext = np.array([[xL-xextra, topo[0,1]]])
    R_ext = np.array([[xR+xextra, topo[topo.shape[0]-1,1]]])
    topo_ext = np.concatenate((L_ext, topo, R_ext))
    topo_bot = np.copy(np.flipud(topo_ext))
    topo_bot[:,1] = topo_bot[:,1]-shdep
    combo = np.concatenate((topo_ext, topo_bot))

    lpoly = mt.createPolygon(combo, isClosed=True, marker=3)

    return topo, lpoly

# ----------Create ANY Layer ----------------------
def createLayer(bot_topo, xextra, thx, marker):

    elev_bot = np.floor(np.min(bot_topo[:,1]))
    xL = np.floor(np.min(bot_topo[:,0]))
    xR = np.ceil(np.max(bot_topo[:,0]))

    L_ext = np.array([[xL-xextra, bot_topo[0,1]]])
    R_ext = np.array([[xR+xextra, bot_topo[bot_topo.shape[0]-1,1]]])
    topo_ext = np.concatenate((L_ext, bot_topo, R_ext))
    topo = np.copy(np.flipud(topo_ext))
    topo[:,1] = topo[:,1]+thx
    combo = np.concatenate((topo_ext, topo))

    markerpos = [topo[0,0], topo[0,1]-thx/2]
    lpoly = mt.createPolygon(combo, isClosed=True, marker=marker, markerPosition=markerpos)
    
    topo2 = np.copy(bot_topo)
    topo2[:,1] = topo2[:,1]+thx

    return topo2, lpoly

# ---------- Create mesh ---------------------
def createMesh(geom, topo, Q, area=None):

    #scheme = pg.DataContainerERT()
    #scheme.setSensorPositions(topo)

    if area is None:
        mesh = mt.createMesh(geom, quality=Q)
    else:
        mesh = mt.createMesh(geom, quality=Q, area=area)

    return mesh
    
class FourLayerMesh:
    def __init__(self, layerthx=[0.5, 0.5, 2.6], efile=None, dep=200, xtra=500, Q=20, zthx=5, outfile=None, area=None):
        self.dep = dep
        self.layerthx = layerthx
        self.Q = Q
        self.outfile = outfile
        self.xtra = xtra
        self.area = area
        self.zthx = zthx
        if efile is None:
            self.efile = DEFAULTPATH+"/data/PH-2018-eloc.txt"
        else:
            self.efile = efile

            
#        print('Creating Background')
        topo, bedrock = createBackground(self.efile, self.xtra, self.dep, sum(self.layerthx))          
        topo, weathered = createLayer(topo, self.xtra, self.layerthx[2], 7)
        topo, subsoil = createLayer(topo, self.xtra, self.layerthx[1], 6)      
        self.topo, topsoil = createLayer(topo, self.xtra, self.layerthx[0], 5)

        self.geom = bedrock+weathered+subsoil+topsoil
        
        self.mesh = createMesh(self.geom, self.topo, self.Q, area=self.area)

        self.mesh3D = pg.meshtools.extrudeMesh(self.mesh, a=np.array([0,zthx]))
        
    def show_geom(self):
        pg.show(self.geom)
        return 
    def show_mesh(self):
        pg.show(self.mesh)
        return 
    def save_mesh_vtk(self, fname):
        self.mesh3D.exportVTK(fname)
        return

# ------------------------------------------
if __name__ == "__main__":

    PH = PHMesh(xpos=300,dip=15,Q=14, xtra=1000)
    #PH.show_geom()
    PH.run_full(showPlot=True)






