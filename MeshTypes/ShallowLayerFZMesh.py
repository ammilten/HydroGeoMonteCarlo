
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

    bpoly = mt.createPolygon(combo, isClosed=True, marker=2)

    return topo, bpoly

# ----------Create Shallow Layer ----------------------
def createLayer(efile, xextra, shdep):
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

# --------- Create Fault Polygon ------------
def createFault(topo, xpos, dip, H, dep, shdep, xtra=500):
    topo2 = np.concatenate((np.array([[topo[0,0]-xtra, topo[0,1]]]), topo, np.array([[topo[topo.shape[0]-1,0]+xtra,topo[topo.shape[0]-1,1]]])))
    topo2[:,1] = topo2[:,1]-shdep
    
    zbot = np.floor(np.min(topo2[:,1])) - dep+shdep

    Z = interpolate.interp1d(topo2[:,0],topo2[:,1])
    zx = Z(xpos)
    Dz = zx-zbot

    Dx = Dz/np.tan(np.deg2rad(dip))
    xbot = xpos + Dx

    Dxbot = H/2/np.sin(np.deg2rad(dip))

    LR = (xbot+Dxbot,zbot)
    LL = (xbot-Dxbot,zbot)

    zfR = zbot + (LR[0]-topo2[:,0])*np.tan(np.deg2rad(dip))
    diffR = interpolate.interp1d(zfR-topo2[:,1],topo2[:,0])
    UR = (float(diffR(0)), float(Z(diffR(0))) )

    zfL = zbot + (LL[0]-topo2[:,0])*np.tan(np.deg2rad(dip))
    diffL = interpolate.interp1d(zfL-topo2[:,1],topo2[:,0])
    UL = (float(diffL(0)), float(Z(diffL(0))) )
    
    idxabove = (topo2[topo2[:,0] > UL[0],0]).argmin() + (topo2[topo2[:,0] < UL[0],0]).shape[0]
    idxbelow = (topo2[topo2[:,0] < UR[0],0]).argmax()

    middles = [(topo[j,0],topo[j,1] - shdep) for j in range(idxabove-1,idxbelow)]
    verts = [LL, UL] + middles + [UR,LR]
    fpoly = mt.createPolygon(verts, isClosed=True, addNodes=0, marker=1)

    return fpoly


# ---------- Create mesh ---------------------
def createMesh(geom, topo, Q, area=None):

    #scheme = pg.DataContainerERT()
    #scheme.setSensorPositions(topo)

    if area is None:
        mesh = mt.createMesh(geom, quality=Q)
    else:
        mesh = mt.createMesh(geom, quality=Q, area=area)

    return mesh
    
class PHMesh:
    def __init__(self, dip=60, H=100, xpos=250, shdep=5, efile=None, dep=200, xtra=500, Q=20, zthx=5, outfile=None, area=None):
        self.dep = dep
        self.dip = dip
        self.H = H
        self.xpos = xpos
        self.shdep = shdep
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
        _, bpoly = createBackground(self.efile, self.xtra, self.dep, self.shdep)   
#        print('Creating Surface Layer')  
        time.sleep(5)
        self.topo, lpoly = createLayer(self.efile, self.xtra, self.shdep)
#        print('Creating Fracture Zone')
        time.sleep(5)
        fpoly = createFault(self.topo, self.xpos, self.dip, self.H, self.dep, self.shdep, xtra=self.xtra)
#        print('Creating Geometry')
        time.sleep(5)
        self.geom = bpoly+lpoly+fpoly
        
#        print('Creating Mesh')
        time.sleep(5)
        self.mesh = createMesh(self.geom, self.topo, self.Q, area=self.area)
#        print('Extruding Mesh')
        time.sleep(5)
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






