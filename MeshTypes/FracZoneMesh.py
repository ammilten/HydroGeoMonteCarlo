
# -----------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate 

import pygimli as pg
import pygimli.meshtools as mt
import pygimli.physics.ert as ert

import pathlib
from os import fspath
DEFAULTPATH = fspath(pathlib.Path(__file__).parent.parent.absolute())

# ----------Create Topo Polygon ----------------------
def createBackground(efile, xextra, botdep):
    # eloc = pd.read_csv(efile,sep='\t',header=None)
    # nelec = eloc.values.shape[0]
    # topo = eloc.values[:,(1,3)]
    
    topo = np.genfromtxt(efile ,delimiter=',',skip_header=True)
    
    e_bot = np.floor(np.min(topo[:,1]))
    xL = np.floor(np.min(topo[:,0]))
    xR = np.ceil(np.max(topo[:,0]))

    L_ext = np.array([[xL-xextra, topo[0,1]]])
    R_ext = np.array([[xR+xextra, topo[topo.shape[0]-1,1]]])
        
    bots = np.array([[R_ext[0,0], e_bot-botdep],[L_ext[0,0],e_bot-botdep]])
    combo = np.concatenate((topo,R_ext,bots,L_ext))

    bpoly = mt.createPolygon(combo, isClosed=True, marker=2)

    return topo, bpoly


# --------- Create Fault Polygon ------------
def createFault(topo, xpos, dip, H, dep, xtra=500):
    topo2 = np.concatenate((np.array([[topo[0,0]-xtra, topo[0,1]]]), topo, np.array([[topo[topo.shape[0]-1,0]+xtra,topo[topo.shape[0]-1,1]]])))
    # topo2[:,1] = topo2[:,1]-shdep
    
    zbot = np.floor(np.min(topo2[:,1])) - dep

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

    middles = [(topo[j,0],topo[j,1]) for j in range(idxabove-1,idxbelow)]
    verts = [LL, UL] + middles + [UR,LR]
    fpoly = mt.createPolygon(verts, isClosed=True, addNodes=0, marker=1)

    return fpoly


# ---------- Create mesh ---------------------
def createMesh(geom, topo, layerz, Q, area=None):
    
    for j in range(layerz.shape[0]):
        for i in range(topo.shape[0]):
            geom.createNode(topo[i,:] + [0.0, layerz[j]])

    if area is None:
        mesh = mt.createMesh(geom, quality=Q)
    else:
        mesh = mt.createMesh(geom, quality=Q, area=area)

    return mesh
    
class FZMesh:
    def __init__(self, dip=60, H=100, xpos=250, efile=None, dep=200, xtra=500, Q=20, zthx=5, deplist=[0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0 ], outfile=None, area=None):
        self.dep = dep
        self.dip = dip
        self.H = H
        self.xpos = xpos
        self.Q = Q
        self.outfile = outfile
        self.xtra = xtra
        self.area = area
        self.zthx = zthx
        self.layerz = -1 * np.array(deplist)
        if efile is None:
            self.efile = DEFAULTPATH+"/data/topo.csv"
        else:
            self.efile = efile

            
        self.topo, bpoly = createBackground(self.efile, self.xtra, self.dep)        
        fpoly = createFault(self.topo, self.xpos, self.dip, self.H, self.dep, xtra=self.xtra)
        self.geom = bpoly+fpoly
        
        self.mesh = createMesh(self.geom, self.topo, self.layerz, self.Q, area=self.area)
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






