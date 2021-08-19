
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
DEFAULTPATH = fspath(pathlib.Path(__file__).parent.absolute())

# ----------Create Topo Polygon ----------------------
def createBackground(efile, xextra, botdep):

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
    
class UniformMesh:
    def __init__(self, efile=None, dep=200, xtra=500, Q=20, deplist=[2.5], zthx=5, outfile=None, area=None):
        self.dep = dep
        self.Q = Q
        self.outfile = outfile
        self.xtra = xtra
        self.area = area
        self.zthx = zthx
        self.layerz = -1 * np.array(deplist)
        if efile is None:
            self.efile = DEFAULTPATH+"/data/PH-2018-eloc.txt"
        else:
            self.efile = efile

            
        self.topo, bpoly = createBackground(self.efile, self.xtra, self.dep)        
        self.geom = bpoly
        
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






