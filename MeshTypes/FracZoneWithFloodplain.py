
# -----------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate 

import pygimli as pg
import pygimli.meshtools as mt
import pygimli.physics.ert as ert

import sys
import time
import pathlib
from os import fspath
DEFAULTPATH = fspath(pathlib.Path(__file__).parent.parent.absolute())

# ----------Create Topo Polygon ----------------------
def createBackground(efile, xextra, botdep, shdep=0):
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
def createFPpoly(efile, fpdep, fplim):

    topo = np.genfromtxt(efile ,delimiter=',',skip_header=True)
    
    e_bot = np.floor(np.min(topo[:,1]))
    iL = np.flipud(np.argwhere(topo[:,0] < fplim[0]))[0][0]
    iR = np.argwhere(topo[:,0] > fplim[1])[0][0]

    topo_fp = topo[iL:iR,:]
    topo_bot = np.copy(np.flipud(topo_fp))
    topo_fp[:,1] = topo_fp[:,1]
    topo_bot[:,1] = topo_bot[:,1]-fpdep
    combo = np.concatenate((topo_fp, topo_bot))

    med = np.int(topo_fp.shape[0]/2) #Setting position for the marker
    markerpos = [topo_bot[med,0],topo_bot[med,1]+fpdep/2]
    
#    lpoly = mt.createPolygon(combo, isClosed=True, markerPosition=markerpos, marker=4)

    return combo, markerpos
    
def createFloodplain(verts, markerpos):
    return mt.createPolygon(verts, isClosed=True, markerPosition=markerpos, marker=4)

# --------- Create Fault Polygon ------------
def createFault(topo, fpverts, xpos, dip, H, dep, shdep=0, xtra=500):
    fpverts_new = np.copy(fpverts)
    
    topo2 = np.concatenate((np.array([[topo[0,0]-xtra, topo[0,1]]]), topo, np.array([[topo[topo.shape[0]-1,0]+xtra,topo[topo.shape[0]-1,1]]])))
    topo2[:,1] = topo2[:,1]-shdep
    
    zbot = np.floor(np.min(topo2[:,1])) - dep+shdep

    Z = interpolate.interp1d(topo2[:,0],topo2[:,1])
    zx = Z(xpos)
    Dz = zx-zbot

    Dx = Dz/np.tan(np.deg2rad(dip))
    xbot = xpos + Dx

    Dxbot = H/2/np.sin(np.deg2rad(dip))

    LR = np.array((xbot+Dxbot,zbot))[np.newaxis,:]
    LL = np.array((xbot-Dxbot,zbot))[np.newaxis,:]

    zfR = zbot + (LR[0,0]-topo2[:,0])*np.tan(np.deg2rad(dip))
    diffR = interpolate.interp1d(zfR-topo2[:,1],topo2[:,0])
    UR = np.array((float(diffR(0)), float(Z(diffR(0))) ))[np.newaxis,:]
    if UR[0,0] > np.min(fpverts[:,0]):
        midind = np.int(fpverts.shape[0]/2)
        botverts = np.flipud(fpverts[midind:,:])
        
        zfR = zbot + (LR[0,0]-botverts[:,0])*np.tan(np.deg2rad(dip))
        diffR = interpolate.interp1d(zfR-botverts[:,1],botverts[:,0])
        
        Zfp = interpolate.interp1d(botverts[:,0],botverts[:,1])
        
        UR = np.array((float(diffR(0)), float(Zfp(diffR(0))) ))
        botverts_new_L = np.flipud(botverts[botverts[:,0]<UR[0],:])
        botverts_new_R = np.flipud(botverts[botverts[:,0]>UR[0],:])
        botverts_new = np.concatenate((botverts_new_R, UR[np.newaxis,:], botverts_new_L))
        fpverts_new = np.concatenate((fpverts[:midind], botverts_new))
        
        inds_left_of_intersect = botverts[:,0] < UR[0]
        
        UR_fp = np.flipud(botverts[inds_left_of_intersect,:])
        UR_total = np.concatenate((UR[np.newaxis,:], UR_fp, fpverts[0,:][np.newaxis,:]))
        UR = np.flipud(UR_total)
        #sys.exit("Cannot handle fracture zone intersecting floodplain yet")

    zfL = zbot + (LL[0,0]-topo2[:,0])*np.tan(np.deg2rad(dip))
    diffL = interpolate.interp1d(zfL-topo2[:,1],topo2[:,0])
    UL = np.array((float(diffL(0)), float(Z(diffL(0))) ))[np.newaxis,:]    
    
    idxabove = (topo2[topo2[:,0] > UL[0,0],0]).argmin() + (topo2[topo2[:,0] < UL[0,0],0]).shape[0]
    idxbelow = (topo2[topo2[:,0] < UR[0,0],0]).argmax()

    middles = np.array([(topo[j,0],topo[j,1] - shdep) for j in range(idxabove-1,idxbelow)])

    #verts = [LL, UL] + middles + [UR,LR]
    verts = np.concatenate((LL,UL,middles,UR,LR))
    fpoly = mt.createPolygon(verts, isClosed=True, addNodes=0, marker=1)

    return fpoly, fpverts_new

# ---------- Create mesh ---------------------
def createMesh(geom, topo, Q, area=None):

    #scheme = pg.DataContainerERT()
    #scheme.setSensorPositions(topo)

    if area is None:
        mesh = mt.createMesh(geom, quality=Q)
    else:
        mesh = mt.createMesh(geom, quality=Q, area=area)

    return mesh
    
class FPFZMesh:
    def __init__(self, dip=60, H=100, xpos=250, fpdep=5, fplim=[450,650], efile=None, dep=200, xtra=500, Q=20, zthx=5, outfile=None, area=None):
        self.dep = dep
        self.dip = dip
        self.H = H
        self.xpos = xpos
        self.fpdep = fpdep
        self.fplim = fplim
        self.Q = Q
        self.outfile = outfile
        self.xtra = xtra
        self.area = area
        self.zthx = zthx
        if efile is None:
            self.efile = DEFAULTPATH+"/data/topo.csv"
        else:
            self.efile = efile

            
#        print('Creating Background')
        self.topo, bpoly = createBackground(self.efile, self.xtra, self.dep)   
#        print('Creating Surface Layer')  
        time.sleep(5)
        fpverts, markerpos = createFPpoly(self.efile, self.fpdep, self.fplim)
#        print('Creating Fracture Zone')
        time.sleep(5)
        fpoly,fpverts_new = createFault(self.topo, fpverts, self.xpos, self.dip, self.H, self.dep, xtra=self.xtra)
#        print('Creating Geometry')
        time.sleep(5)
        lpoly = createFloodplain(fpverts_new, markerpos)
        self.geom = bpoly+fpoly+lpoly
        
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






