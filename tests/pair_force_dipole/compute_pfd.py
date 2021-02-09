import numpy as np
import math
import matplotlib.pyplot as plt
from lammpstools import DumpLoader



def PBC(r1,r2,Lx,Ly,Lz):

    x1 = r1[0]
    y1 = r1[1]
    z1 = r1[2]

    x2 = r2[0]
    y2 = r2[1]
    z2 = r2[2]

    dx = x2-x1
    dy = y2-y1
    dz = z2-z1    
    while (dx > Lx/2.0):
        dx -= Lx
    while (dx < -Lx/2.0):
        dx += Lx

    while (dy > Ly/2.0):
        dy -= Ly
    while (dy < -Ly/2.0):
        dy += Ly
        
    while (dz > Lz/2.0):
        dz -= Lz
    while (dz < -Lz/2.0):
        dz += Lz
        
    return (dx,dy,dz)


def lj_fpair(r,epsilon=1.0):

    return -24*epsilon*(2/r**14-1/r**8)

def over_cut(r1,r2,Lx,Ly,Lz,cutoff=1.1224):

    
    dxs = []
    dys = []
    dzs = []
    drs = []

    dx,dy,dz = PBC(r1,r2,Lx,Ly,Lz)

    dxs.append(dx)
    dys.append(dy)
    dzs.append(dz)
    drs.append(np.sqrt(dx*dx+dy*dy+dz*dz))

    nx = math.floor(2*cutoff/Lx)+1
    ny = math.floor(2*cutoff/Ly)+1
    nz = math.floor(2*cutoff/Lz)+1

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if ( i== 0 and j == 0 and k == 0):
                    continue
                dx1 = dxs[0]+i*Lx
                dx2 = dxs[0]-i*Lx
                dy1 = dys[0]+j*Ly
                dy2 = dys[0]-j*Ly
                dz1 = dzs[0]+k*Lz
                dz2 = dzs[0]-k*Lz

                dr1 = np.sqrt(dx1**2+dy1**2+dz1**2)
                dr2 = np.sqrt(dx2**2+dy2**2+dz2**2)

                if dr1 < cutoff:
                    dxs.append(dx1)
                    dys.append(dy1)
                    dzs.append(dz1)
                    drs.append(dr1)

                if dr2 < cutoff:
                    dxs.append(dx2)
                    dys.append(dy2)
                    dzs.append(dz2)
                    drs.append(dr2)

    return dxs,dys,dzs,drs

def distance_matrices(xs,ys,zs,muxs,muys,muzs,
                      Lx,Ly,Lz):


    N = len(xs)
    dxs = []
    dys = []
    dzs = []
    drs = []
    dmuxs = []
    dmuys = []
    dmuzs = []
    atom_ids = []
    for i in range(N):
        for j in range(N):

            dx,dy,dz,dr = over_cut((xs[i],ys[i],zs[i]),
                                   (xs[j],ys[j],zs[j]),
                                   Lx,Ly,Lz)

            dxs = dxs + dx
            dys = dys + dy
            dzs = dzs + dz
            drs = drs + dr
            dmuxs = dmuxs + [muxs[i]-muxs[j]]
            dmuys = dmuys + [muys[i]-muys[j]]
            dmuzs = dmuzs + [muzs[i]-muzs[j]]
            atom_ids = atom_ids + [[i,j]]*len(dx)

            
    return dxs, dys, dzs, drs,dmuxs,dmuys,dmuzs, atom_ids


def fdote(xs,ys,zs,muxs,muys,muzs,Lx,Ly,Lz,d=3):

    dxs,dys,dzs,drs,dmuxs,dmuys,dmuzs,atom_ids = distance_matrices(xs,ys,zs,
                                                                   muxs,muys,muzs,
                                                                   Lx,Ly,Lz)

    ss = 0.0
    for i in range(len(dxs)):
        if np.abs(drs[i])>1e-12:

            ss += (dxs[i]*dmuxs[i]+dys[i]*dmuys[i]+dzs[i]*dmuzs[i])*lj_fpair(drs[i])

    if d == 3:
        vol = Lx*Ly*Lz
    else:
        vol = Lx*Ly

    ss *= 1.0/(2*d*(d-1)*vol)
    return ss



if __name__=="__main__":

    dumpname = 'dump.lammpstrj'

    dl = DumpLoader(dumpname,pstatus=True)



    Lz = dl.data[0]['z_size']
    Lx,Ly = dl.data[0]['box']

    muxs = dl.data[0]['mux']
    muys = dl.data[0]['muy']
    muzs = dl.data[0]['muz']


    xs = dl.data[0]['x']
    ys = dl.data[0]['y']
    zs = dl.data[0]['z']


    ss = fdote(xs,ys,zs,muxs,muys,muzs,Lx,Ly,Lz,d=3)

    print(ss)
