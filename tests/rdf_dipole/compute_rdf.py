import numpy as np
import math
import matplotlib.pyplot as plt
from lammpstools.dumploader import DumpLoader
from lammpstools.histogramloader import HistogramLoader


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


def over_cut(r1,r2,Lx,Ly,Lz,cutoff=1.5):

    
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

def distance_matrices(xs,ys,zs,Lx,Ly,Lz):


    N = len(xs)
    dxs = []
    dys = []
    dzs = []
    drs = []
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
            atom_ids = atom_ids + [[i,j]]*len(dx)

            
    return dxs, dys, dzs, drs, atom_ids


def rdf(xs,ys,zs,muxs,muys,muzs,Lx,Ly,Lz,
        nbins=5,cutoff=1.5):

    dxs,dys,dzs,drs,atom_ids = distance_matrices(xs,ys,zs,Lx,Ly,Lz)
    
    N = len(drs)

    dr = 1.5/nbins

    rs = np.linspace(0,cutoff,num=nbins,endpoint=False)+0.5*dr

    hists = np.zeros([7,nbins],float)
    outputs = np.zeros([9,nbins],float)
    outputs[0] = rs

    for m in range(N):

        i,j = atom_ids[m]
        dxij = dxs[m]
        dyij = dys[m]
        dzij = dzs[m]
        rij = drs[m]

        if (i == j):
            continue
        
        ibin = int(math.floor(rij/dr))
        
        if ibin >= nbins:
            continue
        
        factor = 1.0

        delmux = muxs[j]-muxs[i]
        delmuy = muys[j]-muys[i]
        delmuz = muzs[j]-muzs[i]
        
        pmux = muxs[i] + muxs[j]
        pmuy = muys[i] + muys[j]
        pmuz = muzs[i] + muzs[j]
            
        udotu = muxs[i]*muxs[j]+muys[i]*muys[j]+muzs[i]*muzs[j]
        udotr = (delmux*dxij + delmuy*dyij + delmuz*dzij)/rij
        pdotr = (pmux*dxij+pmuy*dyij+pmuz*dzij)/rij
        hists[0,ibin] += factor
        hists[1,ibin] += factor*udotu
        hists[2,ibin] += factor*udotr
        hists[3,ibin] += factor*udotu*udotu
        hists[4,ibin] += factor*udotr*udotr
        hists[5,ibin] += factor*udotu*udotr
        hists[6,ibin] += factor*pdotr*pdotr


    constant = 4.0*np.pi/(3.0*Lx*Ly*Lz)

    ncoord = 0.0
    normfac = 31.0
    icount = 32.0
    for ibin in range(nbins):

        rlower = ibin*dr
        rupper = (ibin+1)*dr
        vfrac = constant*(rupper*rupper*rupper-rlower*rlower*rlower)

        if (vfrac*normfac != 0):

            denom = normfac*vfrac*icount
            gr = hists[0,ibin]/denom
            u_u = hists[1,ibin]/denom
            u_r = hists[2,ibin]/denom
            u2_u2 = hists[3,ibin]/denom
            u2_r2 = hists[4,ibin]/denom
            u_r_u_r = hists[5,ibin]/denom
            p2_r2 = hists[6,ibin]/denom

        else:
            
            gr = 0
            u_u = 0
            u_r = 0
            u2_u2 = 0
            u2_r2 = 0
            u_r_u_r = 0
            p2_r2 = 0

        ncoord += gr*vfrac*normfac

        outputs[1,ibin] = gr
        outputs[2,ibin] = u_u
        outputs[3,ibin] = u_r
        outputs[4,ibin] = u2_u2
        outputs[5,ibin] = u2_r2
        outputs[6,ibin] = u_r_u_r
        outputs[7,ibin] = p2_r2
        outputs[8,ibin] = ncoord
    
        
    return outputs



if __name__=="__main__":

    dumpname = 'dump.lammpstrj'
    histname = 'histos.rdf'

    dl = DumpLoader(dumpname,pstatus=True)
    hl = HistogramLoader(histname,pstatus=True)



    Lz = dl.data[0]['z_size']
    Lx,Ly = dl.data[0]['box']

    muxs = dl.data[0]['mux']
    muys = dl.data[0]['muy']
    muzs = dl.data[0]['muz']


    xs = dl.data[0]['x']
    ys = dl.data[0]['y']
    zs = dl.data[0]['z']


    hldat = hl.data[0]

    outputs = rdf(xs,ys,zs,muxs,muys,muzs,Lx,Ly,Lz)


    for i in range(9):
        name = f'c_myhistos[{i+1}]'

        bools = np.isclose(hldat[name],outputs[i],rtol=1e-3)
        if not np.any(bools):
            print('lammps computation is different than python computation!')
        else:
            print('lammps computations match. great success!')
