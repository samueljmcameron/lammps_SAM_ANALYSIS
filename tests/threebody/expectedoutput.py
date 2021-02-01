import numpy as np
from lammpstools import ReadData


def flat_index(ubin,vbin,abin,Npos_bins,Nangle_bins,nskip):
    v_r =  vbin-nskip

    u_r = ubin-nskip
    N0 = Npos_bins-2*nskip+1

    return abin+ (v_r + u_r*N0 - (u_r*(u_r+1))//2)*Nangle_bins


def v_3D(ulower,uupper,vlower,vupper,alower,aupper):

    v0 = (uupper**3-ulower**3)/3.0
    v0 *= (vupper**3-vlower**3)/3.0
    v0 *= (np.cos(alower) - np.cos(aupper))

    return v0

def v_2D(ulower,uupper,vlower,vupper,alower,aupper):

    v0 = (uupper**2-ulower**2)/2.0
    v0 *= (vupper**2-vlower**2)/2.0
    v0 *= (aupper-alower)

    return v0


def compute_threebody(xs,ys,zs,Natoms,npos_bins,nangle_bins,nskip,rc,
                      xL,yL,zL,dimension):
    dr = rc/npos_bins
    dalpha = np.pi/nangle_bins
    hists = np.zeros([nangle_bins,npos_bins,npos_bins],float)                      
    for i in range(Natoms):
        tmpx = xs[i]
        tmpy = ys[i]
        tmpz = zs[i]
        for j in range(Natoms):
            if (i == j):
                continue
            xij = xs[j]-tmpx
            yij = ys[j]-tmpy
            zij = zs[j]-tmpz

            for k in range(Natoms):
                if (k == j) or (k==i):
                    continue

                xik = xs[k]-tmpx
                yik = ys[k]-tmpy
                zik = zs[k]-tmpz

                xjk = xik - xij
                yjk = yik - yij
                zjk = zik - zij

                rij = np.sqrt(xij**2 + yij**2 + zij**2)
                rik = np.sqrt(xik**2 + yik**2 + zik**2)
                rjk = np.sqrt(xjk**2 + yjk**2 + zjk**2)


                ij_bin = int(rij/dr)
                ik_bin = int(rik/dr)
                jk_bin = int(rjk/dr)

                if (ij_bin >= npos_bins or ik_bin >= npos_bins
                    or jk_bin >= npos_bins):
                    continue

                rij_dot_rik = xij*xik + yij*yik + zij*zik
                alpha = np.arccos(rij_dot_rik/(rij*rik))



                alpha_bin = int(alpha/dalpha)
                hists[alpha_bin,ij_bin,ik_bin] += 1.0


    totbins = (npos_bins-2*nskip)*(npos_bins-2*nskip+1)*nangle_bins//2
    output = np.empty([totbins,Natoms],float)



    if dimension == 2:
        vol = xL*yL
    else:
        vol = xL*yL*zL



    const = vol*vol/(Natoms*Natoms*Natoms)

    for ubin in range(nskip,npos_bins-nskip):
        for vbin in range(nskip,npos_bins-ubin):
            for abin in range(nangle_bins):

                ulower = ubin*dr
                uupper = (ubin+1)*dr

                vlower = vbin*dr
                vupper = (vbin+1)*dr

                alower = abin*dalpha
                aupper = (abin+1)*dalpha


                if dimension == 2:
                    vf = v_2D(ulower,uupper,vlower,vupper,
                              alower,aupper)*4*np.pi/const
                else:
                    vf = v_3D(ulower,uupper,vlower,vupper,
                              alower,aupper)*8*np.pi**2/const


                fl = flat_index(ubin,vbin,abin,npos_bins,nangle_bins,nskip)
                output[fl,0] = int(fl+1)
                output[fl,1] = (ulower+uupper)/2.0
                output[fl,2] = (vlower+vupper)/2.0
                output[fl,3] = (alower+aupper)/2.0
                output[fl,4] =  hists[abin,vbin,ubin]/vf

    return output


if __name__ == "__main__":

    def check(dimension):

        RD = ReadData(f"data{dimension}d.3bod","atomic")


        headers,atoms = RD.headers, RD.atoms
        npos_bins = 15
        nangle_bins = 10
        nskip = 3
        rc = 1.3


        Atoms = atoms['Atoms']


        Natoms = headers["atoms"]

        xs = Atoms['x']
        ys = Atoms['y']
        zs = Atoms['z']

        xL = headers['xlo xhi'][1]-headers['xlo xhi'][0]
        yL = headers['ylo yhi'][1]-headers['ylo yhi'][0]
        zL = headers['zlo zhi'][1]-headers['zlo zhi'][0]


        output = compute_threebody(xs,ys,zs,Natoms,npos_bins,nangle_bins,nskip,rc,
                                   xL,yL,zL,dimension)            

        lmpout = np.loadtxt(f"out{dimension}d.3bod",skiprows=4)

        if np.allclose(output,lmpout):
            print(f'{dimension}D computations match. great success!')
        else:
            print(f'{dimension}D computations do not match.')




    check(2)
    check(3)
