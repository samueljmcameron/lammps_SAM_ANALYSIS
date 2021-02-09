/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "compute_pair_force_dipole.h"


#include "atom.h"
#include "error.h"
#include "force.h"
#include "domain.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"


#include <cmath>
#include <cstring>


using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

ComputePairForceDipole::ComputePairForceDipole(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute pair/force/dipole command");

  scalar_flag = 1;
  extscalar = 0;

  // error check

  if (!atom->mu_flag)
    error->all(FLERR,"Compute pair/force/dipole requires atom attribute mu");

}



/* ---------------------------------------------------------------------- */
void ComputePairForceDipole::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"No pair style is defined for compute pair/force/dipole");
  if (force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support compute pair/force/dipole");

  // need an occasional half neighbor list
  // set size to same value as request made by force->pair
  // this should enable it to always be a copy list (e.g. for granular pstyle)

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
  NeighRequest *pairrequest = neighbor->find_request((void *) force->pair);
  if (pairrequest) neighbor->requests[irequest]->size = pairrequest->size;

}

/* ---------------------------------------------------------------------- */

void ComputePairForceDipole::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}


double ComputePairForceDipole::compute_scalar()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,fpair,factor_lj,factor_coul;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **mu = atom->mu;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  double dimension = domain->dimension;
  double inv_volume;

  if (dimension == 3) {
    inv_volume = 1.0/(domain->xprd * domain->yprd * domain->zprd);
  } else {
    inv_volume = 1.0/(domain->xprd * domain->yprd);
  }
  double prefac = inv_volume/(2*dimension*(dimension-1));

  double tmpout;
  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // scalar to store sum in
  double output = 0.0;

  
  // loop over neighbors of my atoms
  // skip if I or J are not in group
  // for newton = 0 and J = ghost atom,
  //   need to insure I,J pair is only output by one proc
  //   use same itag,jtag logic as in Neighbor::neigh_half_nsq()

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itag = tag[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      // itag = jtag is possible for long cutoffs that include images of self

      if (newton_pair == 0 && j >= nlocal) {
        jtag = tag[j];
        if (itag > jtag) {
          if ((itag+jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag+jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp) {
            if (x[j][1] < ytmp) continue;
            if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
          }
        }
      }


      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	// compute fpair
	pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);


	tmpout = fpair*(delx*(mu[i][0]-mu[j][0])+dely*(mu[i][1]-mu[j][1])
			+delz*(mu[i][2]-mu[j][2]));
	output += tmpout;
	
	if (newton_pair || j < nlocal) {
	  output += tmpout;
	}
      }
	  
    }
  }


  double nktv2p = force->nktv2p;

  // sum over all processors
  MPI_Allreduce(&output,&scalar,1,MPI_DOUBLE,MPI_SUM,world);

  scalar = scalar*prefac * nktv2p;
  return scalar;
}



