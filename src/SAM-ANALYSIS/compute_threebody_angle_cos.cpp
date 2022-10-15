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
   Contributing author: Sam Cameron
------------------------------------------------------------------------- */

#include "compute_threebody_angle_cos.h"
#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace MathConst;


static double compute_3Dvf(double ulower, double uupper,
			   double vlower, double vupper,
			   double delalpha, int alpha_bin);

static double compute_2Dvf(double ulower, double uupper,
			   double vlower, double vupper,
			   double delalpha, int alpha_bin );


/* ---------------------------------------------------------------------- */

ComputeThreeBodyAngleCos::ComputeThreeBodyAngleCos(LAMMPS *lmp, int narg, char **arg) :
  ComputeThreeBody(lmp, narg, arg,false)
{
  
  size_array_rows = (npos_bins-2*nskip)*(npos_bins-2*nskip+1)/2;
  size_array_cols = 4;

  memory->create(array,size_array_rows,size_array_cols,"rdf:array");  
  
  dynamic = 0;
  natoms_old = 0;
}

/* ---------------------------------------------------------------------- */

ComputeThreeBodyAngleCos::~ComputeThreeBodyAngleCos()
{
  memory->destroy(hist);
  memory->destroy(histall);
  memory->destroy(array);
}


/* ---------------------------------------------------------------------- */

void ComputeThreeBodyAngleCos::compute_array()
// ======================================================================
// Only counting triplets where all three distances are less than cutoff
// away.
//
// ======================================================================
{
  int i,j,m,ii,jj,inum,jnum,itype,jtype,ipair,jpair;
  int ij_bin,ik_bin,alpha_bin,dum_jk_bin;
  int flat;
  int k,kk,knum,ktype,kpair;
  double xtmp,ytmp,ztmp;
  double xij,yij,zij,xik,yik,zik,rij,rik,rjk,alpha;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double dum_cosalpha;
  double factor_lj,factor_coul;

  double denom;
  
  if (natoms_old != atom->natoms) {
    dynamic = 1;
    natoms_old = atom->natoms;
  }

  // if the number of atoms has changed or we have a dynamic group
  // or dynamic updates are requested (e.g. when changing atom types)
  // we need to recompute some normalization parameters

  if (dynamic) init_norm();

  invoked_array = update->ntimestep;

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero the histogram counts

  for (ij_bin = nskip; ij_bin < (npos_bins-nskip); ij_bin++) 
    for (ik_bin = nskip; ik_bin < npos_bins-ij_bin; ik_bin++) 
      for (alpha_bin = 0; alpha_bin < nangle_bins; alpha_bin ++ ) {

	flat = flat_index(ij_bin, ik_bin, alpha_bin, npos_bins,
			  nangle_bins,nskip);
	hist[flat] = 0;
      }


  // tally the three body
  // both atom i and j must be in fix group
  // itype,jtype must have been specified by user
  // consider I,J as one interaction even if neighbor pair is stored on 2 procs
  // tally I,J pair each time I is central atom, and each time J is central

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;

  // loop over full neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      // if both weighting factors are 0, skip this pair
      // could be 0 and still be in neigh list for long-range Coulombics
      // want consistency with non-charged pairs which wouldn't be in list

      if (factor_lj == 0.0 && factor_coul == 0.0) continue;



      if (!(mask[j] & groupbit)) continue;

      for (kk = 0; kk < jnum; kk++) {
	if (jj == kk) continue;
	k = jlist[kk];
	factor_lj = special_lj[sbmask(k)];
	factor_coul = special_coul[sbmask(k)];
	
	k &= NEIGHMASK;


	if (factor_lj == 0.0 && factor_coul == 0.0) continue;

	if (!(mask[k] & groupbit)) continue;

	
	xij = x[j][0]-xtmp;
	yij = x[j][1]-ytmp;
	zij = x[j][2]-ztmp;

	xik = x[k][0]-xtmp;
	yik = x[k][1]-ytmp;
	zik = x[k][2]-ztmp;

	rij = sqrt(xij*xij+yij*yij+zij*zij);
	rik = sqrt(xik*xik+yik*yik+zik*zik);

	rjk = sqrt((xik-xij)*(xik-xij)
		   +(yik-yij)*(yik-yij)+(zik-zij)*(zik-zij));
	
	dum_cosalpha = (xij*xik + yij*yik + zij*zik)/(rij*rik);

	ij_bin = static_cast<int> (rij*deldistinv);
	ik_bin = static_cast<int> (rik*deldistinv);
	dum_jk_bin = static_cast<int> (rjk*deldistinv);


	// have to add in this check in case r_ij and r_ik are
	// perfectly collinear, and then round off error might
	// make |dumcostheta| slightly greater than one
	if (dum_cosalpha > 1.0) {
	  dum_cosalpha = 1.0;

	} else if (dum_cosalpha < -1.0) {
	  dum_cosalpha = -1.0;
	}
	
	alpha = acos(dum_cosalpha);


	alpha_bin = static_cast<int> (alpha*delalphainv);	

	// when dum_cosalpha = -1 necessary to check for this
	if (alpha_bin == nangle_bins) {
	  alpha_bin = nangle_bins -1;
	}
	
	if (ij_bin >= npos_bins || ik_bin >=npos_bins
	    || dum_jk_bin >= npos_bins) continue;
	if (alpha_bin >= nangle_bins) {
	  printf("alpha = %f, alpha_bin = %d\n",alpha,alpha_bin);
	  error->one(FLERR,"theta > 3.1415 somehow? "
		     "Error in compute threebody.");
	}

	flat = flat_index(ij_bin, ik_bin, alpha_bin, npos_bins,
			  nangle_bins,nskip);
	hist[flat] += 1.0;

      }
    }
  }

  // sum histograms across procs

  MPI_Allreduce(hist,histall,nbin_total,MPI_DOUBLE,MPI_SUM,world);

  // convert counts to g(r) and coord(r) and copy into output array
  // vfrac = fraction of volume in shell m
  // npairs = number of pairs, corrected for duplicates
  // duplicates = pairs in which both atoms are the same

  
  double constant,vfrac,gr,ulower,uupper,vlower,vupper,normfac;
  int flat_index;
  
  if (domain->dimension==2) {
    constant = (domain->xprd*domain->yprd);
    constant = constant*constant;
    
    // the 2 in front of the delalpha below is necessary because
    // delalpha = pi/nangle_bins over the [0,pi] region, but really
    // the g_3bod is defined from 0 to 2pi. So we are double-counting
    // in alpha since we are only letting theta [0,pi].
    constant = 2*2*MY_PI/constant;
    
    normfac = (icount > 0) ? static_cast<double>(icount) : 0.0;
    normfac = normfac*normfac*normfac;

    
    set_array(constant, normfac,compute_2Dvf);
    
  } else  {
    
    constant = domain->xprd*domain->yprd*domain->zprd;
    constant = constant*constant;

    constant = 8*MY_PI*MY_PI/constant;

    normfac = (icount > 0) ? static_cast<double>(icount) : 0.0;

    normfac = normfac*normfac*normfac;
      
    set_array(constant, normfac,compute_3Dvf);
  }

}

void ComputeThreeBodyAngleCos::set_array(double constant, double normfac,
					 double (*vf)(double,double,double,
						      double,double,int))
{
  double ulower,uupper,vlower,vupper;
  double vfrac;

  int alpha_bin,ij_bin,ik_bin;
  double gr,gr_cos;
  
  int flat;
  

    
  for (ij_bin = nskip; ij_bin < (npos_bins-nskip); ij_bin++) {
    
    ulower = ij_bin*deldist;
    uupper = (ij_bin+1)*deldist;
    
    for (ik_bin = nskip; ik_bin < npos_bins-ij_bin; ik_bin++) {
      
      vlower = ik_bin*deldist;
      vupper = (ik_bin+1)*deldist;	  

      gr = 0.0;
      gr_cos = 0.0;
      
      for (alpha_bin = 0; alpha_bin < nangle_bins; alpha_bin ++ ) {

	
	vfrac = constant * vf(ulower,uupper,vlower,vupper,delalpha,
			      alpha_bin);

	flat = flat_index(ij_bin, ik_bin, alpha_bin, npos_bins,
			  nangle_bins,nskip);
	if (vfrac * normfac != 0.0) {
	  gr += histall[flat]/(vfrac *normfac)*delalpha;
	  gr_cos += (histall[flat]/(vfrac *normfac)
		     *cos((alpha_bin + 0.5)*delalpha));
	} else {
	  gr += 0.0;
	  gr_cos += 0.0;
	}

	flat = flat_index_no_theta(ij_bin, ik_bin, npos_bins,nskip);
	array[flat][0] = (ij_bin + 0.5)*deldist;
	array[flat][1] = (ik_bin + 0.5)*deldist;
	array[flat][2] = gr;
	array[flat][3] = gr_cos;
      }
    }
  }

  return;
}


int ComputeThreeBodyAngleCos::flat_index_no_theta(int ubin, int vbin,
						  int Nuv,int nsk)
{

  int udum = ubin - nsk;
  return (vbin-nsk + udum*(Nuv-2*nsk+1)
	  - (udum * (udum+1) ) / 2);
}


double compute_3Dvf(double ulower, double uupper,
		    double vlower, double vupper,
		    double delalpha,int alpha_bin)

{
  double alower,aupper,vf;
  
  alower = cos((alpha_bin+1)*delalpha);
  aupper = cos(alpha_bin*delalpha); 
   
  return ((uupper*uupper*uupper - ulower*ulower*ulower)/3.0
	  * (vupper*vupper*vupper - vlower*vlower*vlower)/3.0
	  * (aupper-alower)
	  / ( sin((alpha_bin + 0.5)*delalpha)*delalpha));
}

double compute_2Dvf(double ulower, double uupper,
		    double vlower, double vupper,
		    double delalpha,int  alpha_bin )
{

  return ( (uupper*uupper - ulower*ulower)/2.0
	   * (vupper*vupper - vlower*vlower)/2.0
	   /(  2));
}



