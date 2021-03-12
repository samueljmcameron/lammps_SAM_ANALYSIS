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

#include "compute_threebody_bessel.h"
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

static double bessk1( double x );
static double Vp_divide_r(double factorlj, double rij);

/* ---------------------------------------------------------------------- */

ComputeThreeBodyBessel::ComputeThreeBodyBessel(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  hist(nullptr), histall(nullptr)
{
  if (narg < 5)
    error->all(FLERR,"Illegal compute threebody/bessel command");

  array_flag = 1;
  extarray = 0;

  npos_bins = utils::inumeric(FLERR,arg[3],false,lmp);
  if (npos_bins < 1) error->all(FLERR,"Illegal compute threebody/bessel command");


  double Pi;
  Pi = utils::numeric(FLERR,arg[4],false,lmp);
  if (Pi < 0) error->all(FLERR,"Illegal compute threebody/bessel command");
  
  sqrtPi2fac = sqrt(Pi/2.0);
  

  // optional args
  // nargpair = # of pairwise args, starting at iarg = 5

  cutflag = 0;


  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"cutoff") == 0) {
      cutoff_user = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (cutoff_user <= 0.0) cutflag = 0;
      else cutflag = 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal compute threebody/bessel command");
  }

  if (force->newton_pair) {
    error->all(FLERR,"compute threebody/bessel command is incompatible with "
	       "newton pair being on.");
  }

  // pairwise args, fix it to be always one set of pairs for now,
  // but might change this later to deal with density differences


  int ntypes = atom->ntypes;
  if (ntypes != 1) {
    error->all(FLERR,"Cannot compute threebody/bessel for system with multiple "
	       "atom types.");
  }


  
  size_array_rows = npos_bins;
  size_array_cols = 5;
  
  memory->create(hist,4,npos_bins,"rdf:hist");
  memory->create(histall,4,npos_bins,"rdf:histall");    

  memory->create(array,size_array_rows,size_array_cols,"rdf:array");  

  dynamic = 0;
  natoms_old = 0;
}

/* ---------------------------------------------------------------------- */

ComputeThreeBodyBessel::~ComputeThreeBodyBessel()
{
  memory->destroy(hist);
  memory->destroy(histall);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeThreeBodyBessel::init()
{

  if (!force->pair && !cutflag)
    error->all(FLERR,"Compute threebody/bessel requires a pair style be defined "
               "or cutoff specified");

  if (cutflag) {
    double skin = neighbor->skin;
    mycutneigh = cutoff_user + skin;

    double cutghost;            // as computed by Neighbor and Comm
    if (force->pair)
      cutghost = MAX(force->pair->cutforce+skin,comm->cutghostuser);
    else
      cutghost = comm->cutghostuser;

    if (mycutneigh > cutghost)
      error->all(FLERR,"Compute threebody/bessel cutoff exceeds ghost atom range - "
                 "use comm_modify cutoff command");
    if (force->pair && mycutneigh < force->pair->cutforce + skin)
      if (comm->me == 0)
        error->warning(FLERR,"Compute threebody cutoff less than neighbor "
		       "cutoff - "
                       "forcing a needless neighbor list build");

    
    deldist = cutoff_user / npos_bins;

  } else deldist = force->pair->cutforce / npos_bins;

  deldistinv = 1.0/deldist;


  for (int i = 0; i < npos_bins; i++)
    array[i][0] = (i+0.5) * deldist;

  // initialize normalization, finite size correction, and changing atom counts
  // ...not sure what this does

  natoms_old = atom->natoms;
  dynamic = group->dynamic[igroup];
  if (dynamic_user) dynamic = 1;
  init_norm();

  // need an occasional half neighbor list
  // if user specified, request a cutoff = cutoff_user + skin
  // skin is included b/c Neighbor uses this value similar
  //   to its cutneighmax = force cutoff + skin
  // also, this NeighList may be used by this compute for multiple steps
  //   (until next reneighbor), so it needs to contain atoms further
  //   than cutoff_user apart, just like a normal neighbor list does

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
  if (cutflag) {
    neighbor->requests[irequest]->cut = 1;
    neighbor->requests[irequest]->cutoff = mycutneigh;
  }
  // need full neighbor list
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  
}

/* ---------------------------------------------------------------------- */

void ComputeThreeBodyBessel::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeThreeBodyBessel::init_norm()
{
  int i,j,m;

  // count atoms of each type that are also in group

  const int nlocal = atom->nlocal;
  const int ntypes = atom->ntypes;
  const int * const mask = atom->mask;

  //here
  typecount = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) typecount++;



  icount = typecount;

  jcount = typecount;

  duplicates = typecount;


  int scratch;
  MPI_Allreduce(&icount,&scratch,1,MPI_INT,MPI_SUM,world);
  icount = scratch;
  MPI_Allreduce(&jcount,&scratch,1,MPI_INT,MPI_SUM,world);
  jcount = scratch;
  MPI_Allreduce(&duplicates,&scratch,1,MPI_INT,MPI_SUM,world);
  duplicates = scratch;
  

}

/* ---------------------------------------------------------------------- */

void ComputeThreeBodyBessel::compute_array()
// ======================================================================
// Only counting triplets where all three distances are less than cutoff
// away.
//
// ======================================================================
{
  int i,j,m,ii,jj,inum,jnum,itype,jtype,ipair,jpair;
  int ij_bin,ik_bin,alpha_bin,dum_jk_bin;
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

  for (i = 0; i < npos_bins; i++) {
    hist[0][i] = 0;
    hist[1][i] = 0;
    hist[2][i] = 0;
    hist[3][i] = 0;
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

	
      xij = x[j][0]-xtmp;
      yij = x[j][1]-ytmp;
      zij = x[j][2]-ztmp;
      rij = sqrt(xij*xij+yij*yij+zij*zij);

      
      hist[0][ij_bin] += bessk1(sqrtPi2fac*rij);
      hist[1][ij_bin] += bessk1(sqrtPi2fac*rij)*Vp_divide_r(factor_lj,rij)*rij;

      for (kk = 0; kk < jnum; kk++) {
	if (jj == kk) continue;
	k = jlist[kk];
	factor_lj = special_lj[sbmask(k)];
	factor_coul = special_coul[sbmask(k)];
	
	k &= NEIGHMASK;


	if (factor_lj == 0.0 && factor_coul == 0.0) continue;

	if (!(mask[k] & groupbit)) continue;



	xik = x[k][0]-xtmp;
	yik = x[k][1]-ytmp;
	zik = x[k][2]-ztmp;


	rik = sqrt(xik*xik+yik*yik+zik*zik);

	rjk = sqrt((xik-xij)*(xik-xij)
		   +(yik-yij)*(yik-yij)+(zik-zij)*(zik-zij));
	
	dum_cosalpha = (xij*xik + yij*yik + zij*zik)/(rij*rik);

	ij_bin = static_cast<int> (rij*deldistinv);
	ik_bin = static_cast<int> (rik*deldistinv);
	dum_jk_bin = static_cast<int> (rjk*deldistinv);


	// have to add in this check in case rij and rik are
	// perfectly collinear, and then round off error might
	// make |dumcostheta| slightly greater than one
	if (dum_cosalpha > 1.0) {
	  dum_cosalpha = 1.0;

	} else if (dum_cosalpha < -1.0) {
	  dum_cosalpha = -1.0;
	}
	

	// check if the configuration can be counted consistently
	// (i.e. all three configs will be counted independent of
	// the labels i,j,k).
	if (ij_bin >= npos_bins || ik_bin >=npos_bins
	    || dum_jk_bin >= npos_bins) continue;

	hist[2][ij_bin] += bessk1(sqrtPi2fac*rik)*dum_cosalpha;
	hist[3][ij_bin] += bessk1(sqrtPi2fac*rik)*dum_cosalpha*Vp_divide_r(factor_lj,rij)*rij;


      }

    }
  }

  // sum histograms across procs

  MPI_Allreduce(hist[0],histall[0],4*npos_bins,MPI_DOUBLE,MPI_SUM,world);

  // convert counts to k(r) and coord(r) and copy into output array
  // vfrac = fraction of volume in shell m
  // npairs = number of pairs, corrected for duplicates
  // duplicates = pairs in which both atoms are the same

  double constant,vfrac,twobod_tmp,threebod_tmp;
  double rlower,rupper,normfac;
  double twobod_ncoord, threebod_ncoord;
  double rhoN;


  if (domain -> dimension == 2) {
    constant = MY_PI / (domain->xprd*domain->yprd);




    normfac = (icount > 0) ? static_cast<double>(jcount)
                - static_cast<double>(duplicates)/icount : 0.0;

    rhoN = icount*normfac/(domain->xprd*domain->yprd);

    
    twobod_ncoord = 0.0;
    threebod_ncoord = 0.0;
    for (ij_bin = 0; ij_bin < npos_bins; ij_bin++) {
      rlower = ij_bin*deldist;
      rupper = (ij_bin+1)*deldist;
      vfrac = constant * (rupper*rupper - rlower*rlower);
      if (vfrac * normfac != 0.0) {
	
	twobod_tmp = histall[0][ij_bin] / (vfrac * normfac * icount);
	threebod_tmp = histall[2][ij_bin] / (vfrac * normfac * icount);
	
      } else {
	twobod_tmp = 0.0;
	threebod_tmp = 0.0;
      }

      if (icount != 0) {
	twobod_ncoord += histall[1][ij_bin]/rhoN;
	threebod_ncoord += histall[3][ij_bin]/rhoN;
      }
      array[ij_bin][1] = twobod_tmp;
      array[ij_bin][2] = twobod_ncoord;
      array[ij_bin][3] = threebod_tmp;
      array[ij_bin][4] = threebod_ncoord;
    }
  } else {
    error->all(FLERR,"Cannot use compute threebody/bessel in a 3D system.");
  }
    
}


static double bessi1( double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
   } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (exp(ax)/sqrt(ax));
   }
   return x < 0.0 ? -ans : ans;
}


static double bessk1( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double y,ans;

   if (x <= 2.0) {
      y=x*x/4.0;
      ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
         +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
         +y*(-0.110404e-2+y*(-0.4686e-4)))))));
   } else {
      y=2.0/x;
      ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
         +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
         +y*(0.325614e-2+y*(-0.68245e-3)))))));
   }
   return ans;
}


static double Vp_divide_r(double factorlj, double rij)
{
  if (rij > pow(2,1./6.)) return 0;
  double r2inv = 1/(rij*rij);
  double r6inv = r2inv*r2inv*r2inv;
  return -24*r6inv*r2inv*(2*r6inv - 1);
}
