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

#include "compute_threebody.h"
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
			   double /* delalpha */,
			   int /* alpha_bin */);


/* ---------------------------------------------------------------------- */

ComputeThreeBody::ComputeThreeBody(LAMMPS *lmp, int narg, char **arg,
				   bool alloc_array) :
  Compute(lmp, narg, arg),
  hist(nullptr), histall(nullptr)
{
  if (narg < 5)
    error->all(FLERR,"Illegal compute three body command");

  array_flag = 1;
  extarray = 0;

  npos_bins = utils::inumeric(FLERR,arg[3],false,lmp);
  if (npos_bins < 1) error->all(FLERR,"Illegal compute three body command");

  nangle_bins = utils::inumeric(FLERR,arg[4],false,lmp);
  if (nangle_bins < 1) error->all(FLERR,"Illegal compute three body command");


  // optional args
  // nargpair = # of pairwise args, starting at iarg = 5

  cutflag = 0;
  nskip = 0;


  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"cutoff") == 0) {
      cutoff_user = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (cutoff_user <= 0.0) cutflag = 0;
      else cutflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"skip") == 0) {
      nskip = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute threebody command");
  }

  if (force->newton_pair) {
    error->all(FLERR,"compute threebody command is incompatible with "
	       "newton pair being on.");
  }

  // pairwise args, fix it to be always one set of pairs for now,
  // but might change this later to deal with density differences

  delalpha = MY_PI/nangle_bins;

  delalphainv = 1.0/delalpha;


  int ntypes = atom->ntypes;
  if (ntypes != 1) {
    error->all(FLERR,"Cannot compute threebody for system with multiple "
	       "atom types.");
  }


  
  size_array_rows = (npos_bins-2*nskip)*(npos_bins-2*nskip+1)*nangle_bins/2;
  size_array_cols = 4;
  nbin_total = size_array_rows;

  memory->create(hist,nbin_total,"rdf:hist");
  memory->create(histall,nbin_total,"rdf:histall");    
    
  if (alloc_array) {

    memory->create(array,size_array_rows,size_array_cols,"rdf:array");  

  }
  dynamic = 0;
  natoms_old = 0;
}

/* ---------------------------------------------------------------------- */

ComputeThreeBody::~ComputeThreeBody()
{
  memory->destroy(hist);
  memory->destroy(histall);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeThreeBody::init()
{

  if (!force->pair && !cutflag)
    error->all(FLERR,"Compute threebody requires a pair style be defined "
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
      error->all(FLERR,"Compute threebody cutoff exceeds ghost atom range - "
                 "use comm_modify cutoff command");
    if (force->pair && mycutneigh < force->pair->cutforce + skin)
      if (comm->me == 0)
        error->warning(FLERR,"Compute threebody cutoff less than neighbor "
		       "cutoff - "
                       "forcing a needless neighbor list build");

    
    deldist = cutoff_user / npos_bins;

  } else deldist = force->pair->cutforce / npos_bins;

  deldistinv = 1.0/deldist;


  // since this file is going to be huge, don't include coordinates in it.

  //  for (int i = 0; i < nbin; i++)
  //array[i][0] = (i+0.5) * delr;

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

  auto req = neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
  if (cutflag) req->set_cutoff(mycutneigh);
  
}

/* ---------------------------------------------------------------------- */

void ComputeThreeBody::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeThreeBody::init_norm()
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

void ComputeThreeBody::compute_array()
// ======================================================================
// Only counting triplets where all three distances are less than cutoff
// away.
//
// ======================================================================
{
  int i,j,m,ii,jj,inum,jnum,itype,jtype,ipair,jpair;
  int ij_bin,ik_bin,alpha_bin,dum_jk_bin;
  int k,kk,knum,ktype,kpair;
  int flat;
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

void ComputeThreeBody::set_array(double constant, double normfac,
				 double (*vf)(double,double,double,
					      double,double,int))
{
  double ulower,uupper,vlower,vupper;
  double vfrac;

  int alpha_bin,ij_bin,ik_bin;
  double gr;

  int flat;
  

    
  for (ij_bin = nskip; ij_bin < (npos_bins-nskip); ij_bin++) {
    
    ulower = ij_bin*deldist;
    uupper = (ij_bin+1)*deldist;
    
    for (ik_bin = nskip; ik_bin < npos_bins-ij_bin; ik_bin++) {
      
      vlower = ik_bin*deldist;
      vupper = (ik_bin+1)*deldist;	  

      
      for (alpha_bin = 0; alpha_bin < nangle_bins; alpha_bin ++ ) {

	
	vfrac = constant * vf(ulower,uupper,vlower,vupper,delalpha,
			      alpha_bin);

	flat = flat_index(ij_bin, ik_bin, alpha_bin, npos_bins,
			  nangle_bins,nskip);
	
	if (vfrac * normfac != 0.0) {
	  gr = histall[flat]/(vfrac *normfac);
	} else {
	  gr = 0.0;
	}

	array[flat][0] = (ij_bin + 0.5)*deldist;
	array[flat][1] = (ik_bin + 0.5)*deldist;
	array[flat][2] = (alpha_bin + 0.5)*delalpha;

	array[flat][3] = gr;
      }
    }
  }

  return;
}

int ComputeThreeBody::flat_index(int ubin, int vbin, int abin,
				 int Nuv, int Na, int nsk)
{

  int udum = ubin - nsk;
  return abin + (vbin-nsk + udum*(Nuv-2*nsk+1)
		 - (udum * (udum+1) ) / 2)*Na;
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
	  * (aupper-alower));
}

double compute_2Dvf(double ulower, double uupper,
		    double vlower, double vupper,
		    double delalpha,
		    int /* alpha_bin */)
{

  return ( (uupper*uupper - ulower*ulower)/2.0
	   * (vupper*vupper - vlower*vlower)/2.0*delalpha);
}



