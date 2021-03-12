/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(threebody/bessel,ComputeThreeBodyBessel)

#else

#ifndef LMP_COMPUTE_THREEBODY_BESSEL_H
#define LMP_COMPUTE_THREEBODY_BESSEL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeThreeBodyBessel : public Compute {
 public:
  ComputeThreeBodyBessel(class LAMMPS *, int, char **);
  ~ComputeThreeBodyBessel();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();

 protected:
  int npos_bins;                 // # of radial bins
  int cutflag;                   // user cutoff flag

  double sqrtPi2fac;             // sqrt(Pi/2.) where Pi = D_r sigma^2/D_t is input by user

  double deldist,deldistinv;     // bin width and its inverse for distance
  double cutoff_user;            // user-specified cutoff
  double mycutneigh;             // user-specified cutoff + neighbor skin
  double **hist;                 // histogram bins
  double **histall;              // summed histogram bins across all procs

  int typecount;
  int icount,jcount;
  int duplicates;

  class NeighList *list; // half neighbor list
  void init_norm();
  bigint natoms_old;

 private:

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute threebody/bessel command is incompatible with 
   newton pair being on.

Self-explanatory.

E: Cannot compute threebody/bessel for system with multiple
   atom types.

Self-explanatory.

E: Compute threebody/bessel requires a pair style be defined or cutoff specified

Self-explanatory.

E: Compure threebody/bessel cutoff exceeds ghost atom range - use comm_modify cutoff command

UNDOCUMENTED

W: Compute threbody/bessel cutoff less than neighbor cutoff - forcing a needless neighbor list build

UNDOCUMENTED


*/
