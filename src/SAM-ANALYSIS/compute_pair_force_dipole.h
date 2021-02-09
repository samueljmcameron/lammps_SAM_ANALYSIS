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

ComputeStyle(pair/force/dipole,ComputePairForceDipole)

#else

#ifndef LMP_COMPUTE_PAIR_FORCE_DIPOLE_H
#define LMP_COMPUTE_PAIR_FORCE_DIPOLE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePairForceDipole : public Compute {
 public:
  ComputePairForceDipole(class LAMMPS *, int, char **);
  void init();
  void init_list(int, class NeighList *);
  double compute_scalar();

 private:
  class NeighList *list; // half neighbor list
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute pair/force/dipole requires atom attribute mu

UNDOCUMENTED

E: No pair style is defined for compute pair/force/dipole

UNDOCUMENTED

E: Pair style does not support compute pair/force/dipole

UNDOCUMENTED

*/
