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

ComputeStyle(threebody/angle/cos,ComputeThreeBodyAngleCos)

#else

#ifndef LMP_COMPUTE_THREEBODY_ANGLE_COS_H
#define LMP_COMPUTE_THREEBODY_ANGLE_COS_H

#include "compute_threebody.h"

namespace LAMMPS_NS {

class ComputeThreeBodyAngleCos : public ComputeThreeBody {
 public:
  ComputeThreeBodyAngleCos(class LAMMPS *, int, char **);
  ~ComputeThreeBodyAngleCos();

  void compute_array();

 private:
  
  void set_array(double, double,
		 double (*)(double,double,double,double,double,int));

  int flat_index_no_theta(int , int , int ,int);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute rdf requires a pair style be defined or cutoff specified

UNDOCUMENTED

E: Compure rdf cutoff exceeds ghost atom range - use comm_modify cutoff command

UNDOCUMENTED

W: Compute rdf cutoff less than neighbor cutoff - forcing a needless neighbor list build

UNDOCUMENTED

U: Compute rdf requires a pair style be defined

Self-explanatory.

*/
