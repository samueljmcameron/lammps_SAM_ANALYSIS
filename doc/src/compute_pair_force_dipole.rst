.. index:: compute pair/force/dipole

compute pair/force/dipole command
=================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID pair/force/dipole

* ID, group-ID are documented in :doc:`compute <compute>` command
* pair = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all pair/force/dipole

Description
"""""""""""

Compute the scalar

.. math::


   = \frac{1}{2d(d-1)V}\sum_{i=1}^{N'}\sum_{j\neq i}^{N'}
   \mathbf{F}_{ij}\cdot\mathbf{e}_{ij}


which is relevant for computing active brownian particle pressure.
:math:`\mathbf{e}_i` is the dipole moment of particle :math:`i`.
When periodic boundary conditions are
used, N' necessarily includes periodic image (ghost) atoms outside the
central box, and the position and force vectors of ghost atoms are thus
included in the summation.  When periodic boundary conditions are not
used, N' = N = the number of atoms in the system.


Output info
"""""""""""

This compute calculates an "intensive" global scalar as specified above.

The scalar value will be in pressure :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute pair <compute_pair>`, :doc:`compute pair/local <compute_pair_local>`

Default
"""""""

none
