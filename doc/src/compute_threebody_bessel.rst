.. index:: compute threebody/bessel

compute threebody/bessel command
================================


Syntax
""""""

.. parsed-literal::

   compute ID group-ID threebody/bessel Nposbins Pi keyword/value ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* Nposbin = number of bins in radial position r
* Pi = constant to be specified in window function.
* zero or more keyword/value pairs may be appended
* keyword = *cutoff* 

  .. parsed-literal::

       *cutoff* value = Rcut
         Rcut = cutoff distance for threebody computation (distance units)
	 

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all threebody/bessel 100 3.0 cutoff 16

Description
"""""""""""

Define a computation that calculates the correlation function
:math:`K^{(3)}(r)`, where :math:`r` is the distance
between particles 1 and 2. This compute is still a three-body
calculation because of it's definition,

.. math::

   K^{(3)}(r) \equiv
   = \frac{V^3}{N^3}\sum_{i=1}^N\sum_{j\neq i}^N\sum_{k\neq j; k\neq i}^N
   \bigg<\delta^2(\mathbf{r}-\mathbf{r}_{ij})K_1(\sqrt{Pi/2}r_{ik})
   \cos\theta_{jk}\bigg>.

where :math:`K_1(x)` is the modified bessel function of the second kind
at integer order 1. This compute is only available in 2D system.

   
The correlation functions are calculated in histogram form by binning
pairwise distances and angles. :math:`r` takes *Nposbins*
equi-spaced values from 0.0 to the maximum
force cutoff defined by the :doc:`pair_style <pair_style>`
command or the cutoff distance *Rcut* specified via the *cutoff* keyword.
The bins are of uniform size in radial or angular distance.  Thus a
single radial bin
encompasses a thin shell of distances in 3d and a thin ring of distances in
2d, and a single angular bin encompasses 
:math:`\cos(\alpha)-\cos(\alpha+d\alpha)` in 3d and
:math:`d\alpha` in 2D, in accordance with standard spherical-polar and
polar coordinates.

See the :doc:`compute rdf <compute_rdf>` command for details on how
binning occurs, and how to use a different cutoff than that defined in
the :doc:`pair_style <pair_style>` command via :doc:`rerun <rerun>`
and :doc:`comm_modify cutoff <comm_modify>`.


The simplest way to output the results of the compute rdf calculation
to a file is to use the :doc:`fix ave/time <fix_ave_time>` command, for
example:

.. code-block:: LAMMPS

   compute myTHREEBODY all threebody/bessel 100 3.0 
   fix 1 all ave/time 1 1 1 c_myTHREEBODY[*] file tmp.3bod mode vector

Output info
"""""""""""

Compute style *threebody/bessel* outputs a global array with the number
of rows = *Nposbins*, and two columns (one for the positions :math:`r`,
and one for :math:`K^{(3)}(r)`.


See the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The array values calculated by this compute are all "intensive".

The :math:`K^{(3)}` column of array values are normalized
numbers >= 0.0. 

Restrictions
""""""""""""

This compute may only be used for 2D systems.

See the restrictions section of :doc:`compute rdf <compute_rdf>`
for a discussion on force cutoffs and assumptions of homogeneity.

.. note::

   compute rdf can handle dynamic groups and systems where atoms
   are added or removed, but this causes that certain normalization
   parameters need to be re-computed in every step and include collective
   communication operations. This will reduce performance and limit
   parallel efficiency and scaling. For systems, where only the type
   of atoms changes (e.g. when using :doc:`fix atom/swap <fix_atom_swap>`),
   you need to explicitly request the dynamic normalization updates
   via :doc:`compute_modify dynamic yes <compute_modify>`

Related commands
""""""""""""""""

:doc:`compute rdf <compute_rdf>`, :doc:`fix ave/time <fix_ave_time>`,
:doc:`compute_modify <compute_modify>`, :doc:`compute adf <compute_adf>`

Default
"""""""

The keyword defaults are *cutoff* = 0.0 (use the pairwise force cutoff)
and *skip* = 0.
