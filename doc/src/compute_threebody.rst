.. index:: compute threebody
.. index:: compute threebody/angle/cos

compute threebody command
=========================


compute threebody/angle/cos command
===================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID style_name Nposbins Nanglebins keyword/value ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* style_name = *threebody* or *threebody/angle/cos*
* Nposbin = number of bins in radial positions u and v
* Nanglebins = number of bins in angular distribution
* zero or more keyword/value pairs may be appended
* keyword = *cutoff* or *skip*

  .. parsed-literal::

       *cutoff* value = Rcut
         Rcut = cutoff distance for threebody computation (distance units)
       *skip* value = nskip
         nskip = number of bins to skip (starting from bin 0).
	 

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all threebody 100 100
   compute 1 all threebody 100 100 cutoff 5.0
   compute 1 all threebody/angle/cos 100 100 skip 20 cutoff 16.0

Description
"""""""""""

Define a computation that calculates the three-body correlation function
:math:`g^{(3)}(u,v,\theta_{uv})`, where :math:`u` is the distance
between particles 1 and 2, :math:`v` is the distance between
particles 1 and three, and :math:`\alpha_{uv}` is the angle between
the vectors which have distances :math:`u` and :math:`v`. The domain of
:math:`u` and :math:`v` are restricted (both) restricted to the same
*cutoff* value. The domain of :math:`\theta_{uv}\in[0,\pi]`. In full
generality, the three-body correlation function :math:`\tilde{g}_g^{(3)}`
depends on nine coordinates (in 3D), and is defined (for an infinite
system) as

.. math::

   \tilde{g}_g^{(3)}( \mathbf{x},\mathbf{y},\mathbf{z} )
   = \frac{V^3}{N^3}\sum_{i=1}^N\sum_{j\neq i}^N\sum_{k\neq j; k\neq i}^N
   \bigg<\delta^d(\mathbf{x}-\mathbf{r}_i)\delta^d(\mathbf{y}-\mathbf{r}_j)
   \delta^d(\mathbf{z}-\mathbf{r}_k)\bigg>

but in this compute we assume
translational and rotational invariance, i.e. for an
arbitrary vector :math:`\mathbf{a}`, we assume

.. math::

   \tilde{g}_g^{(3)}(\mathbf{x}+\mathbf{a},\mathbf{y}+\mathbf{a},\mathbf{z}+\mathbf{a})
   = \tilde{g}_g^{(3)}(\mathbf{x},\mathbf{y},\mathbf{z}).

Since this vector can be any value, choose it to be :math:`-\mathbf{x}`, so
that

.. math::

   \tilde{g}^{(3)}(\mathbf{0},\mathbf{y}-\mathbf{x},\mathbf{z}-\mathbf{x})
   \equiv g^{(3)}_g(\mathbf{u},\mathbf{v})

which therefore allows one to write

.. math::
   
   g^{(3)}_g(\mathbf{u},\mathbf{v})
   = \frac{V^3}{N^3}\sum_{i=1}^N\sum_{j\neq i}^N\sum_{k\neq j; k\neq i}^N
   \bigg<\delta^d(\mathbf{r}_i)\delta^d(\mathbf{u}-\mathbf{r}_j)
   \delta^d(\mathbf{v}-\mathbf{r}_k)\bigg>\\
   = \frac{V^2}{N^3}\sum_{i=1}^N\sum_{j\neq i}^N\sum_{k\neq j; k\neq i}^N
   \bigg<\delta^d(\mathbf{u}-\mathbf{r}_{ij})
   \delta^d(\mathbf{v}-\mathbf{r}_{ik})\bigg>.

Finally, it is assumed that the system is rotationally invariant, so that

.. math::

   g^{(3)}_g(\mathbf{u},\mathbf{v}) \equiv g^{(3)}(u,v,\alpha_{uv})

where :math:`\cos\alpha_{uv} = \mathbf{u}\cdot\mathbf{v}/(uv)`.

The final simplification depends on the dimensions of the system. For a
two dimensional system, it is fairly straightforward as rotations commute,
so :math:`\alpha_{uv} = \theta_v - \theta_u`, the difference between polar
angles of the two vectors (as measured in some arbitrary coordinate system).
So, one know that for a rotationally invariant system an arbitrary
angle :math:`\beta` can be added onto the polar angles, leading to

.. math::

   g^{(3)}_{2d}(u,v,\alpha_{uv})
   = \frac{V^2}{N^3}\sum_{i=1}^N\sum_{j\neq i}^N\sum_{k\neq j; k\neq i}^N
   \bigg<\frac{\delta(u-r_{ij})\delta(\theta_u+\beta-\theta_{ij})
   \delta(v-r_{ik})\delta(\theta_v+\beta-\theta_{ik})}{r_{ij}r_{ik}}\bigg>\\
   = \frac{V^2}{N^3}\sum_{i=1}^N\sum_{j\neq i}^N\sum_{k\neq j; k\neq i}^N
   \bigg<\frac{\delta(u-r_{ij})
   \delta(v-r_{ik})\delta(\alpha_{uv}-\alpha_{jk})}{2\pi uv}\bigg>.


.. note::
   
   In 2D, :math:`\alpha_{uv}>0` configurations
   are equivalent to :math:`\alpha_{uv}<0` configurations, and so only the
   domain :math:`\alpha_{uv}\in[0,\pi]` of the three body correlation
   function is saved (the :math:`\alpha_{uv}<0` are mapped onto the
   equivalent :math:`\alpha_{uv}>0`, effectively double-counting these
   configurations - this is remedied by dividing the counts by a factor of
   two). 

For the three dimensional case the rotational invariance manifests
in a somewhat more complicated way, as :math:`\alpha_{uv}`
is not simply the difference in polar or
azimuthal angle in an arbitrary reference frame. We must then
write :math:`\mathbf{v}` in terms of perpendicular and parallel components of
:math:`\mathbf{u}`, so that
      
.. math::

   \mathbf{v}=v(\cos\alpha_{uv}\hat{\mathbf{u}}_{\parallel}
   +\sin\alpha_{uv}\cos\gamma_{uv}\hat{\mathbf{u}}_{\perp,1}
   +\sin\alpha_{uv}\sin\gamma_{uv}\hat{\mathbf{u}}_{\perp,2})

which means that one can write the three body correlation function as

.. math::
   
   g^{(3)}_{3d}(u,v,\alpha_{uv})
   = \frac{V^2}{N^3}\sum_{i=1}^N\sum_{j\neq i}^N\sum_{k\neq j; k\neq i}^N
   \bigg<\frac{\delta(u-r_{ij})\delta(\theta_u-\theta_{ij})
   \delta(\phi_u-\phi_{ij})
   \delta(v-r_{ik})\delta(\alpha_{uv}-\alpha_{ik})
   \delta(\gamma_{uv}-\gamma_{ik})}
   {r_{ij}^2\sin(\theta_{ij})r_{ik}^2\sin(\alpha_{ik})}\bigg>\\
   = \frac{V^2}{N^3}\sum_{i=1}^N\sum_{j\neq i}^N\sum_{k\neq j; k\neq i}^N
   \bigg<\frac{\delta(u-r_{ij})\delta(v-r_{ik})
   \delta(\alpha_{uv}-\alpha_{ik})}{8\pi^2u^2v^2\sin(\alpha_{uv})}\bigg>.

   
   
The correlation functions are calculated in histogram form by binning
pairwise distances and angles. :math:`u` and :math:`v` take *Nposbins*
equi-spaced values from 0.0 to the maximum
force cutoff defined by the :doc:`pair_style <pair_style>`
command or the cutoff distance *Rcut* specified via the *cutoff* keyword.
:math:`\alpha_{uv}` takes *Nanglebins* equi-spaced values from 0 to
:math:`\pi`.
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

.. note::

   Due to the amount of data storage required for this (three variable)
   function and the lack of demand, there is no option to specify
   different correlation functions for different particle pairs as in
   :doc:`compute rdf <compute_rdf>`. Such an extension would be not to
   difficult to do if it was desired.

.. note::

   Due to the finite cutoff of the radial variables :math:`u` and :math:`v`,
   all angles :math:`\alpha_{uv}\in[0,\pi]` are only counted when
   :math:`|\mathbf{u}-\mathbf{v}|\leq r_c`, where :math:`r_c` is the
   cutoff (specified by the :doc:`pair_style <pair_style>` or the
   *cutoff* keyword). This is equivalent to requiring
   :math:`u+v\leq r_c` since both :math:`u>0` and :math:`v>0`.

The simplest way to output the results of the compute rdf calculation
to a file is to use the :doc:`fix ave/time <fix_ave_time>` command, for
example:

.. code-block:: LAMMPS

   compute myTHREEBODY all threebody 100 100
   fix 1 all ave/time 1 1 1 c_myTHREEBODY[*] file tmp.3bod mode vector

Output info
"""""""""""

Compute style *threebody* outputs a global array with the number of rows =
:math:`(N_p-2n_s)(N_p-2n_s+1)N_a/2`, where :math:`N_p` = *Nposbins*,
:math:`n_s` = *nskip* (which is zero if not specified) and
:math:`N_a` = *Nanglebins*. The reason for this strange number is
that the domain of validity for this function is when :math:`u+v\leq r_c`
(as noted above), and so to avoid wasting space only values of the
three body correlation function within this range of validity are output.
The number of columns = 1, holding the values of the three body
correlation function. For example, if *skip* has value :math:`n_s`,
*cutoff* has value :math:`r_c`, then
the position spacings are :math:`dr= r_c/(N_p-1)`, and the output
will be a single, flattened (column) array with index given by

.. math::

   i = j_{\alpha}+ \bigg(j_v-nskip + (j_u-nskip)*(Nposbins-2nskip+1)
   - \frac{(j_u-nskip)*(j_u-nskip+1)}{2}\bigg)*Nanglebins

where the unflattened array would have indices
:math:`(j_u,j_v,j_{\alpha})` for the :math:`(u,v,\alpha)` variables.
This output may seem complex,
but it allows for more compact storage of this data, while also
avoiding the domain where the calculation is invalid.

Compute style *threebody/angle/cos* outputs a global array where
:math:`g^{(3)}(u,v,\alpha_{uv})` is integrated over the angle coordinate,
i.e.
      
.. math::
   2\int_0^{\pi}g^{(3)}_{2d}(u,v,\alpha_{uv})f(\cos\alpha_{uv})
   d\alpha_{uv},\\
   \int_0^{\pi} g^{(3)}_{3d}(u,v,\alpha_{uv})f(\cos\alpha_{uv})
   \sin\alpha_{uv} d\alpha_{uv}

where the angular factor :math:`f(\cos\alpha_{uv})=1` for the first
output column and :math:`f(\cos\alpha_{uv})=\cos\alpha_{uv}` for
the second output column.

The angle integral is computed via summing over the *Nanglebins*.
The resulting arrays are therefore only dependent on the two
position coordinates :math:`u` and :math:`v`, and so the output
reflects this (but is otherwise the same as described for the
*threebody* compute style).

See the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The array values calculated by this compute are all "intensive".

The :math:`g^{(3)}` column of array values are normalized
numbers >= 0.0. 

Restrictions
""""""""""""

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
