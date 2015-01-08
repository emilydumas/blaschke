## blaschke

Given a complex polynomial, compute the Blaschke metric and affine frame field of the associated affine sphere in R^3.

The Blaschke metric is computed by applying the forward Euler method to the discretization of Wang's equation (i.e. the vortex equation) on a polar grid.

The frame field is then computed by integrating the affine structure equations using an ODE solver from *scipy.integrate*.

### Example

Compute the vertices of the polygon corresponding to z^2 dz^3.  The result should be projectively equivalent to the regular pentagon.  Note that the polynomial is specified by its tuple of coefficients, z^2 = [0, 0, 1].

    python ngon.py --vertices 1 0 1

The vertices are written to stdout, one per line.

### Known issues

The vertices are computed by a clustering algorithm, which is an attempt to limit the influence of some numerical instability in the ODE integration step.

This should be fixed by integrating the osculation map instead.

### Authors

David Dumas <david@dumas.io> and Michael Wolf <mwolf@math.rice.edu>
