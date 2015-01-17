## blaschke

For any complex polynomial P, this package allows one to compute the
Blaschke metric and affine frame field of the complete hyperbolic
affine sphere in R^3 with Pick differential P(z)dz^3.

The projectivization of this affine sphere is the interior of a convex
polygon in RP^2 (see [our paper][affine-paper]).  The main program in this
package, *ngon.py*, uses the affine frame field to compute the
vertices of this polygon.

[affine-paper]:  http://arxiv.org/abs/1407.8149

The Blaschke metric is computed by applying the forward Euler method
to the discretization of Wang's equation on a rectangular grid.

The frame field is then computed by integrating the affine structure
equations using an ODE solver from `scipy.integrate`.

### Requirements

* Python 2.7+
* Scipy
* Numpy


### Usage examples

The main program is `ngon.py` which takes a polynomial as input and
writes a list of vertices as output.

* Compute the vertices of the polygon corresponding to z^2 dz^3:

        python ngon.py --vertices 1 0 1
  
  By default the polynomial is specified by listing its coefficients
  starting with the constant term, hence z^2 corresponds to (0, 0, 1).

  The vertices are written to stdout, one per line.  The result should
  be projectively equivalent to the regular pentagon.

* As above, but also compute the interior points of the polygon that
  correspond to the points 1+j and 1-j in the complex plane (where "j"
  is the imaginary unit, in python notation):

        python ngon.py --vertices 1 0 1 --images 1+1j 1-1j

  Note that "--images" can be followed by a list of any number of
  points.
  
* Compute the vertices and developed images of the Pick zeros for a
  monic polynomial with zeros at 1+j, 1-j, and -j:
    
        python ngon.py --vertices --root-images --roots -- 1+1j 1-1j -1j  
  
  The "--roots" option signals that the complex numbers given are the
  roots rather than the coefficients.
  
  (Numbers starting with "-" confuse the option parser, so a double
  hyphen is used to signal that there are no further options on the
  command line.)


### Status

The code is still under active development and is of alpha quality.  Comments, questions, suggestions, and bug reports are welcomed.

### Applications: Moduli spaces and the fence conjecture

This package was created to explore the mapping between polynomial and
polygon moduli spaces that is the subject of
[this preprint][affine-paper].

Early numerical experiments with this package led to the formulation of the
["fence conjecture"][fence-conjecture] for the Pick zeros of a convex
polygon in the projective plane.

[fence-conjecture]: http://dumas.io/fence-conjecture/

### Acknowledgement

This material is based upon work supported by the National Science
Foundation. Any opinions, findings, and conclusions or recommendations
expressed in this material are those of the author and do not
necessarily reflect the views of the National Science Foundation.

### Authors

David Dumas <david@dumas.io> and Michael Wolf <mwolf@math.rice.edu>
