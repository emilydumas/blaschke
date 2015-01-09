#import numpy as np
#import math
from scipy import sparse
import sys


def disclap(g):
    '''Make a discrete laplacian matrix for a square mesh

    Args:
      g = SquareGrid object representing the mesh
    '''

    sys.stderr.write('PDE: Generating discrete laplacian matrix\n')
    idxs = 1.0 / (g.dx*g.dx)
    idys = 1.0 / (g.dy*g.dy)

    m = sparse.dok_matrix( (g.N,g.N) )

    # Interior points -- centered finite difference 2nd deriv
    for i in range(1,g.nx-1):
        for j in range(1,g.ny-1):
            # d/dx^2
            m[g.toidx(i,j),g.toidx(i-1,j)] += idxs
            m[g.toidx(i,j),g.toidx(i+1,j)] += idxs
            m[g.toidx(i,j),g.toidx(i,j)] += -2.0*idxs
            # d/dy^2
            m[g.toidx(i,j),g.toidx(i,j+1)] += idys
            m[g.toidx(i,j),g.toidx(i,j-1)] += idys
            m[g.toidx(i,j),g.toidx(i,j)] += -2.0*idys

    return m.tocsr()

def discdx(g):
    '''Make a discrete d/dx matrix for a square mesh

    Args:
      g = SquareGrid object representing the mesh
    '''

    sys.stderr.write('PDE: Generating discrete d/dx matrix\n')
    idx = 1.0 / g.dx

    m = sparse.dok_matrix( (g.N,g.N) )

    # Interior points -- centered finite difference 2nd deriv
    for i in range(1,g.nx-1):
        for j in range(0,g.ny):
            m[g.toidx(i,j),g.toidx(i-1,j)] += -0.5*idx
            m[g.toidx(i,j),g.toidx(i+1,j)] += 0.5*idx

    return m.tocsr()

def discdy(g):
    '''Make a discrete d/dx matrix for a square mesh

    Args:
      g = SquareGrid object representing the mesh
    '''

    sys.stderr.write('PDE: Generating discrete d/dy matrix\n')
    idy = 1.0 / g.dy

    m = sparse.dok_matrix( (g.N,g.N) )

    # Interior points -- centered finite difference 2nd deriv
    for i in range(0,g.nx):
        for j in range(1,g.ny-1):
            m[g.toidx(i,j),g.toidx(i,j-1)] += -0.5*idy
            m[g.toidx(i,j),g.toidx(i,j+1)] += 0.5*idy

    return m.tocsr()


def _moduletest():
    import squaregrid
    gr = squaregrid.SquareGrid(3.0,5)
    print gr
    print disclap(gr).todense()

if __name__=='__main__':
    _moduletest()

