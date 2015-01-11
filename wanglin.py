import numpy as np
from scipy import sparse,apply_along_axis,vectorize
from scipy.sparse.linalg import spsolve
import disclap
import math
import sys

def bfunc(a):
    cs = a[0]
    u = a[1]
    return (2.0*math.exp(u) - 4.0*math.exp(-2.0*u)*cs)

def cfunc(a):
    cs = a[0]
    u = a[1]
    return (-2.0*math.exp(u) - 8.0*math.exp(-2.0*u)*cs)

def smoothed_step(t):
    if t>1.0:
        return 1.0
    elif t<0.0:
        return 0.0
    else:
        return 0.5*(1.0 - math.cos(math.pi*t))

# In what follows, cs means "norm C squared", where C is the cubic differential

def uinitfunc(req, z, cs):
    '''Smoothed flat metric for initial guess; is exactly the flat metric when |z| > req
    cs = |C|^2'''
    rho = smoothed_step(1.0 - (abs(z) / req))
    return (1.0/3.0)*math.log(2.0*rho*math.exp(-cs*cs)+2.0*cs)

class WangLinearization(object):
    """Linearization of Wang equation on a polar grid"""
    def __init__(self,c,grid,thresh=0.0000001,maxiter=30,req=None):
        self.grid = grid
        self.thresh = thresh
        self.maxiter = maxiter
        self.lapmat = disclap.disclap(grid)
        if req == None:
            self.req = 0.9*self.grid.r
        else:
            self.req = req

        self.c = c
        self.cvec = vectorize(c)(self.grid.zv)
        self.csqvec = abs(self.cvec*self.cvec)
        self.u = self._iterate()

    def op_bvec(self,uvec):
        return apply_along_axis(bfunc,0,(self.csqvec,uvec)) - (self.lapmat * uvec)

    def op_mulvec(self,uvec):
        return apply_along_axis(cfunc,0,(self.csqvec,uvec))

    def todiag(self,v):
        return sparse.spdiags(v, [0], self.grid.N, self.grid.N )

    def linmat(self,uvec):
        return self.lapmat + self.todiag(self.op_mulvec(uvec))

    def step(self,u):
        udot = spsolve(self.linmat(u),self.op_bvec(u))
        return udot

    def _iterate(self):
        self.uzero = vectorize(uinitfunc)(self.req,self.grid.zv,self.csqvec)
        u = np.copy(self.uzero)
        n = 0
        delta = max(abs(self.op_bvec(u)))
        sys.stderr.write('PDE: GOAL=%f\n' % self.thresh)
        while (n < self.maxiter) and (delta > self.thresh):
            udot = self.step(u)
            u = u + 0.9*udot
            delta = max(abs(self.op_bvec(u)))
            udotnorm = max(abs(udot))
            n = n + 1
            sys.stderr.write('PDE: error=%f  delta=%f\n' % (delta,udotnorm))
        return u

def _moduletest():
    import squaregrid
    def c(z):
        return 2.0*z*z*z - 1j*z + 0.2
    gr = squaregrid.SquareGrid(3.0,31)
    W = WangLinearization(c,gr)
    j = int(gr.ny / 2)
    for i in range(gr.nx):
        k = gr.toidx(i,j)
        z = gr.zv[k]
        print 'u(%g%+gi) = \t%f  (diff from uzero is %f)' % (z.real,z.imag,W.u[k],W.u[k]-W.uzero[k])

if __name__=='__main__':
    _moduletest()

