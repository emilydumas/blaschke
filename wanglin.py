from polarlap import PolarLaplacian
from scipy import sparse,apply_along_axis,vectorize
from scipy.sparse.linalg import spsolve
from itertools import partial
import math
import cmath
import sys

def bfunc(a):
    cs = a[0]
    u = a[1]
    return (2.0*math.exp(u) - 4.0*math.exp(-2.0*u)*cs)

def cfunc(a):
    cs = a[0]
    u = a[1]
    return (-2.0*math.exp(u) - 8.0*math.exp(-2.0*u)*cs)

# In what follows, cs means "norm C squared", where C is the cubic differential

def uinitfunc(rfrac, cs):
    '''Smoothed flat metric for initial guess; is exactly the flat metric on the boundary.
    rfrac = (r/rmax),
    cs = |C|^2'''
    rho = 1.0 - (rfrac)**2
    return (1.0/3.0)*math.log(2.0*rho*math.exp(-cs*cs)+2.0*cs)

class WangLinearization(PolarLaplacian):
    """Linearization of Wang equation on a polar grid"""
    def __init__(self,rmax,nr,nt=None,thresh=0.0000001,maxiter=30):
        PolarLaplacian.__init__(self,rmax,nr,nt)
        self.thresh = thresh
        self.maxiter = maxiter

    def compute(self,c,*args,**kwargs):
        self.cvec = vectorize(c)(self.cplxvec)
        self.csqvec = abs(self.cvec*self.cvec)
        self.u = self.iterate(*args,**kwargs)

    def op_bvec(self,uvec):
        return apply_along_axis(bfunc,0,(self.csqvec,uvec)) - (self.plapmat * uvec)

    def op_mulvec(self,uvec):
        return apply_along_axis(cfunc,0,(self.csqvec,uvec))

    def linmat(self,uvec):
        return self.plapmat + self.todiag(self.op_mulvec(uvec))

    def step(self,u):
        udot = spsolve(self.linmat(u),self.op_bvec(u))
        return udot

    def iterate(self,uzero=None):
        if uzero != None:
            u = uzero
        else:
            u = vectorize(uinitfunc)((1/self.rmax)*self.rvec,self.csqvec)
        n = 0
        delta = max(abs(self.op_bvec(u)))
        sys.stderr.write('PDE: GOAL=%f\n' % self.thresh)
        while (n < self.maxiter) and (delta > self.thresh):
            udot = self.step(u)
            u = u + 0.7*udot
            delta = max(abs(self.op_bvec(u)))
            n = n + 1
            sys.stderr.write('PDE: delta=%f\n' % delta)
        return u

def _moduletest():
    def c(z):
        return 2.0*z*z*z
    gr = WangLinearization(8.0,25,50)
    gr.compute(c)
    for r,k in gr.diameter_rk(0):
        print r, gr.u[k]

if __name__=='__main__':
    _moduletest()

