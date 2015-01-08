from polargrid import PolarGrid
import math
from scipy import sparse,vectorize
import sys

def _safeinv(x):
    if x != 0:
        return(1.0/x)
    else:
        return(0.0)

safeinv = vectorize(_safeinv)

class PolarLaplacian(PolarGrid):
    """Differentiation matrices for polar grid"""
    def __init__(self,rmax,nr,nt=None):
        PolarGrid.__init__(self,rmax,nr,nt)
        sys.stderr.write("PDE:Generating matrices...")
        self.irvec = safeinv(self.rvec)
        self.drmat = self.make_drmat()
        self.d2rmat = self.make_d2rmat()
        self.dtmat = self.make_dtmat()
        self.d2tmat = self.make_d2tmat()
        self.plapmat = self.make_plapmat()
        sys.stderr.write("done.\n")

    def __str__(self):
        return "<PolarLaplacian(rmax=%d,nr=%d,nt=%d) at 0x%x>" % (self.rmax,self.nr,self.nt,id(self))

    def todiag(self,v):
        return sparse.spdiags(v, [0], self.N, self.N )

    def make_drmat(self):
        m = sparse.dok_matrix((self.N,self.N))
        epsilon = 1.0 / self.dr
        
        # df/dr always zero at origin; iterate over other points
        for i in range(1,self.nr-1):
            for j in range(0,self.nt):
                m[self.toidx(i,j),self.toidx(i-1,j)] += -0.5*epsilon
                m[self.toidx(i,j),self.toidx(i+1,j)] += 0.5*epsilon

        return m.tocsr()

    def make_d2rmat(self):
        m = sparse.dok_matrix((self.N,self.N))
        esq = 1.0 / (self.dr*self.dr)
        
        # Regular points -- standard finite difference 2nd deriv
        for i in range(1,self.nr-1):
            for j in range(0,self.nt):
                m[self.toidx(i,j),self.toidx(i-1,j)] += esq
                m[self.toidx(i,j),self.toidx(i+1,j)] += esq
                m[self.toidx(i,j),self.toidx(i,j)] += -2.0*esq
        # Special handling for the origin
        # (average over circle of radius dr minus value at (0,0))
        for j in range(self.nt):
            m[self.toidx(0,0),self.toidx(1,j)] += 2.0*esq / float(self.nt)
        m[self.toidx(0,0),self.toidx(0,0)] += -2.0*esq

        return m.tocsr()

    def make_dtmat(self):
        m = sparse.dok_matrix((self.N,self.N))
        epsilon = 1.0/self.dt
        
        # df/dtheta always zero at the origin; iterate over other points
        for i in range(1,self.nr-1):
            for j in range(0,self.nt):
                m[self.toidx(i,j),self.toidx(i,(j+1)%self.nt)] += 0.5*epsilon
                m[self.toidx(i,j),self.toidx(i,(j-1)%self.nt)] += -0.5*epsilon
 
        return m.tocsr()


    def make_d2tmat(self):
        m = sparse.dok_matrix((self.N,self.N))
        esq = 1.0/(self.dt*self.dt)
        
        # d^2f/dtheta^2 always zero at the origin; iterate over other points
        for i in range(1,self.nr-1):
            for j in range(0,self.nt):
                m[self.toidx(i,j),self.toidx(i,(j+1)%self.nt)] += esq
                m[self.toidx(i,j),self.toidx(i,(j-1)%self.nt)] += esq
                m[self.toidx(i,j),self.toidx(i,j)] += -2.0*esq
 
        return m.tocsr()

    def make_plapmat(self):
        invr = self.todiag(self.irvec)
        invr2 = invr * invr
        return (self.d2rmat + invr*self.drmat + invr2*self.d2tmat)

    def makevec(self,f):
        """Make vector from complex-valued function"""
        return vectorize(f)(self.cplxvec)

def _moduletest():
    gr = PolarLaplacian(5.0,5,5)
    print gr
    print gr.plapmat.todense()

if __name__=='__main__':
    _moduletest()

