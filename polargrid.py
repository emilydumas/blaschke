import math
from scipy import array,dot,column_stack

def _safeinv(x):
    if x != 0:
        return(1.0/x)
    else:
        return(0.0)

class PolarGrid(object):
    """Geometry for polar grid functions in flat vectors"""
    def __init__(self,rmax,nr,nt=None):
        self.rmax = float(rmax)
        self.nr = nr
        if nt:
            self.nt = nt
        else:
            self.nt = nr
        self.N = 1 + nt*(nr-1) 
        self.dr = self.rmax / float(self.nr-1)
        self.dt = 2.0*math.pi / float(self.nt)

        self.rvec = array( [self.rcoord(*self.togrid(k)) for k in range(self.N)] )
        self.tvec = array( [self.tcoord(*self.togrid(k)) for k in range(self.N)] )
        self.polar_coordvec = column_stack((self.rvec,self.tvec))
        self.cart_coordvec = array( [ self.from_polar(*p) for p in self.polar_coordvec ] )
        self.cplxvec = dot(self.cart_coordvec,[1.0,1.0j])

    def __str__(self):
        return "<PolarGrid(rmax=%d,nr=%d,nt=%d) at 0x%x>" % (self.rmax,self.nr,self.nt,id(self))

    def from_polar(self,r,t):
        return (r*math.cos(t), r*math.sin(t))

    def toidx(self,i,j):
        if i==0:
            return 0
        else:
            return 1 + j + (i-1)*self.nt

    def togrid(self,k):
        if k==0:
            return 0,0
        else:
            return 1+int((k-1)/self.nt), (k-1) % self.nt

    def rcoord(self,i,j):
        return float(i)*self.dr

    def tcoord(self,i,j):
        return float(j)*self.dt

    def coord(self,i,j):
        return (self.rcoord(i,j), self.tcoord(i,j))

    def diameter_rk(self,j):
        jcomp = (j + (self.nt / 2)) % self.nt
        for i in xrange(self.nr-1,0,-1):
            yield -self.rvec[self.toidx(i,jcomp)],self.toidx(i,jcomp)
        for i in xrange(self.nr):
            yield self.rvec[self.toidx(i,j)],self.toidx(i,j)

    def diameter_indices(self,j):
        jcomp = (j + (self.nt / 2)) % self.nt
        for i in xrange(self.nr-1,0,-1):
            yield self.toidx(i,jcomp)
        for i in xrange(self.nr):
            yield self.toidx(i,j)

def _moduletest():
    gr = PolarGrid(3.0,6,6)
    print gr
    for k in range(gr.N):
        print gr.togrid(k)
    print [gr.tvec[k] for k in gr.diameter_indices(0)]
    print [gr.rvec[k] for k in gr.diameter_indices(0)]
    print gr.cart_coordvec
    print gr.cplxvec

if __name__=='__main__':
    _moduletest()

