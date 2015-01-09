import numpy as np

class SquareGrid(object):
    """Geometry for rectangular grid centered at origin functions stored as flat vectors"""
    def __init__(self,r,n):
        '''Args:
        r = inradius
        n = number of mesh points in each direction
        '''
        self.n = int(n)
        self.nx = self.ny = n

        self.N = n*n
        self.r = float(r)
        self.xmin = self.ymin = -self.r
        self.xmax = self.ymax = self.r
        
        self.dx = (self.xmax - self.xmin) / float(n-1)
        self.dy = (self.xmax - self.xmin) / float(n-1)

        self.x = np.linspace(self.xmin,self.xmax,num=n,endpoint=True)
        self.y = np.linspace(self.ymin,self.ymax,num=n,endpoint=True)
        
        self.xm, self.ym = np.meshgrid(self.x,self.y)

        self.xv = self.xm.ravel()
        self.yv = self.ym.ravel()

        self.zm = self.xm + 1j*self.ym
        self.zv = self.xv + 1j*self.yv


    def __str__(self):
        return "<SquareGrid(r=%f,n=%d),id=0x%x>" % (self.r,self.n,id(self))

    def toidx(self,i,j):
        return i + j*self.n

    def fromidx(self,k):
        return k % self.n, int(k/self.n)

def _moduletest():
    gr = SquareGrid(3.0,6)
    print gr
    for k in range(gr.N):
        print gr.fromidx(k)
    print gr.zv

if __name__=='__main__':
    _moduletest()

