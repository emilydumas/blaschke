from wanglin import WangLinearization
from scipy import sin,cos,arange,interpolate,array,vectorize,matrix,real
import functools
import sys
import cmath
from indicator import percentdone

try:
    from scipy.integrate import complex_ode
except ImportError:
    from cxode import complex_ode

def pmat(a):
    u = a[0]
    ux = a[1]
    uy = a[2]
    cr = a[3]
    ci = a[4]
    return matrix([ [0.0,1.0,0.0], 
                    [0,0.5*(ux - 1j*uy),(cr+1j*ci)*cmath.exp(-u)],
                    [0.5*cmath.exp(u),0.0,0.0] ])

def qmat(a):
    u = a[0]
    ux = a[1]
    uy = a[2]
    cr = a[3]
    ci = a[4]
    return matrix([ [0.0,0.0,1.0], 
                    [0.5*cmath.exp(u),0.0,0.0],
                    [0.0,(cr-1j*ci)*cmath.exp(-u),0.5*(ux + 1j*uy)] ])

def coefmat(theta,a):
    return cmath.exp(1j*theta)*pmat(a) + cmath.exp(-1j*theta)*qmat(a)

def odeA(theta,ucf,t,y):
    m = coefmat(theta,ucf(t))
    ya = matrix(y)
    ya.shape = (3,3)
    return (m * ya).ravel()

def phase(z):
    a = cmath.phase(z)
    if a < 0:
        return 2.0*cmath.pi+a
    else:
        return a

class BlaschkeMetric(WangLinearization):
    def compute(self,c,*args,**kwargs):
        WangLinearization.compute(self,c,*args,**kwargs)
        self._ut = self.dtmat * self.u
        self._ur = self.drmat * self.u
        self.ux = cos(self.tvec)*self._ur - sin(self.tvec)*self.irvec*self._ut
        self.uy = sin(self.tvec)*self._ur + cos(self.tvec)*self.irvec*self._ut
        self.yinit = array([1,0,0,0,0.5,0.5j,0,0.5,-0.5j])
        
    def diamvec(self,j):
        """vector of values (r,u,ux,uy,c)"""
        return array([ (r,self.u[k],self.ux[k],self.uy[k],self.cvec[k].real, self.cvec[k].imag) for r,k in self.diameter_rk(j) ])

    def theta(self,j):
        return self.tvec[j]

    def diamfunc(self,j):
        v = self.diamvec(j)
        return interpolate.interp1d(v[:,0],v[:,1:],"quadratic",axis=0)

    def rayint(self,j,r,step=0.01,tol=0.00001):
        ucf = self.diamfunc(j)  # ucf = "u and c function; gives (u,ux,uy,creal,cimag)"
        odef = functools.partial(odeA,self.theta(j),ucf)
        solver = complex_ode(odef)
        solver.set_integrator("vode",method="adams",with_jacobian=False,first_step=(step*r),rtol=tol,nsteps=5000)
        solver.set_initial_value(self.yinit,0.0)
        solver.integrate(r)
        return self.coord_normalize(solver.y)

    def coord_normalize(self,y):
        return (y[0].real, y[1].real, y[2].real)

    def getboundary(self,c,r=None,stride=2,step=0.01,tol=0.0001):
        BlaschkeMetric.compute(self,c)
        if r==None:
            r = 0.65 * self.rmax
        nrays = int(self.nt / stride)
        sys.stderr.write('ODE: integrating %d rays (stride=%d)\n' % (nrays,stride))
        prog = percentdone(nrays,stream=sys.stderr,prefix='ODE: rays')
        vlist = []
        for j in xrange(0,self.nt,stride):
            res = self.rayint(j,r,step=step,tol=tol)
            x,y,z = res
            vlist.append( array( (y/x,z/x) ) )
            prog.step()
        return vlist

    def getimages(self,zlist,step=0.01,tol=0.0001):
        # BAD PRACTICE: This function assumes a previous call to BlaschkeMetric.compute!
        vlist = []
        for z in zlist:
            r,t = abs(z), phase(z)
            if r > 0.65*self.rmax:
                sys.stderr.write('ODE: Warning: requested interior point is near the boundary; accuracy will be poor.\n')
            j = int(t / self.dt) # Nearest theta that has been computed
            sys.stderr.write('ODE: Using r=%f, t=%f, j=%d (of %d) for z=%f+j%f\n' % (r,t,j,self.nt,z.real,z.imag))
            res = self.rayint(j,r,step=step,tol=tol)
            if result=="point":
                x,y,z = res
                vlist.append( array( (y/x,z/x) ) )
            else:
                vlist.append(res)
        return vlist

def _moduletest():
    def c1(z):
        return z*z+0j
    bl = BlaschkeMetric(6.0,60,180,thresh=0.00001)
    bl.compute(c1)
    for j in range(bl.nt):
        x,y,z=bl.rayint(j,4.0)
        print y/x, z/x

if __name__=='__main__':
    _moduletest()
