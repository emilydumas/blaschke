import numpy as np
import disclap
from wanglin import WangLinearization
from scipy import sin,cos,arange,interpolate,array,vectorize,matrix,real
import functools
import sys
from indicator import percentdone

from scipy.integrate import ode

# CONVENTIONS

# Infinitesimal invariants of the affine sphere are passed around as "a-vectors",
# which are 5-tuples consisting of
#   (u(p), du/dx(p), du/dy(p), Re(C(p)), Im(C(p)))

# The affine frame field is considered as a 3x3 matrix in which the ROWS are
#   f, f_x, f_y

def coefmat_dx(a):
    u = a[0]
    ux = a[1]
    uy = a[2]
    cr = a[3]
    ci = a[4]

    emu = np.exp(-u)

    return matrix([ [0.0, 1.0, 0.0],
                    [np.exp(u), cr*emu + 0.5*ux, -ci*emu - 0.5*uy],
                    [0.0, -ci*emu + 0.5*uy, -cr*emu + 0.5*ux] ])

def coefmat_dy(a):
    u = a[0]
    ux = a[1]
    uy = a[2]
    cr = a[3]
    ci = a[4]

    emu = np.exp(-u)

    return matrix([ [0.0, 0.0, 1.0],
                    [0.0, -ci*emu + 0.5*uy, -cr*emu + 0.5*ux],
                    [np.exp(u), -cr*emu - 0.5*ux, ci*emu + 0.5*uy] ])


def coefmat(theta,a):
    return cos(theta)*coefmat_dx(a) + sin(theta)*coefmat_dy(a)

def odeA(theta,pcdata,t,y):
    m = coefmat(theta,pcdata(t * np.exp(1j*theta)))
    ya = matrix(y)
    ya.shape = (3,3)
    return (m * ya).ravel() # flatten 3x3 matrix to 9-vector

def phase(z):
    a = np.angle(z)
    if a < 0:
        return 2.0*np.pi+a
    else:
        return a

class BlaschkeMetric(WangLinearization):
    def compute(self,c,*args,**kwargs):
        WangLinearization.compute(self,c,*args,**kwargs)
        self.ux = disclap.discdx(self.grid) * self.u
        self.uy = disclap.discdy(self.grid) * self.u
        self.yinit = np.eye(3).ravel()
        
        self.uinterp = interpolate.RectBivariateSpline(self.grid.x,self.grid.y,np.transpose(self.u.reshape((self.grid.nx,self.grid.ny))))
        self.uxinterp = interpolate.RectBivariateSpline(self.grid.x,self.grid.y,np.transpose(self.ux.reshape((self.grid.nx,self.grid.ny))))
        self.uyinterp = interpolate.RectBivariateSpline(self.grid.x,self.grid.y,np.transpose(self.uy.reshape((self.grid.nx,self.grid.ny))))

    def pcdata(self,z):
        cval = self.c(z)
        return np.array( [self.uinterp(z.real,z.imag),
                          self.uxinterp(z.real,z.imag),
                          self.uyinterp(z.real,z.imag),
                          cval.real,
                          cval.imag] )

    def rayint(self,r,theta,step=0.01,tol=0.00001):
        odef = functools.partial(odeA,theta,self.pcdata)
        solver = ode(odef)
        solver.set_integrator("vode",method="adams",with_jacobian=False,first_step=(step*r),rtol=tol,nsteps=5000)
        solver.set_initial_value(self.yinit,0.0)
        solver.integrate(r)
        return solver.y.reshape((3,3))


    def getboundary(self,c,n,r=None,theta0=0.0,step=0.01,tol=0.0001,result='affine'):
        thetas = np.linspace(theta0, theta0 + 2*np.pi, num=n, endpoint=False)
        BlaschkeMetric.compute(self,c)
        if r==None:
            r = 0.9 * self.grid.r
        sys.stderr.write('ODE: integrating %d rays\n' % n)
        prog = percentdone(n,stream=sys.stderr,prefix='ODE: rays')
        vlist = []
        for theta in thetas:
            res = self.rayint(r,theta,step=step,tol=tol)
            if result == 'affine':
                vlist.append( res[0,1:] / res[0,0] )
            elif result == 'homog':
                vlist.append( res[0,:] )
            elif result == 'frame':
                vlist.append( res )
            prog.step()
        return vlist

    # def getimages(self,zlist,step=0.01,tol=0.0001):
    #     # BAD PRACTICE: This function assumes a previous call to BlaschkeMetric.compute!
    #     vlist = []
    #     for z in zlist:
    #         r,t = abs(z), phase(z)
    #         if r > 0.65*self.rmax:
    #             sys.stderr.write('ODE: Warning: requested interior point is near the boundary; accuracy will be poor.\n')
    #         j = int(t / self.dt) # Nearest theta that has been computed
    #         sys.stderr.write('ODE: Using r=%f, t=%f, j=%d (of %d) for z=%f+j%f\n' % (r,t,j,self.nt,z.real,z.imag))
    #         res = self.rayint(j,r,step=step,tol=tol)
    #         if result=="point":
    #             x,y,z = res
    #             vlist.append( array( (y/x,z/x) ) )
    #         else:
    #             vlist.append(res)
    #     return vlist

def _moduletest():
    import squaregrid
    def c(z):
        return z*z
    gr = squaregrid.SquareGrid(5.0,100)
    B = BlaschkeMetric(gr)
    vlist = B.getboundary(c,50)
    for v in vlist:
        print v[0],v[1]

if __name__=='__main__':
    _moduletest()
