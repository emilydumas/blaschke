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
    DEFAULT_RET_TYPE='affine'
    def __init__(self,*args,**kwargs):
        WangLinearization.__init__(self,*args,**kwargs)

        # Prepare for ODE integration
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

    def prepare_output(self,F,return_type='affine'):
        if return_type == 'affine':
            return F[0,1:] / F[0,0]
        elif return_type == 'homog':
            return F[0,:]
        elif return_type == 'frame':
            return F

    def integrate_polar(self,r=None,theta=0.0,step=0.01,tol=0.00001,return_type='affine'):
        # If r not given, use circle 90% of radius of grid incircle
        if r==None:
            r = 0.9 * self.grid.r
        odef = functools.partial(odeA,theta,self.pcdata)
        solver = ode(odef)
        solver.set_integrator("vode",method="adams",with_jacobian=False,first_step=(step*r),rtol=tol,nsteps=50000)
        solver.set_initial_value(self.yinit,0.0)
        solver.integrate(r)
        return self.prepare_output(solver.y.reshape((3,3)), return_type)

    def integrate_point(self,z,step=0.01,tol=0.00001,return_type=DEFAULT_RET_TYPE):
        return self.integrate_polar(r=abs(z),theta=np.angle(z),step=step,tol=tol,return_type=return_type)

    def integrate_rays(self,n,r=None,theta0=0.0,step=0.01,tol=0.0001,return_type=DEFAULT_RET_TYPE):
        thetalist = np.linspace(theta0, theta0 + 2*np.pi, num=n, endpoint=False)
        return [ self.integrate_polar(r=r,theta=theta,step=step,tol=tol,return_type=return_type) for theta in thetalist ]

    def integrate_points(self,zlist,step=0.01,tol=0.00001,return_type=DEFAULT_RET_TYPE):
        return [ self.integrate_point(z,step=step,tol=tol,return_type=return_type) for z in zlist ]

def _moduletest():
    import squaregrid
    def c(z):
        return z*z
    gr = squaregrid.SquareGrid(5.0,100)
    B = BlaschkeMetric(c,gr)
    vlist = B.integrate_rays(50)
    for v in vlist:
        print v[0],v[1]

if __name__=='__main__':
    _moduletest()
