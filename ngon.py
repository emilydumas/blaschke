#!/usr/bin/python
'''Compute the polygon corresponding to an affine sphere with a given polynomial Pick differential'''

_VERSION = '0.1.1'

# DEFAULT SOLVER PARAMETERS
PDE_THRESH=5e-7
PDE_NMESH=50
PDE_RMAX=6.0
ODE_NRAY=50
ODE_RMAX=3.0
ODE_THRESH=0.00001

# OPTIONS CURRENTLY NOT EXPOSED ON COMMAND LINE
_ODE_RSTEP=0.005

import sys
import argparse

parser = argparse.ArgumentParser(description='Computes (degree+3)-gon corresponding to polynomial.\n')
parser.add_argument('coefs',type=complex,nargs='+',help='coefficients, constant term first, complex() format')
parser.add_argument('--version', action='version', version='%(prog)s '+_VERSION)

group_basic = parser.add_argument_group('basic options')
group_basic.add_argument('-q','--quiet',action='store_true',help='do not put comment lines (date/time, command line, etc.) in output')
group_basic.add_argument('-t','--theta',type=float,default=0.0,help='multply by exp(i theta)')
group_basic.add_argument('-b','--boundary',action='store_true',help='output full boundary')
group_basic.add_argument('-v','--vertices',action='store_true',help='output the vertices')
group_basic.add_argument('-r','--roots',action='store_true',help='instead of coefficients, roots are given')
group_basic.add_argument('-z','--zero',action='store_true',help='roots sum to zero (only meaningful if -r | --roots given)')

group_basic.add_argument('--images',type=complex,nargs='+',help='compute images of additional points under developing map (each point given in complex() format)')
group_basic.add_argument('--root-images',action='store_true',help='compute the images of the roots of the polynomial (requires --roots)')

group_advanced = parser.add_argument_group('advanced options')
group_advanced.add_argument('--nmesh',type=int,default=PDE_NMESH,help='Number of mesh points on each axis')
group_advanced.add_argument('--nray',type=int,default=ODE_NRAY,help='Number of boundary points to compute')
group_advanced.add_argument('--rmax',type=float,default=PDE_RMAX,help='Max radius for Blaschke metric computation')
group_advanced.add_argument('--rint',type=float,default=ODE_RMAX,help='Max radius for ODE integration')
group_advanced.add_argument('--pde-epsilon',type=float,default=PDE_THRESH,help='Tolerance for PDE solver')
group_advanced.add_argument('--ode-epsilon',type=float,default=ODE_THRESH,help='Tolerance for ODE solver')

args = parser.parse_args()

if not args.boundary and not args.vertices:
    args.vertices = True

if args.root_images and not args.roots:
    parser.error(message='Option --roots required if --root-images is specified')

theta = args.theta

roots = None
if args.roots:
    # convert roots to polynomial
    roots = args.coefs
    from operator import mul
    from itertools import combinations
    if args.zero:
        kmin = 1
        roots.append(-1.0*sum(roots))
    else:
        kmin = 0
    degree = len(roots)
    coefs = []
    for k in range(degree,kmin,-1):
        c = 0
        for subset in combinations(roots,k):
            c += reduce(mul, subset)
        if k%2 == 1:
            c = -c
        coefs.append(c)
    if args.zero:
        coefs.append(0.0)
    coefs.append(1.0)
else:
    coefs = args.coefs
    degree = len(coefs)-1

sys.stderr.write('degree=%d coefs=%s\n' % (degree,coefs))
nvert = degree+3

import numpy as np
from indicator import percentdone
from squaregrid import SquareGrid
from blaschkemet import BlaschkeMetric
from functools import partial

def poly(theta,coefs,z):
    result = 0j
    for c in reversed(coefs):
        result = result*z + c
    return np.exp(1j * theta)*result

# Compute Blaschke metric, frame field, approximate boundary

gr = SquareGrid(args.rmax,args.nmesh)
c = partial(poly,theta,coefs)
bl = BlaschkeMetric(c,gr)

# Output

if not args.quiet:
    from email.Utils import formatdate
    print '#',formatdate(localtime=True)
    print '#',' '.join(sys.argv)
    print '#'

if args.boundary:
    if not args.quiet:
        print '# BOUNDARY'
    ptlist = bl.integrate_rays(n=args.nray,r=args.rint,step=_ODE_RSTEP,tol=args.ode_epsilon,return_type='affine')
    for p in ptlist:
        print '%f %f' % tuple(p)
    print ''

if args.vertices:
    if not args.quiet:
        print '# VERTICES'
    star_angle = np.angle(coefs[-1])/float(degree+3)
    framelist = bl.integrate_rays(n=nvert,theta0=star_angle,r=args.rint,step=_ODE_RSTEP,tol=args.ode_epsilon,return_type='frame')
    for F in framelist:
        w,v = np.linalg.eig(np.transpose(F))
        i = np.argmax(abs(w))
	# TODO: Investigate why x,y sometimes have tiny imaginary parts
        x,y = v[1,i] / v[0,i], v[2,i] / v[0,i]
	print x.real, y.real
    print ''

if args.root_images:
    if not args.quiet:
        print '# ROOT IMAGES'
    ptlist = bl.integrate_points(roots,step=_ODE_RSTEP,tol=args.ode_epsilon,return_type='affine')
    for p in ptlist:
        print '%f %f' % tuple(p)
    print ''

if args.images:
    if not args.quiet:
        print '# POINT IMAGES'
    ptlist = bl.integrate_points(args.images,step=_ODE_RSTEP,tol=args.ode_epsilon,return_type='affine')
    for p in ptlist:
        print '%f %f' % tuple(p)
    print ''

