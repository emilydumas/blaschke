#!/usr/bin/python
'''Compute the polygon corresponding to an affine sphere with a given polynomial Pick differential'''

# DEFAULT SOLVER PARAMETERS

PDE_THRESH=5e-7
PDE_NR=50
PDE_NTHETA=150
PDE_RMAX=6.0
ODE_RAYSTRIDE=2
ODE_RSTEP=0.005
ODE_RMAX=3.0
ODE_TOL=0.00001


import sys
import argparse

parser = argparse.ArgumentParser(description='Computes (degree+3)-gon corresponding to polynomial.\n')
parser.add_argument('-t','--theta',type=float,default=0.0,help='multply by exp(i theta)')
parser.add_argument('-b','--boundary',action='store_true',help='output full boundary')
parser.add_argument('-v','--vertices',action='store_true',help='output the vertices')
parser.add_argument('-r','--roots',action='store_true',help='instead of coefficients, roots are given')
parser.add_argument('-z','--zero',action='store_true',help='roots sum to zero (only meaningful if -r | --roots given)')
parser.add_argument('coefs',type=complex,nargs='+',help='coefficients, constant term first, complex() format')
parser.add_argument('--develop',type=complex,nargs='+',help='other points to develop, complex() format')
parser.add_argument('--nr',type=int,default=PDE_NR,help='Number of radius steps')
parser.add_argument('--ntheta',type=int,default=PDE_NTHETA,help='Number of theta steps')
parser.add_argument('--rmax',type=float,default=PDE_RMAX,help='Maximum integration radius')

args = parser.parse_args()

if not args.boundary and not args.vertices:
    args.vertices = True

theta = args.theta

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

PDE_NR = args.nr
PDE_NTHETA = args.ntheta
PDE_RMAX = args.rmax

sys.stderr.write('degree=%d coefs=%s\n' % (degree,coefs))
nvert = degree+3

from cmath import exp
from indicator import percentdone
from blaschkemet import BlaschkeMetric
from clusterng import cluster
from functools import partial

def poly(theta,coefs,z):
    result = 0j
    for c in reversed(coefs):
        result = result*z + c
    return exp(1j * theta)*result

bl = BlaschkeMetric(PDE_RMAX,PDE_NR,PDE_NTHETA,thresh=PDE_THRESH)
c = partial(poly,theta,coefs)

ptlist = bl.getboundary(c,r=0.65*PDE_RMAX,stride=ODE_RAYSTRIDE,step=ODE_RSTEP,tol=ODE_TOL)

if args.boundary:
    for p in ptlist:
        print '%f %f' % tuple(p)

if args.vertices:
    sys.stderr.write('Clustering.\n')
    verts = cluster(ptlist,nvert)
    if len(verts) != nvert:
        raise ValueError('Wrong number of clusters:  Expected %d, found %d' % (nvert, len(verts))) 
    if args.boundary:
        print ''
    for v in verts:
        print '%f %f' % tuple(v)

if args.roots:
    vlist = bl.getimages(roots,step=ODE_RSTEP,tol=ODE_TOL)
    print ''
    for v in vlist:
        print '%f %f' % tuple(v)

if args.develop:
    vlist = bl.getimages(args.develop,step=ODE_RSTEP,tol=ODE_TOL)
    print ''
    for v in vlist:
        print '%f %f' % tuple(v)
