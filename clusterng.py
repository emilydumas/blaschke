"""Try to find an n-gon given a set of points clustered around its vertices."""

import sys
from math import sqrt,atan2
from scipy import array

def dist2(p,q):
    x = p[0] - q[0]
    y = p[1] - q[1]
    return sqrt(x*x+y*y)

def mean(L):
    m = None
    n = 0
    for p in L:
        n = n + 1
        if m == None:
            m = array(p)
        else:
            m = m + array(p)
    return tuple((1.0/float(n)) * m)

def circulate_left(L,k):
    return L[k:]+L[:k]

def pairdist(L):
    return [ dist2(L[i],L[(i+1)%len(L)]) for i in range(len(L)) ]

def nearlist(L,p,epsilon):
    return [ q for q in L if dist2(p,q) < epsilon ]

def farlist(L,p,epsilon):
    return [ q for q in L if dist2(p,q) > epsilon ]

def anglefrom(center,pt):
    return atan2( pt[1]-center[1], pt[0]-center[0] )

def cluster(vlist,n=5):
    """Identify clusters, assuming there are n of them evenly spaced in list"""

    centroid = mean(vlist)

    # In case a cluster includes index 0, rotate to that first distance is the max
    dlist = pairdist(vlist)
    maxindex = dlist.index(max(dlist))
    vlist = circulate_left(vlist,maxindex)
    dlist = circulate_left(dlist,maxindex)
    # Now all clusters are in the "interior" of vlist

    corners = []
    while len(corners) < n:
        mindist = min(dlist)
        inrad = pow(mindist,0.7)
        outrad = 15.0*pow(mindist,0.7)

        k = dlist.index(mindist)
        mid = mean( [vlist[k],vlist[k+1]] )
        cluster = nearlist(vlist,mid,inrad)
        corners.append(mean(cluster))

        vlist = farlist(vlist,mid,outrad)
        dlist = pairdist(vlist)

        sys.stderr.write('  cluster size=%d min=%f inrad=%f outrad=%f center=%s\n' % (len(cluster),mindist,inrad,outrad,str(corners[-1])))

    corners.sort(key=lambda p:anglefrom(centroid,p))
    return corners
