#!/usr/bin/python
import sys

def birapport(p0,p1,p2,p3,p4):
    '''Cross ratio of lines [p0p2], [p1p2], [p3p2], [p4p2]'''
    x0,y0 = p2  # base point

    x1,y1 = p0  # for the rest, use the cyclic order
    x2,y2 = p1
    x3,y3 = p3
    x4,y4 = p4

    m1 = (y1-y0)/(x1-x0)
    m2 = (y2-y0)/(x2-x0)
    m3 = (y3-y0)/(x3-x0)
    m4 = (y4-y0)/(x4-x0)

    return -(m3-m2)*(m1-m4) / ((m2-m1)*(m4-m3))

def parser(f):
    for line in f:
        fields = line.split()
        if len(fields)==2:
            try:
                x,y = float(fields[0]), float(fields[1])
                yield (x,y)
            except Exception:
                pass

def polygon_5chain_invts(verts):
    '''Compute cross ratios of lines from each vertex to its four closes neighbors'''
    nverts = len(verts)
    verts = verts + verts
    return [ birapport(*verts[i:i+5]) for i in range(nverts) ]

if __name__=='__main__':
    verts = [p for p in parser(sys.stdin)]
    xlist = polygon_5chain_invts(verts)
    for x in xlist:
        print x,

