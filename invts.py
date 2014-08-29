#!/usr/bin/python
import sys

def birapport(p0,p1,p2,p3,p4):
    x0,y0 = p0
    x1,y1 = p1
    x2,y2 = p2
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

verts = [p for p in parser(sys.stdin)]
nvert = len(verts)

verts = verts + verts

xlist = [ birapport(*verts[i:i+5]) for i in range(nvert) ]
for x in xlist:
    print x,

