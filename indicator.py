import sys

class percentdone(object):
    def __init__(self,N,stream=sys.stdout,prefix=''):
        self.N = N
        self.f = stream
        self.i = 0
        self.tic = 0
        self.prefix=prefix

    def step(self):
        if (self.i == 0):
            self.output()
        self.set(self.i+1)

    def set(self,i):
        self.i = i
        newtic = int(10.0*float(self.i) / float(self.N))
        if newtic != self.tic:
            self.tic = newtic
            self.output()

    def output(self):
        self.f.write('%s %3.1f%% done\n' % (self.prefix, 10.0*self.tic))
        
