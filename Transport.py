import matplotlib.pyplot as plt # 2D plotting library (www.matplotlib.org)
import numpy as np # library for efficient numerical computing with arrays (www.numpy.org)

import pygimli as pg # *G*eophysical *I*nversion and *M*odeling *Li*brary (www.pygimli.org)
import pygimli.meshtools as mt

from pygimli.solver import solveFiniteVolume as solveFV

def solveTransport(mesh, D, v, S, times):
    #u = tracer concentration c [kg/l]
    c1 = solveFV(mesh, a = D, scheme="PS", vel = v, f = S, times = times)
    c2 = solveFV(mesh, a = D, scheme="PS", vel = v, f = 0, times = times, u0 = c1[-1])
    c = np.vstack([c1,c2])
    return c