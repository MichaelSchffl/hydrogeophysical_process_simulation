import matplotlib.pyplot as plt
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt


def solveStreampot(mesh, sigma, bc, f):
    #u = streaming potential phi
    phi = pg.solve(mesh, a = sigma, f = f, bc = bc)
    E = pg.solver.grad(mesh, -phi)
    fig, axes = plt.subplots(figsize =(17,11.5))
    axes.tick_params(labelsize=15)
    #pg.show(mesh, E, ax = axes, hold=True)
    #mt.nodeDataToCellData(mesh1, phi)
    pg.show(mesh, phi*1000, label = r"streaming potential " r"$\varphi$ [mV]", ax=axes, cMap = 'RdBu_r', nLevs = 11)
    return phi,E