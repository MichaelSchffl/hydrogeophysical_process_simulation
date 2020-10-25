import matplotlib.pyplot as plt
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt


def solveDarcy(mesh, K, bc):
    #u = h (hydraulic head), v = q (specific flux), differential equation:
    h = pg.solve(mesh, a = K, bc = bc) #m
    #Return the discrete interpolated gradient u for a given scalar field h.
    fig,ax = plt.subplots()
    ax.tick_params(labelsize=15)
    pg.show(mesh, h, ax=ax, label="hydraulic head distribution h [m]", nLevs= 9)
    fig.savefig('images/hydr_grad2.png', bbox_inches='tight')
    return(h)