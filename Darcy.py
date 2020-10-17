import matplotlib.pyplot as plt
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt


def solveDarcy(mesh, K, bc):
    #u = h (hydraulic head), v = q (specific flux), differential equation:
    h = pg.solve(mesh, a = K, bc = bc) #m
    #Return the discrete interpolated gradient u for a given scalar field h.
    fig,ax = plt.subplots(1,1,figsize=(17, 11.5))
    ax.tick_params(labelsize=15)
    pg.show(mesh, h, ax=ax, label="hydraulic head distribution h [m]", nLevs= 9)
    return(h)