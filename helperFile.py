import numpy as np 
from sympy import *
from sympy.physics.vector import *


def indexed_variable(Sym,i):
    return Symbol(''.join([Sym,'_',str(i)]))

def drawCircles(angle,radius):
    nbSample=200
    angleSamples = np.linspace( 0 , angle , nbSample )
    x = radius * np.cos( angleSamples ) 
    y = radius * np.sin( angleSamples ) 
    return x,y