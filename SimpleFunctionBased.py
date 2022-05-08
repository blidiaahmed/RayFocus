# this is the other version 
# using a simple model based on basic functions

import numpy as np


#start by computing the expression
from sympy import *
nb=6

width=15
X=[i for i in range(-nb,nb+1)]
X=[width*i/(nb*2) for  i in X]

Y=[i*i for i in X]
def indexed_variable(Sym,i):
    return Symbol(''.join([Sym,'_',str(i)]))
SY=[indexed_variable('y',i) for i in range(len(X))]
SX=[indexed_variable('x',i) for i in range(len(X))]

from sympy.physics.vector import *
R = ReferenceFrame('R')
BY=[X[i]*R.x+Y[i]*R.y for i in range(len(X))]
B=[X[i]*R.x+SY[i]*R.y for i in range(len(X))]
BSX=[SX[i]*R.x+SY[i]*R.y for i in range(len(X))]

N=[-(B[1]-B[0])^R.z]+[-(B[i+1]-B[i-1])^R.z for i in range(1,len(X)-1)]+[-(B[len(X)-1]-B[len(X)-2])^R.z]
NY=[-(BY[1]-BY[0])^R.z]+[-(BY[i+1]-BY[i-1])^R.z for i in range(1,len(X)-1)]+[-(BY[len(X)-1]-BY[len(X)-2])^R.z]
NY=[NY[i]/sqrt(NY[i]&NY[i]) for i in range(len(NY))]
NSX=[-(BSX[1]-BSX[0])^R.z]+[-(BSX[i+1]-BSX[i-1])^R.z for i in range(1,len(X)-1)]+[-(B[len(X)-1]-B[len(X)-2])^R.z]
#never change the S coefficient below
S=R.y
A=19*R.y

#
Exp=[N[i]^(sqrt((A-B[i])&(A-B[i]))*S+A-B[i]) for i in range(len(X))]
#


ExpX=[NSX[i]^(sqrt((A-BSX[i])&(A-BSX[i]))*S+A-BSX[i]) for i in range(len(X))]

MinimisationExp=sum([dot(i,i) for i in Exp])

#compute minimisation
from scipy.optimize import minimize
from sympy import symbols, lambdify
import matplotlib.pyplot as plt
f=lambdify(SY, MinimisationExp)
def rosen(x):
    return f(*x)
p=[lambdify(SY, Exp[i]&R.z) for i in range(len(Exp))]
def objective(x):
    return [i(*x) for i in p]
res = minimize(rosen, np.array(Y), method='nelder-mead',
               options={'xatol': 1e-8, 'disp': True,"maxiter": 500000})
# drawing the normals
Ysol=res.x

Bsol=[X[i]*R.x+Ysol[i]*R.y for i in range(len(X))]
Nsol=[-(Bsol[1]-Bsol[0])^R.z]+[-(Bsol[i+1]-Bsol[i-1])^R.z for i in range(1,len(X)-1)]+[-(Bsol[len(X)-1]-Bsol[len(X)-2])^R.z]


