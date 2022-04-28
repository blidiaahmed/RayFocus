# this is the other version 
# using a simple model based on basic functions

import numpy as np
from BezierFunction import *
#start by computing the expression
from sympy import *
nb=15
deg=4
width=1.25
X=[i for i in range(-nb,nb+1)]
X=[width*i/(nb) for  i in X]
BBc0=[i*i*10 for i in range(deg+1)]
Y=[sum(BBc0[j]*Bezier(j,deg,X[i]) for j in range(len(BBc0))) for i in range(len(X))]
def indexed_variable(Sym,i):
    return Symbol(''.join([Sym,'_',str(i)]))

SX=[indexed_variable('x',i) for i in range(len(X))]
BBc=[indexed_variable('c',i) for i in range(deg+1)]
SY=[sum(BBc[j]*Bezier(j,deg,X[i]) for j in range(len(BBc))) for i in range(len(X))]
from sympy.physics.vector import *
R = ReferenceFrame('R')
BY=[X[i]*R.x+Y[i]*R.y for i in range(len(X))]
B=[X[i]*R.x+SY[i]*R.y for i in range(len(X))]
BSX=[SX[i]*R.x+SY[i]*R.y for i in range(len(X))]

N=[-(B[1]-B[0])^R.z]+[-(B[i+1]-B[i-1])^R.z for i in range(1,len(X)-1)]+[-(B[len(X)-1]-B[len(X)-2])^R.z]
NY=[-(BY[1]-BY[0])^R.z]+[-(BY[i+1]-BY[i-1])^R.z for i in range(1,len(X)-1)]+[-(BY[len(X)-1]-BY[len(X)-2])^R.z]
NY=[NY[i]/sqrt(NY[i]&NY[i]) for i in range(len(NY))]
NSX=[-(BSX[1]-BSX[0])^R.z]+[-(BSX[i+1]-BSX[i-1])^R.z for i in range(1,len(X)-1)]+[-(B[len(X)-1]-B[len(X)-2])^R.z]

S=R.y#+1*R.x
S/=sqrt(S&S)
concentrationPoint=0.75
A=concentrationPoint*R.y#-0*R.x


#
Exp=[N[i]^(sqrt((A-B[i])&(A-B[i]))*S+A-B[i]) for i in range(len(X))]
#


ExpX=[NSX[i]^(sqrt((A-BSX[i])&(A-BSX[i]))*S+A-BSX[i]) for i in range(len(X))]


MinimisationExp=100*sum([dot(i,i) for i in Exp])+1*sum(i**2 for i in BBc )+5*indexed_variable('c',0)**2
#+20*(indexed_variable('c',1)-concentrationPoint)**2

#compute minimisation
from scipy.optimize import minimize
from sympy import symbols, lambdify
import matplotlib.pyplot as plt
f=lambdify(BBc, MinimisationExp)
def rosen(x):
    return f(*x)
p=[lambdify(BBc, Exp[i]&R.z) for i in range(len(Exp))]
def objective(x):
    return [i(*x) for i in p]
res = minimize(rosen, np.array(BBc0), method='nelder-mead',
               options={'xatol': 1e-8, 'disp': True,"maxiter": 500})
# drawing the normals
BBcSol=res.x



Ysol=[sum(BBcSol[j]*Bezier(j,deg,X[i]) for j in range(len(BBcSol))) for i in range(len(X))]

Bsol=[X[i]*R.x+Ysol[i]*R.y for i in range(len(X))]
Nsol=[-(Bsol[1]-Bsol[0])^R.z]+[-(Bsol[i+1]-Bsol[i-1])^R.z for i in range(1,len(X)-1)]+[-(Bsol[len(X)-1]-Bsol[len(X)-2])^R.z]


