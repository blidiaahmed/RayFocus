# this is the other version 
# using a simple model based on basic functions

import numpy as np
from BezierFunction import *
#start by computing the expression
from sympy import *
from helperFile	import *
# عدد العينات النقطية
nb=20
# درجة كثيرات الحدود
deg=3
# عرض الحزام الضوئي
width=1.5

# العينات النقطية5
X=[i for i in range(-nb,nb+1)]
X=[width*i/(nb) for  i in X]

# نقاط التحكم الابتدائية
BBc0=[i*i*10 for i in range(deg+1)]

#نقاط التحكم المتحركة
BBc=[indexed_variable('c',i) for i in range(deg+1)]
SY=[sum(BBc[j]*Bezier(j,deg,X[i]) for j in range(len(BBc))) for i in range(len(X))]
from sympy.physics.vector import *

R = ReferenceFrame('R')
B=[X[i]*R.x+SY[i]*R.y for i in range(len(X))]

N=[-(B[1]-B[0])^R.z]+[-(B[i+1]-B[i-1])^R.z for i in range(1,len(X)-1)]+[-(B[len(X)-1]-B[len(X)-2])^R.z]

#اتجاه الأشعة الواردة S
S=R.y#+1*R.x
S/=sqrt(S&S)
concentrationPoint=0.65

#َ A نقطة التركيز
A=concentrationPoint*R.y#-0*R.x


#
Exp=[N[i]^(sqrt((A-B[i])&(A-B[i]))*S+A-B[i]) for i in range(len(X))]
#


MinimisationExp=100*sum([dot(i,i) for i in Exp])+1*sum(i**2 for i in BBc )+500*indexed_variable('c',0)**2
#+20*(indexed_variable('c',1)-concentrationPoint)**2

#compute minimisation
from scipy.optimize import minimize
from sympy import symbols, lambdify
import matplotlib.pyplot as plt
f=lambdify(BBc, MinimisationExp)
def rosen(x):
    return f(*x)

res = minimize(rosen, np.array(BBc0), method='nelder-mead',
               options={'xatol': 1e-8, 'disp': True,"maxiter": 10000})
# drawing the normals
BBcSol=res.x



Ysol=[sum(BBcSol[j]*Bezier(j,deg,X[i]) for j in range(len(BBcSol))) for i in range(len(X))]

Bsol=[X[i]*R.x+Ysol[i]*R.y for i in range(len(X))]
Nsol=[-(Bsol[1]-Bsol[0])^R.z]+[-(Bsol[i+1]-Bsol[i-1])^R.z for i in range(1,len(X)-1)]+[-(Bsol[len(X)-1]-Bsol[len(X)-2])^R.z]

