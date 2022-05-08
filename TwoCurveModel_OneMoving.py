# this is the other version 
# using a simple model based on basic functions

import numpy as np
from BezierFunction import *
#start by computing the expression
from sympy import *
from helperFile	import *
# عدد العينات النقطية
nb=10
# درجة كثيرات الحدود
deg=3
# عرض الحزام الضوئي
width=20.5

# العينات النقطية
X=[i+5 for i in range(nb+1)]
X=[width*i/(nb) for  i in X]

# نقاط التحكم الابتدائية
BBc0=[-(i*i/500+10) +20 for i in range(deg+1)]
BBd0=[-(i*i/500) +20 for i in range(deg+1)]

#نقاط التحكم المتحركة

BBd=[indexed_variable('d',i) for i in range(deg+1)]

#فواصل العينات نقطية من الدالة بدلالة نقاط التحكم

SYd=[sum(BBd[j]*Bezier(j,deg,X[i]) for j in range(len(BBd))) for i in range(len(X))]

SYc0=[sum(BBc0[j]*Bezier(j,deg,X[i]) for j in range(len(BBc0))) for i in range(len(X))]
SYd0=[sum(BBd0[j]*Bezier(j,deg,X[i]) for j in range(len(BBd))) for i in range(len(X))]

from sympy.physics.vector import *

R = ReferenceFrame('R')

#عينات نقطية من الدالة بدلالة نقاط التحكم

Bd=[X[i]*R.x+SYd[i]*R.y for i in range(len(X))]

Bc0=[X[i]*R.x+SYc0[i]*R.y for i in range(len(X))]
Bd0=[X[i]*R.x+SYd0[i]*R.y for i in range(len(X))]

#المرايا
mirrors0=[(Bd0[i+1]-Bc0[i]) for i in range(len(X)-1)]
mirrors=[(Bd[i+1]-Bc0[i]) for i in range(len(X)-1)]

#الشعاع الناظم المرافق للنقاط
#N=[-(Bc[1]-Bc[0])^R.z]+[-(Bc[i+1]-Bc[i-1])^R.z for i in range(1,len(X)-1)]+[-(Bc[len(X)-1]-Bc[len(X)-2])^R.z]
N=[-mirrors[i]^R.z for i in range(len(X)-1)]
#N0=[-(Bd0[i]-Bc0[i+1])^R.z for i in range(0,nb-1)]+[-(Bd0[i+1]-Bc0[i])^R.z for i in range(nb-1,len(X)-2)]
N0=[-mirrors0[i]^R.z for i in range(len(X)-1)]

N0=[i/sqrt(i&i) for i in N0]


#اتجاه الأشعة الواردة S
S=R.y#+1*R.x
S/=sqrt(S&S)
concentrationPoint=0#min(SYc0)-10

#َ A نقطة التركيز
A=concentrationPoint*R.y#-*R.x

#
Exp=[N[i]^(sqrt((A-Bc0[i])&(A-Bc0[i]))*S+A-Bc0[i]) for i in range(1,len(X)-1)]
#

ControlPoints=BBd.copy()
#ControlPoints.extend(BBd)

MinimisationExp=sum([dot(i,i) for i in Exp])#+9000000*sum(i**2 for i in ControlPoints )#+500*indexed_variable('c',0)**2
#+20*(indexed_variable('c',1)-concentrationPoint)**2

#compute minimisation
from scipy.optimize import minimize
from sympy import symbols, lambdify
import matplotlib.pyplot as plt



f=lambdify(ControlPoints, MinimisationExp)
def rosen(x):
    return f(*x)

InitialControlPoints=BBd0.copy()

res = minimize(rosen, np.array(InitialControlPoints), method='nelder-mead',
               options={'xatol': 1e-8, 'disp': True,"maxiter": 100000})

# تمثيل الحلول هندسيا
#نقاط التحكم 
BBdSol=res.x[range(len(BBd))]

#تراتيب المنحنى
Ysolc=[sum(BBc0[j]*Bezier(j,deg,X[i]) for j in range(len(BBc0))) for i in range(len(X))]
Ysold=[sum(BBdSol[j]*Bezier(j,deg,X[i]) for j in range(len(BBdSol))) for i in range(len(X))]
# نقاط من المنحنى
BSolc=[X[i]*R.x+Ysolc[i]*R.y for i in range(len(X))]
BSold=[X[i]*R.x+Ysold[i]*R.y for i in range(len(X))]
#المرايا
mirrorsSol=[(BSold[i+1]-BSolc[i]) for i in range(len(X)-1)]
#النواظم
Nsol=[-mirrorsSol[i]^R.z for i in range(len(X)-1)]

