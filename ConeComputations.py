
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
pi=3.14159265359
def ConeAngle(R,r,h):
    return 2*pi*sqrt((R+r)**2-h**2)/R+r
def Slop(X1,X2):
    return 1/((X2[0]-X1[0])/(X2[1]-X1[1]))
def bVal(X1,X2):
    return X2[1]-Slop(X1,X2)*X2[0]
    
def VoidRadius(X1,X2):
    return sqrt( X1[0]**2+(X1[1]-bVal(X1,X2))**2)
def couroneRadius(X1,X2):
    return sqrt((X1[0]-X2[0])**2+(X1[1]-X2[1])**2)
def alpha(X1,X2):
    #angle in radian
    VoidRad=VoidRadius(X1,X2)
    r=couroneRadius(X1,X2)
    #h=sqrt((VoidRad+r)**2-X2[1]**2)
    return  2*pi *X2[0]/(VoidRad+r)
