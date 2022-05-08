
#In this file we organie the computations in classes

from helperFile import  *
import numpy as np
from BezierFunction import *
from sympy import *
from sympy.physics.vector import *

class parameterClass:
    """these are the parameters"""
    def __init__(self,pd=3,SN=20,SW=1.25,XCP=0.75,YCP=0):
        self.polynomialDegree=pd
        self.SamplesNumber=SN
        self.SectionWidth=SW
        self.XconcentrationPoint=XCP
        self.YconcentrationPoint=YCP
        

        

class Samples:
    def __init__(self,Param=parameterClass()):
        self.parameters=Param
        self.XSamples=[self.parameters.SectionWidth*i/self.parameters.SamplesNumber 
                       for  i in range(self.parameters.SamplesNumber+1)]
class Bezier:   
    def InitialiseCoefficients(self,BBc0):
        self.splineCoefficients=BBc0
        
    def __init__(self,sam=Samples()):
        self.sample=sam
        Pdeg=self.sample.parameters.polynomialDegree
        XSam=self.sample.XSamples
        self.splineCoefficientsSymboles=[indexed_variable('c',i) for i in range(Pdeg+1)]
        self.InitialiseCoefficients([i*i*10 for i in 
                                range(self.sample.parameters.polynomialDegree+1)])
        



