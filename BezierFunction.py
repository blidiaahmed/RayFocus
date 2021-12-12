#define the bezier functions
import scipy.special

def Bezier(i,n,x):
    return scipy.special.binom(n, i)*(x**i)*(1-x)**(n-i)
