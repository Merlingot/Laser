import numpy as np
from numpy import sqrt, exp, linspace, trapz, absolute, argmin, argmax, zeros, real, log
from numpy import pi
from scipy.integrate import quad
import matplotlib.pyplot as plt
from sympy import polylog, re


def f(E, mu, T):
    kb = 1.38064852e-23
    return 1/(exp((E-mu)/(kb*T))+1)

def fdelec(eta):
    def fd12(eps, eta):
        return eps**(1/2)/(1+exp(eps-eta))
    resultat=zeros(len(eta)); err=zeros(len(eta))
    for i in range(len(eta)):
        resultat[i], err[i]=quad(fd12, 0, np.inf, args=(eta[i]))
    return resultat, err


def fdpolylog(eta):
    resultat=zeros(len(eta)); err=zeros(len(eta))
    for i in range(len(eta)):
        resultat[i] = - re(sqrt(pi)/2*polylog(2/3, -exp(eta[i]) ))
    return resultat, 0


def u_approx(eta):
    eps = 3*sqrt(pi/2)*( eta+2.13 + (absolute(eta-2.13)**2.4 + 9.6)**(5/12) )**(-3/2)
    return 1/(exp(-eta) + eps)

def approx_eta(u):
    return log(u)/(1-u**2) + ((3*sqrt(pi)*u/4)**(2/3))/( 1 + (0.24+1.08*(3*sqrt(pi)*u/4)**(2/3))**(-2) )
