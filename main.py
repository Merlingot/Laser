import numpy as np
from numpy import sqrt, exp, linspace, trapz, absolute, argmin, argmax, zeros, real, log
from numpy import pi
from scipy.integrate import quad
import matplotlib.pyplot as plt
from papanoel import fdelec, fdpolylog, approx_eta, f
from scipy import interpolate
from scipy.optimize import curve_fit
from classes import *
import matplotlib as mpl
from gain import gain
# Figures:
from matplotlib import rc
rc('text', usetex=True)
vis=True

### PARAMETRES -------------------------------------------------------------
T = 273.15 + 25
doping=0
l=Laser()
m=Materiau(ratio_mc=0.07, ratio_mv=0.5, Eg0=1.519, A=5.4e-4, B=204, T=0, doping=doping)
### -----------------------------------------------------------------------

# COURBE I-V ----------------------
Id = 1e-3; Vf = 2.2; I0 = Id/(exp(Vf*e/kb/T)-1)
Voltage = linspace(2.2, 2.45);
def I(V,T):
    return I0*(exp(e*V/(kb*T))-1)
Idiode = I(Voltage, T)
plt.plot()
plt.ylabel('Courant [A]')
plt.xlabel('Voltage [V]')
plt.plot(Voltage, Idiode)
plt.show(block=vis)
plt.savefig('fig/figIV.eps', format='eps')
# --------------------------------

## COURBE DE GAIN ----------------------------------------------------
delta_n_cm = np.linspace(1.5, 2.5, num=10)*1e18 # [cm-3]
delta_n = delta_n_cm*1e2**3 #[m-3]
gainMax = gain(l, m, delta_n, doping, T, vis=True)
# COURBE GAIN MAX :
def gammaP(deltaN, B, b):
    # b = delta_N_T*a
    return B*deltaN - b
popt, pcov = curve_fit(gammaP, delta_n[gainMax>0], gainMax[gainMax>0])
B, deltanT =  popt[0], popt[1]/popt[0]

color = plt.cm.coolwarm(np.linspace(0, 1,2))
mpl.rcParams['axes.prop_cycle'] = plt.cycler('color', color)
plt.plot()
plt.xlabel('Densité de porteurs injectée [cm$^{-3}$]')
plt.ylabel('Gain maximal $\gamma_{max}$ [cm$^{-1}$]')
plt.plot(delta_n_cm, gainMax*1e-2, label='$\gamma_{max}$')
plt.plot(delta_n_cm, gammaP(delta_n, *popt)*1e-2, '--', label='Régression linéaire')
plt.legend()
plt.show(block=vis)
plt.savefig('fig/gainMax.eps', format='eps')

## COURBE LI -----------------------------------------------------------
vecT = np.array([-25, 0, 25, 50, 75]) + 273.15
#
# for Ti in vecT :
#     m.updateT(Ti)
#     delta_ni, gainMaxi = gain(l,m, delta_n, doping, ti,  vis=True)









#
