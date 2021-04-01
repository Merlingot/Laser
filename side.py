from numpy import linspace, exp
import matplotlib as mpl
import matplotlib.pyplot as plt
from const import *
import numpy as np

# COURBE I-V ----------------------
Id = 1e-3; Vf = 2.2; I0 = Id/(exp(Vf*e/kb/temperature)-1)
Voltage = linspace(2.2, 2.45);
def I(V,T):
    return I0*(exp(e*V/(kb*T))-1)
Idiode = I(Voltage, temperature)
plt.figure()
color = plt.cm.coolwarm(linspace(0, 1,1))
mpl.rcParams['axes.prop_cycle'] = plt.cycler('color', color)
plt.plot()
plt.ylabel('Courant [A]')
plt.xlabel('Voltage [V]')
plt.plot(Voltage, Idiode)
plt.savefig('fig/figIV.eps', format='eps')
plt.show()
# --------------------------------


# COURBE LI -----------------------------------------------------------
# vecT = np.array([15, 20, 25, 30, 35]) + 273.15
vecT = np.array([25]) + 273.15
vecB=np.zeros(len(vecT)); vec_ntr=np.zeros(len(vecT))
vecIth=np.zeros(len(vecT))
vecg=np.zeros([len(delta_n),len(vecT)])
for i in range(len(vecT)) :
    ti=vecT[i]
    m.updateT(ti)
    gainMaxi,  *rest = calcul_gain(delta_n, l, m, ti)
    popti, pcovi = curve_fit(gammaP, delta_n[gainMaxi>0], gainMaxi[gainMaxi>0])
    vecB[i], vec_ntr[i] =  popti[0], popti[1]/popti[0] #bon
    vecIth[i]=(l.alpha_r/(vecB[i]*l.Gamma) + vec_ntr[i] )*e*l.vol/l.tau_r
    vecg[:,i]=gainMaxi



## EFFICACITÉ DE PENTE
def effPente(I, a, b ):
    return a*I + b
def Pout(laser, I, Ith):
    pout = laser.eta_ext*h*laser.nu0/e*(I-Ith)
    return pout

vecPout= np.zeros([len(vecI),len(vecT)])
legLI=[None]*len(vecT)
pente= np.zeros(len(vecT)) #EFFICACITÉ DE PENTE
for j in range(len(vecT)):
    vecPout[:,j] = Pout(l, vecI, vecIth[j])
    legLI[j]='T={:.0f} C'.format(vecT[j]-273.15)
    poptj, pcovj = curve_fit(effPente, vecI, vecPout[:,j])
    pente[j] = poptj[0]

plt.figure()
color = plt.cm.coolwarm(np.linspace(0, 1,len(vecT)))
mpl.rcParams['axes.prop_cycle'] = plt.cycler('color', color)
plt.plot(vecI, vecPout*1e3)
plt.xlabel('Courant [A]')
plt.ylabel('Puissance [mW]')
plt.legend(legLI)
plt.ylim(bottom=0)
# plt.xlim(left=0.4)
# plt.savefig('fig/LI.eps', format='eps')
plt.show()
plt.close()
# Enregistrer l'efficacité de pente:
l.eff_pente=pente[0] #Meme valeur peu importe la T
