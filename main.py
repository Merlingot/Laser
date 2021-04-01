import numpy as np
from numpy import sqrt, exp, linspace, trapz, absolute, argmin, argmax, zeros, real, log, sign, sin
from numpy import pi
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib as mpl
from scipy.integrate import Radau, solve_ivp
from scipy.signal import square

from papanoel import fdelec, fdpolylog, approx_eta, f
from classes import *
from gain import *
from const import *

# Figures:
from matplotlib import rc
rc('text', usetex=True)
vis=True
# col1='tab:red'
# col2 = 'tab:blue'
color = plt.cm.cool(np.linspace(0, 1,3))
mpl.rcParams['axes.prop_cycle'] = plt.cycler('color', color)
col1=color[2]; col2=color[1]
#------------------------
# PARAMETRES ----------------------
doping =-1.7e18*1e2**3 #[m^-3] Positif pour dopé P, négatif pour dopé N
# doping=0
pipi=200/3.623*850e-9
caca=300/3.623*850e-9
# pipi=9.6e-5
# caca=250/2/3.655*850e-9
l=Laser(caca, pipi)
m=Materiau(ratio_mc=0.07, ratio_mv=0.5, Eg0=1.519, A=5.4e-4, B=204, T=temperature, doping=doping)
vecI = np.array([50,100,150,200,250])*1e-3 #A
delta_n=injection(vecI, l)
delta_n_cm = delta_n*1e-2**3 # [cm-3]
# -------------------------------
## COURBE DE GAIN @T=25--------------
gainMax = calcul_du_gain(vecI, 150e-3, l, m, temperature)


## COURBE DE GAIN MAX @T=25--------------
def gammaP(deltaN, B, b):
    return B*deltaN - b
popt, pcov = curve_fit(gammaP, delta_n[gainMax>0]+m.n0, gainMax[gainMax>0])
# Enregistrer B et delta_n_tr @T=25C
l.B, l.delta_n_tr =  popt[0], popt[1]/popt[0]-m.n0

# plt.plot()
# plt.xlabel('Densité de porteurs [cm$^{-3}$]')
# plt.ylabel('Gain maximal $\gamma_{max}$ [cm$^{-1}$]')
# plt.plot(delta_n_cm, gainMax*1e-2, label='$\gamma_{max}$', color=col1)
# plt.plot(delta_n_cm, gammaP(delta_n+m.n0, *popt)*1e-2, '--', label='Régression linéaire', color=col2)
# plt.legend()
# plt.savefig('fig/gainMax.eps', format='eps')
# plt.show()
# plt.close()

# Enregistrer les Threshold @T=25C
l.I_th = (l.alpha_r/(l.B*l.Gamma) + l.delta_n_tr )*e*l.vol/l.tau_r
l.delta_n_th = l.alpha_r/(l.B*l.Gamma) + l.delta_n_tr
l.gamma_th=l.B*(l.delta_n_th-l.delta_n_tr)

# # COURBE LI -----------------------------------------------------------
# vecT = np.array([15, 20, 25, 30, 35]) + 273.15
# # vecT = np.array([25]) + 273.15
# vecB=np.zeros(len(vecT)); vec_ntr=np.zeros(len(vecT))
# vecIth=np.zeros(len(vecT))
# vecg=np.zeros([len(delta_n),len(vecT)])
# for i in range(len(vecT)) :
#     ti=vecT[i]
#     m.updateT(ti)
#     gainMaxi,  *rest = calcul_gain(delta_n, l, m, ti)
#     popti, pcovi = curve_fit(gammaP, delta_n[gainMaxi>0], gainMaxi[gainMaxi>0])
#     vecB[i], vec_ntr[i] =  popti[0], popti[1]/popti[0] #bon
#     vecIth[i]=(l.alpha_r/(vecB[i]*l.Gamma) + vec_ntr[i] )*e*l.vol/l.tau_r
#     vecg[:,i]=gainMaxi
#
# ## EFFICACITÉ DE PENTE
# def effPente(I, a, b ):
#     return a*I + b
# def Pout(laser, I, Ith):
#     pout = laser.eta_ext*h*laser.nu0/e*(I-Ith)
#     return pout
#
# vecPout= np.zeros([len(vecI),len(vecT)])
# legLI=[None]*len(vecT)
# pente= np.zeros(len(vecT)) #EFFICACITÉ DE PENTE
# for j in range(len(vecT)):
#     vecPout[:,j] = Pout(l, vecI, vecIth[j])
#     legLI[j]='T={:.0f} C'.format(vecT[j]-273.15)
#     poptj, pcovj = curve_fit(effPente, vecI, vecPout[:,j])
#     pente[j] = poptj[0]
#
# plt.figure()
# color = plt.cm.cool(np.linspace(0, 1,len(vecT)))
# mpl.rcParams['axes.prop_cycle'] = plt.cycler('color', color)
# plt.plot(vecI, vecPout*1e3)
# plt.xlabel('Courant [A]')
# plt.ylabel('Puissance [mW]')
# plt.legend(legLI)
# plt.ylim(bottom=0)
# plt.savefig('fig/LI.eps', format='eps')
# plt.show()
# plt.close()
# # Enregistrer l'efficacité de pente:
# l.eff_pente=pente[0] #Meme valeur peu importe la T




# # MODES QUI LASENT et SPECTRE @ T=25 -----------------------------------------

# laser=l; mat=m; T=temperature; I150=150e-3
#
# vec_I=np.array([I150])
# delta_n = injection(vec_I, laser)
# labellegend=[None]
# for i in range(len(vec_I)):
#     labellegend[i] = "{:10.4g} mA".format( vec_I[i]*1e3 ) #mA
# gainMax, gainMaxLoc, gamma0, nu  = calcul_gain(delta_n, laser, mat, T, vis=False, leg=labellegend)
#
# deltanuEV=laser.isl/(pi/(laser.alpha_r*laser.l))*h/e
# ISL_eV = h*laser.isl/e
# rang = nu[gamma0[:,0]>=laser.gamma_th]
# b = rang[-1] - rang[1]
# M = int(b/laser.isl)
# # MODES POSSIBLES @150mA
# nu0s = h*laser.nu0/e + ISL_eV*np.arange(-(M-1)/2, (M-1)/2+1, 1)
# gamma_m = np.zeros(len(nu0s))
#
# fig=plt.figure()
# plt.plot(nu*h/e, gamma0*1e-2, color=col1) #cm-1
# plt.legend(labellegend)
# plt.axhline(y=laser.gamma_th*1e-2, linestyle = '-', color='black', linewidth=0.75) # gamma_th cm-1
# for i in range(len(nu0s)) :
#     nui=nu0s[i]
#     loc = find_nearest(nu*h/e, nui)
#     ampli =  gamma0[loc,0]
#     gamma_m[i] = ampli
#     plt.plot(nu*h/e, 1e-2*ampli*Lorentienne(deltanuEV, nui, nu*h/e), color='black', linewidth=0.75)
# plt.xlabel('Énergie [eV]')
# plt.ylabel('Gain $\gamma$ [cm$^{-1}$]')
# plt.savefig('fig/modes.eps', format='eps')
# plt.show()
#
#
#
# # SPECTRE D'ÉMISSION @150mA
# somme = np.sum( gamma_m*c*laser.beta/(laser.n_GaAs*laser.tau_r) * (1/laser.tau_p + gamma_m*laser.Gamma*c/laser.n_GaAs )**(-1) )
# n = I150*laser.eta/(e*laser.vol)*( 1/laser.tau_r +  somme )**(-1)
# nph_m = laser.beta*n/laser.tau*( 1/laser.tau_p*1.1 - gamma_m*laser.Gamma*c/laser.n_GaAs )**(-1)
#
#
# P_m = nph_m*h*(nu0s*e/h)*laser.vol/laser.tau_m
# color = plt.cm.coolwarm(np.linspace(0, 1, 5))
# mpl.rcParams['axes.prop_cycle'] = plt.cycler('color', color)
#
# fig=plt.figure()
# plt.semilogy(nu*h/e, P_m[int((len(nu0s)-1)/2)]*1e3*Lorentienne(deltanuEV, nu0s[int((len(nu0s)-1)/2)], nu*h/e), color=col1)
# for i in range(len(nu0s)) :
#     nui=nu0s[i]
#     if i==(len(nu0s)-1)/2:
#         pass
#     else :
#         plt.semilogy(nu*h/e, P_m[i]*1e3*Lorentienne(deltanuEV, nui, nu*h/e), color='black', linewidth=0.75)
# plt.legend(['Mode central', 'Mode latéraux'])
# plt.xlabel('Énergie [eV]')
# plt.ylabel('Puissance [mW]')
# plt.savefig('fig/spectre.eps', format='eps')
# plt.show()


## MODULATION ------------------------------------------------

def square_wave(I0, Is, f, t):
    return I0 + Is*0.5*( 1. + square(2.*np.pi*f*t + np.pi))
def square_wave2(I0, Is, f, t):
    return I0 + Is*0.5*( 1. + square(2.*np.pi*f*t ))
y0=[0,0]
t0=0
F=np.array([1, 5, 10])*1e9/2
tboundS=np.array([20, 50, 100])

# nph0=l.eta*l.n_GaAs/(e*l.vol*l.alpha_r*c/l.Gamma)*(i-l.I_th)
# n0=((l.eta*i/(e*l.vol) + l.B*l.delta_n_tr*nph0*c/l.n_GaAs )/(1./l.tau +l.B*c*nph0/l.n_GaAs))
# y0=[absolute(doping)+injection(i, l), nph0]


i=200e-3
kk=100e-3
for j in range(len(F)):
    f = F[j]
    # print(f)
    t_bound=tboundS[j]/f
    # f=10e8/2
    # t_bound=10/f
    tt = np.linspace(0.,t_bound,num=500000)

    def ezpez(t, y):
        I = square_wave(i, i*kk, f, t)
        N, N_ph = y[0], y[1]
        return np.array([l.eta*I/(e*l.vol)-N/l.tau-l.B*(N-l.delta_n_tr)*c*N_ph/l.n_GaAs,
        l.B*(N-l.delta_n_tr)*c*N_ph*l.Gamma/l.n_GaAs-N_ph/l.tau_p+l.beta*N/l.tau_r])

    u = solve_ivp(ezpez, (t0,t_bound), y0, method=Radau, dense_output=True, t_eval=tt)
    data1 = u.y[0]*1e-2**3
    data2 = u.y[1]*h*l.nu0/l.tau_m*l.vol*1e3
    data3 = square_wave(i, i*kk, f, u.t)
    # plt.plot(u.t*1e9, data2, color=col2)

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('temps (ns)')
    ax1.set_ylabel('Puissance [mW]', color=col2)
    ax1.plot(u.t, data2, color=col2)
    ax1.tick_params(axis='y', labelcolor=col2)


    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel('Courant [mA]', color=col1)
    ax2.plot(u.t, data3*1e3, color=col1, linewidth=0.5)
    ax2.tick_params(axis='y', labelcolor=col1)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    # plt.savefig('fig/modul_{}.eps'.format(j), format='eps')
    plt.show() #mW

    t_bound=3/f
    tt = np.linspace(0.,t_bound,num=10000)
    u = solve_ivp(ezpez, (t0,t_bound), [u.y[0][-1], u.y[1][-1]], method=Radau, dense_output=True, t_eval=tt)
    data1 = u.y[0]*1e-2**3
    data2 = u.y[1]*h*l.nu0/l.tau_m*l.vol*1e3
    data3 = square_wave(i, i*kk, f, u.t)
    # if j==2 :
        # data3 = square_wave2(i, i*kk, f, u.t)
    # plt.plot(u.t*1e9, data2, color=col2)

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('temps (ns)')
    ax1.set_ylabel('Puissance [mW]', color=col2)
    ax1.plot(u.t, data2, color=col2)
    ax1.tick_params(axis='y', labelcolor=col2)


    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel('Courant [mA]', color=col1)
    ax2.plot(u.t, data3*1e3, color=col1, linewidth=0.5)
    ax2.tick_params(axis='y', labelcolor=col1)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig('fig/modul_steady_{}.eps'.format(j), format='eps')
    plt.show() #mW
