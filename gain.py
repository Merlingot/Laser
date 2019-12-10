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

from matplotlib import rc
rc('text', usetex=True)


def calcul_gain(delta_n, l, m, T):
    # ELECTroNNNNNNNS --------------------------------------------------------
    ## Calcul des populations a l'equilibre
    ni = sqrt(4*(2*pi*kb*T/h**2)**3*(m.mc*m.mv)**(3/2)*exp(-m.Eg0/(kb*T)));
    if m.doping > 0:
        m.p0 = m.doping
        m.n0 = ni**2/m.p0
    elif doping < 0:
        m.n0 = m.doping;
        m.p0 = ni**2/m.n0;
    else:
        m.n0 = ni
        m.p0 = ni
    # Ajout de la population de pompe
    n = m.n0+delta_n; p = m.p0+delta_n

    # Pour plus tard
    E = np.linspace(-1, 2, num=100)*e #nu (nm) (Plage d'Énergie pour tests)
    Efc=zeros([len(n)]); Efv=zeros([len(n)]); # Quasi niveaux de Fermi
    fc=zeros([len(E),len(n)]); fv=zeros([len(E),len(n)]); # fv et fc
    labellegend=[None]*len(n)

    #---------------------- methode integrale
    # Trouver le quasi niveau de fermi (ELECTRONS) pour tout n:
    ## Resoudre integrale de fermi pour obtenir FD(eta)
    #Range different pour les electrons et les trous !
    Efc_test=linspace(0.9, 1.2, num=1000)*m.Ec; Efv_test=linspace(0.1, -0.3, num=1000)*m.Ec;
    etaE = (Efc_test - m.Ec)/(kb*T); etaP =(m.Ev-Efv_test)/(kb*T)
    fd_etaE, errE= fdelec(etaE) #intégration numérique
    # fd_etaE, errE= fdpolylog(etaE)
    fd_etaP, errP= fdelec(etaP) #intégration numérique
    # fd_etaP, errP= fdpolylog(etaP)

    # Obtenir relations pour n et quasi niveau
    gc=(2*m.mc)**(3/2)/(2*pi**2*hbar**3); Nc=(kb*T)**(3/2)*gc; nE=Nc*fd_etaE;
    gv=(2*m.mv)**(3/2)/(2*pi**2*hbar**3); Nv=(kb*T)**(3/2)*gv; nP=Nv*fd_etaP

    efE = m.Ec + kb*T*etaE; efP = m.Ev - kb*T*etaP
    tckE=interpolate.splrep(nE, efE) #quasi niveau (Electrons) en fct de n
    tckP=interpolate.splrep(nP, efP) #quasi niveau (Trous) en fct de p

    # Calcul des Quasi niveaux Ef_c,v et des f_c,v(E)
    for i in range(len(n)):
        Efc[i]=interpolate.splev(n[i], tckE, der=0)
        Efv[i]=interpolate.splev(p[i], tckP, der=0)
        fc[:,i]=f(E, Efc[i] ,T);
        fv[:,i]=f(E, Efv[i], T)
        labellegend[i] = "{:10.4g} cm-3".format( delta_n[i]*1e-2**3 ) #cm^-3
    # --------- fin methode intégrale


    # INTERACTIONS AVEC DES PHOTONS
    nu = np.linspace(1.42, 1.50, num=10000)*e/h #nu (nm)
    # Énergies concernées par la transition
    E2 = m.Ec + m.mr/m.mc*(h*nu-m.Eg); E1=E2-h*nu
    # Densité d'état jointe
    rho_j= np.nan_to_num(np.real((2*m.mr)**(3/2)*(h*nu-m.Eg)**(1/2)/(pi*hbar**2)))

    gainMax = zeros(len(n)); gainMaxLoc=zeros(len(n), dtype=int)
    gamma0=zeros([len(nu),len(n)]);
    fg_vec=zeros([len(nu),len(n)]); fc=zeros([len(nu),len(n)]); fv=zeros([len(nu),len(n)]);
    for i in range(len(n)):
        #Fermi
        fc[:,i]=f(E2,Efc[i], T=T)
        fv[:,i]=f(E1,Efv[i], T=T)
        fg_vec[:,i]=fc[:,i]-fv[:,i]
        # Gain
        gamma0[:,i] = real(((c/nu)/l.n_GaAs)**2/(8*pi*l.tau_r)*rho_j*fg_vec[:,i])
        # Calcul du gain max
        index = argmax(gamma0[:,i]); gainMaxLoc[i] = index;
        gainMax[i]= gamma0[index, i]


    # Trouver la courbe du gain max @ 850 nm
    nu0=c/l.lambda0
    loc_gain850 = np.argmin(absolute(nu[gainMaxLoc] - nu0)) #Index de la courbe de gain max
    gamma850 = gamma0[loc_gain850]


    ## FIGURE GAIN
    color = plt.cm.coolwarm(np.linspace(0, 1,len(delta_n)))
    mpl.rcParams['axes.prop_cycle'] = plt.cycler('color', color)
    fig=plt.figure()
    plt.plot(nu*h/e, gamma0*1e-2) #cm-1
    plt.plot(nu[gainMaxLoc]*h/e, gainMax*1e-2, 'k.')
    plt.legend(labellegend)
    plt.axvline(x=h*c/(l.lambda0*e), linestyle = ':', color='black', linewidth=0.75)
    plt.xlabel('Énergie [eV]')
    plt.ylabel('Gain $\gamma$ [cm$^{-1}$]')
    plt.ylim(bottom=-30)
    plt.show()
    plt.savefig('fig/gain.eps', format='eps')

    return gainMax

# test:
# calcul_gain(delta_N, laser, mat, temp)











#
