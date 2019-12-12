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

def calcul_gain(delta_n, laser, mat, T, vis=False, leg=None, loc850=None):
    mat.updateT(T)

    # ELECTroNNNNNNNS --------------------------------------------------------
    ## Calcul des populations a l'equilibre
    ni = sqrt(4*(2*pi*kb*T/h**2)**3*(mat.mc*mat.mv)**(3/2)*exp(-mat.Eg0/(kb*T)));
    if mat.doping > 0:
       mat.p0 = mat.doping
       mat.n0 = ni**2/mat.p0
    elif mat.doping < 0:
       mat.n0 =-mat.doping;
       mat.p0 = ni**2/mat.n0;
    else:
       mat.n0 = ni
       mat.p0 = ni
    # Ajout de la population de pompe
    n=mat.n0+delta_n; p=mat.p0+delta_n

    # Pour plus tard
    E = np.linspace(-1, 2, num=100)*e #nu (nm) (Plage d'Énergie pour tests)
    Efc=zeros([len(n)]); Efv=zeros([len(n)]); # Quasi niveaux de Fermi
    fc=zeros([len(E),len(n)]); fv=zeros([len(E),len(n)]); # fv et fc
    labellegend=[None]*len(n)

    #---------------------- methode integrale
    # Trouver le quasi niveau de fermi (ELECTRONS) pour tout n:
    ## Resoudre integrale de fermi pour obtenir FD(eta)
    #Range different pour les electrons et les trous !
    Efc_test=linspace(0.9, 1.2, num=1000)*mat.Ec; Efv_test=linspace(0.1, -0.3, num=1000)*mat.Ec;
    etaE = (Efc_test -mat.Ec)/(kb*T); etaP =(mat.Ev-Efv_test)/(kb*T)
    fd_etaE, errE= fdelec(etaE) #intégration numérique
    # fd_etaE, errE= fdpolylog(etaE)
    fd_etaP, errP= fdelec(etaP) #intégration numérique
    # fd_etaP, errP= fdpolylog(etaP)

    # Obtenir relations pour n et quasi niveau
    gc=(2*mat.mc)**(3/2)/(2*pi**2*hbar**3); Nc=(kb*T)**(3/2)*gc; nE=Nc*fd_etaE;
    gv=(2*mat.mv)**(3/2)/(2*pi**2*hbar**3); Nv=(kb*T)**(3/2)*gv; nP=Nv*fd_etaP

    efE =mat.Ec + kb*T*etaE; efP =mat.Ev - kb*T*etaP
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
    nu = np.linspace(1.42, 1.50, num=1000)*e/h #nu (nm)
    # Énergies concernées par la transition
    E2 =mat.Ec +mat.mr/mat.mc*(h*nu-mat.Eg); E1=E2-h*nu
    # Densité d'état jointe
    rho_j= np.nan_to_num(np.real((2*mat.mr)**(3/2)*(h*nu-mat.Eg)**(1/2)/(pi*hbar**2)))

    gainMax = zeros(len(n)); gainMaxLoc=zeros(len(n), dtype=int)
    gamma0=zeros([len(nu),len(n)]);
    fg_vec=zeros([len(nu),len(n)]); fc=zeros([len(nu),len(n)]); fv=zeros([len(nu),len(n)]);
    for i in range(len(n)):
        #Fermi
        fc[:,i]=f(E2,Efc[i], T=T)
        fv[:,i]=f(E1,Efv[i], T=T)
        fg_vec[:,i]=fc[:,i]-fv[:,i]
        # Gain
        gamma0[:,i] = real(((c/nu)/laser.n_GaAs)**2/(8*pi*laser.tau_r)*rho_j*fg_vec[:,i])
        # Calcul du gain max
        index = argmax(gamma0[:,i]); gainMaxLoc[i] = index;
        gainMax[i]= gamma0[index, i]

    ## FIGURE GAIN
    if vis:
        color = plt.cm.winter(np.linspace(0, 1,len(delta_n)))
        mpl.rcParams['axes.prop_cycle'] = plt.cycler('color', color)
        fig=plt.figure()
        if loc850:
            plt.plot(nu*h/e, gamma0[:,0:loc850]*1e-2, '--') #cm-1
            plt.plot(nu*h/e, gamma0[:, loc850]*1e-2) #cm-1
            if loc850 < gamma0.shape[1]-1:
                plt.plot(nu*h/e, gamma0[:,loc850+1:]*1e-2, '--') #cm-1
        plt.plot(nu[gainMaxLoc]*h/e, gainMax*1e-2, 'k.')
        if leg:
            plt.legend(leg)
        else:
            plt.legend(labellegend)
        plt.axvline(x=h*c/(laser.lambda0*e), linestyle = ':', color='black', linewidth=0.75)
        plt.xlabel('Énergie [eV]')
        plt.ylabel('Gain $\gamma$ [cm$^{-1}$]')
        plt.ylim(bottom=-30)
        # plt.savefig('fig/gain.eps', format='eps')
        plt.show()
        plt.close()

    return gainMax, gainMaxLoc, gamma0, nu




def injection(I, laser):
    """ Densité de porteurs de charge injectée"""
    return laser.eta*I*laser.tau/(e*laser.vol)


def penis(vec_I, I150, laser, mat, T):
    #Fait juste calculer le gain
    loc=find_nearest(vec_I, I150)
    delta_n = injection(vec_I, laser)
    labellegend=[None]*len(vec_I)
    for i in range(len(vec_I)):
        labellegend[i] = "{:10.4g} mA".format( vec_I[i]*1e3 ) #mA
    gainMax, gainMaxLoc, gamma0, nu  = calcul_gain(delta_n, laser, mat, T, vis=True, leg=labellegend, loc850=loc)

    return gainMax


def Lorentienne(dv, v0, v):
    L= 1/pi*( dv/2)*((v-v0)**2 + (dv/2)**2)**(-1)
    return L/L[L.argmax()]
    # return 1*((v-v0)**2 + 1)**(-1)


def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return idx


#
