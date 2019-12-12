from numpy import pi, exp, log
from const import *


class Materiau:
    def __init__(self, ratio_mc, ratio_mv, Eg0, A, B, T, doping):
        self.mv=ratio_mv*m0 #Masses effective trous [kg]
        self.mc=ratio_mc*m0 #Masses effective electrons [kg]
        self.mr=1/(1/self.mc+1/self.mv) #[kg]
        # Energies
        self.A = A
        self.B= B
        self.Eg0=Eg0*e                    #Gap @t=0 [J]
        self.Eg=(self.Eg0 - (A*e*T**2)/(B+T)) #E_gap(T) [J]
        self.Ev=0                       #Bande de valence [J]
        self.Ec=self.Ev+self.Eg         #Bande de conduction [J]
        # densite de porteurs [m^-3]
        self.doping=doping
        self.n0=0
        self.p0=0

    def updateT(self, temp):
        self.Eg=(self.Eg0 - (self.A*e*temp**2)/(self.B+temp))

class Laser:
    def __init__(self, caca, pipi):
        # Parametres de l'enonce
        self.tau=2e-9 # [s]Temps recombinaison e-trous
        self.eta = 0.99 # Efficacite quantique interne
        self.tau_r=self.tau/self.eta
        self.lambda0=850e-9 # [m] Longeur d'onde d'operation
        self.nu0=c/self.lambda0
        self.R1=0.99 #Reflectivite du mirroir R1
        self.beta=1e-4 #Facteur d'emission spontanee
        self.alpha_int=10*1e2 #[m-1] Pertes internes
        self.n_GaAs=3.655
        self.n_AlGaAs=3.3620
        self.R2=((self.n_GaAs-1)/(self.n_GaAs+1))**2
        # self.R2=0.8

        #Dimensions
        # f=0.9
        self.l=caca#longueur
        # self.lambda0/(self.n_GaAs)
        self.w=pipi#largeur
        self.d=2*1.45e-7
        self.vol=self.d*self.w*self.l
        #Intervalle spectral libre
        self.isl = c/(2*self.l*self.n_GaAs)
        # Efficacité quantique externe
        self.alpha_m1 = 1/(2*self.l)*log(1/self.R1)
        self.alpha_m2 = 1/(2*self.l)*log(1/self.R2)
        self.alpha_r = self.alpha_m1 + self.alpha_m2 + self.alpha_int
        self.eta_ext = self.eta*self.alpha_m2/self.alpha_r

        # Facteur de confinement
        G0 = 2*pi**2*(self.n_GaAs**2 - self.n_AlGaAs**2)*(self.d/850e-9)**2
        self.Gamma=G0/(1+G0)

        # Angle de divergence
        self.theta_par = self.lambda0/self.w #parallèle
        self.theta_perp = self.lambda0/self.d #perpendiculaire

        # Temps de vie photon :
        self.tau_p=self.n_GaAs/(c*self.alpha_r)
        self.tau_m=self.n_GaAs/((self.alpha_m1+self.alpha_m2)*c)

        # Transparency @ T=25 C:
        self.B=None
        self.delta_n_tr=None

        # Threshold @ T=25C
        self.I_th = None
        self.delta_n_th=None
        self.gamma_th=None

        # Efficacité de pente:
        self.eff_pente=None
