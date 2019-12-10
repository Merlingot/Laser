from numpy import pi, exp, log
hbar = 1.05457180014e-34
c=3e8
kb = 1.38064852e-23
me = 9.10938356e-31
h=hbar*2*pi
e=1.60218e-19
m0=9.109e-31 #kg (Free electron mass)



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
    def __init__(self):
        # Parametres de l'enonce
        self.tau=2e-9 # [s]Temps recombinaison e-trous
        self.eta = 0.99 # Efficacite quantique interne
        self.tau_r=self.tau/self.eta
        self.lambda0=850e-9 # [m] Longeur d'onde d'operation
        self.nu0=c/self.lambda0
        self.R1=0.99 #Reflectivite du mirroir R1
        self.beta=1e-4 #Facteur d'emission spontanee
        self.alpha_int=10/1e-2 #[m-1] Pertes internes
        self.n_GaAs=3.67
        self.n_AlGaAs=3.4117 #À vérifier
        self.R2=((self.n_GaAs-1)/(self.n_GaAs+1))**2
        #Dimensions
        self.d=1e-7 #llongueur
        self.w=1e-7 #largeur
        self.l=0.2e-6 #épaisseur
        self.vol=self.d*self.w*self.l

        # Efficacité quantique externe
        alpha_m1 = 1/(2*self.l)*log(1/self.R1)
        alpha_m2 = 1/(2*self.l)*log(1/self.R2)
        self.alpha_r = alpha_m1 + alpha_m2 + self.alpha_int
        self.eta_ext = self.eta*alpha_m2/self.alpha_r

        # Facteur de confinement
        G0 = 2*pi**2*(self.n_GaAs**2 - self.n_AlGaAs**2)*(self.l/850e-9)**2
        self.Gamma0=G0/(1+G0)
