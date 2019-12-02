from numpy import pi, exp
hbar = 1.05457180014e-34
c=3e8
kb = 1.38064852e-23
me = 9.10938356e-31
h=hbar*2*pi
e=1.60218e-19
m0=9.109e-31 #kg (Free electron mass)
T=25+273.15 # [K]



class Materiau:
    def __init__(self, ratio_mc, ratio_mv, Eg0, A, B, T, doping):
        self.mv=ratio_mv*m0 #Masses effective trous [kg]
        self.mc=ratio_mc*m0 #Masses effective electrons [kg]
        self.mr=1/(1/self.mc+1/self.mv) #[kg]
        # Energies
        self.Eg0=Eg0*e                    #Gap @t=0 [J]
        self.Eg=(self.Eg0 - (A*e*T**2)/(B+T)) #E_gap(T) [J]
        self.Ev=0                       #Bande de valence [J]
        self.Ec=self.Ev+self.Eg         #Bande de conduction [J]
        # densite de porteurs [m^-3]
        self.doping=doping
        self.n0=0
        self.p0=0

class Laser:
    def __init__(self):
        # Parametres de l'enonce
        self.tau=2e-9 # [s]Temps recombinaison e-trous
        self.eta = 0.99 # Efficacite quantique interne
        self.tau_r=self.tau/self.eta
        self.lambda0=850e-9 # [m] Longeur d'onde d'operation
        self.nu0=c/self.lambda0
        self.R=0.99 #Reflectivite du mirroir R1
        self.beta=1e-4 #Facteur d'emission spontanee
        self.pertes=10/1e-2 #[m-1] Pertes internes
        self.n_ref=1.2640 #GaAs
        #Dimensions
        self.d=1e-7; self.w=5e-6; self.L=5e-6
        self.Vol=self.d*self.w*self.L
