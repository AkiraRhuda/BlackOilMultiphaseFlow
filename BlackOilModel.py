import numpy as np
import unitsconverter
import warnings

class GasPhase:
    """
    Classe para calcular compressibilidade, viscosidade e massa específica na fase gás.

    Parâmetros:
    ----------
    zcorrelation : str
        Informa a correlação a ser usada para calcular o fator Z - suporta somente Papay (1985)
    μcorrelation : str
        Informa a correlação a ser usada para calcular a viscosidade do gás - suporta Lee et al (1966) e Dempsey (1965)
    P : float
        Pressao absoluta em psia
    Ppc : float
        Pressao pseudocritica em psia
    Ppr : float
        Pressao pseudoreduzida em psia
    T : float
        Temperatura em R
    Tpc : float
        Temperatura pseudocritica em R
    Tpr : float
        Temperatura pseudoreduzida em R
    """
    
    def __init__(self, zcorrelation, μcorrelation, dg=None,P=None, Ppc=None, Ppr=None, T=None, Tpc=None, Tpr=None):
        if dg is None:
            raise ValueError("Densidade relativa do gás (dg) precisa ser informada.")
        else:
            self.dg = dg
        self.P, self.Ppc, self.Ppr, self.T, self.Tpc, self.Tpr = P, Ppc, Ppr, T, Tpc, Tpr
        self.zcorrelation, self.μcorrelation = zcorrelation, μcorrelation
        self.Cpr, self.Cg, self.Z, self.ρ = None, None, None, None
        self.initialproperties()
        self.zfactcorrelselector()
        self.ρ_g()
        self.μfactcorrelselector()
        self.compressisoterm()

    # initial properties #

    def initialproperties(self):
        if self.dg is not None and self.Ppc is None and self.Tpc is None:
            if self.dg < 0.75:
                self.Ppc = 677 + 15*self.dg - 37.5*self.dg**2
                self.Tpc = 168 + 325*self.dg - 12.5*self.dg**2
            else:
                self.Ppc = 706 - 51.7*self.dg - 11.1*self.dg**2
                self.Tpc = 187 + 330*self.dg - 71.5*self.dg**2
        if self.Ppr is None and self.Tpr is None:
            self.Ppr = self.P/self.Ppc
            self.Tpr = self.T/self.Tpc
        if self.dg is None:
            raise ValueError("Densidade relativa do gás (dg) precisa ser informada.")
        self.Mg = self.dg * 28.96 # lb / lbmol
        self.Psc = 14.696  # psi
        self.Tsc = 519.67  # °R (60°F)
    
    def ρ_g(self):
        R = 10.73 # psia·ft³/ (lb·mol·°R)
        self.ρ = (self.P * self.Mg)/(self.Z * R * self.T)
        
    def Bg(self):
        """
        Retorna Bg em ft³/scf.

        P e Psc em psi.
        T e Tsc em Rankine.
        """
        return self.Psc/(self.Tsc) * self.Z * self.T/self.P # SERA QUE É EM FAHRENHEIT MESMOOO ?????????? unitsconverter.Temperature(self.T, 'R', 'F')
        
    def zfactcorrelselector(self):
        if self.zcorrelation == "Papay":
            self.papay()
        else:
            raise Exception('Correlação não encontrada')

    def μfactcorrelselector(self):
        if self.μcorrelation == "Lee":
            self.lee()
        elif self.μcorrelation == "Dempsey":
            self.dempsey()
        else:
            raise Exception('Correlação não encontrada')
            
    # Z factor correlations #

    def papay(self, Ppr=None):
        if Ppr is None:
            self.Z = 1 - (3.53*self.Ppr)/(10**(0.9813*self.Tpr)) + (0.274*self.Ppr**2)/(10**(0.8157*self.Tpr))
        else:
            return 1 - (3.53*Ppr)/(10**(0.9813*self.Tpr)) + (0.274*Ppr**2)/(10**(0.8157*self.Tpr))

    # μ correlations #
    
    def lee(self):
        x_v = 3.448 + (986.4/self.T) + (0.01009*self.Mg)
        y_v = 2.4 - (0.2*x_v)
        k_v = ((9.379+0.0160*self.Mg) * (self.T**1.5)) / (209.2 + (19.26*self.Mg) + self.T)

        factor = x_v * ((self.ρ/62.4)**y_v)
        self.μ = 1e-4 * k_v * np.exp(factor)
        
    def dempsey(self):
        self.μΔgN2, self.μΔgCO2, self.μΔgH2S = 0, 0, 0

        μgncorrigida = (1.709e-5 - 2.062e-6*self.dg)*unitsconverter.Temperature(self.T, 'R', 'F') + 8.188e-3-6.15e-3*np.log10(self.dg)
        
        self.μg = μgncorrigida + self.μΔgN2 + self.μΔgCO2 + self.μΔgH2S
        
        a0, a1, a2, a3, a4, a5 = -2.4621, 2.9705, -2.8626e-1, 8.0542e-3, 2.8086, -3.4980
        a6, a7, a8, a9, a10, a11 = 3.6037e-1, -1.0443e-2, -7.9339e-1, 1.3964, 1.4914e-1, 4.4102e-3
        a12, a13, a14, a15 = 8.3939e-2, -1.8641e-1, 2.0336e-2, -6.0958e-4

        lnexpression = a0 + a1*self.Ppr + a2*self.Ppr**2 + a3*self.Ppr**3 + self.Tpr*(a4 + a5*self.Ppr + a6*self.Ppr**2
                    + a7*self.Ppr**3) + self.Tpr**2*(a8 + a9*self.Ppr + a10*self.Ppr**2 +a11*self.Ppr**3) \
                    + self.Tpr**3*(a12 + a13*self.Ppr + a14*self.Ppr**2 + a15*self.Ppr**3)
                    
        self.μ = np.exp(lnexpression)*self.μg/self.Tpr

    def compressisoterm(self):
        Ppr = self.Ppr
        dx = 0.5 * 10 ** (-6)
        self.Cpr = 1/self.Ppr - 1/self.Z * (self.papay(Ppr+dx)-self.papay(Ppr))/dx
        self.Cg = self.Cpr/self.Ppc

    def output(self):
        # print(f'Massa específica do gás: {self.ρ:.4f} lbm/ft³')
        # print(f'Compressibilidade do gás: {self.Cg:.4e} psi⁻¹')
        # print(f'Viscosidade do gás: { self.μ:.4f} cp')
        return self # self.ρ, self.Cg, self.μ

class OilPhase:
    """
    Classe para calcular compressibilidade, viscosidades e massa específica na fase óleo.

    Parâmetros:
    ----------
    soluratiocorrelation : str
        Informa a correlação a ser usada para calcular a razão de solubilidade - suporta somente Standing (1947)
    oilvolcorrelation : str
         Informa a correlação a ser usada para calcular a compressibilidade isotérmica do óleo e o fator volume-formação do óleo - suporta somente Standing (1947)
    oilviscosityselector : str
        Informa a correlação a ser usada para calcular as viscosidades do óleo morto e saturado - suporta somente Brill e Beggs (1974)
    do : float
        densidade relativa do óleo
    dg : float
        densidade relativa do gás
    P : float
        Pressao absoluta
    Pb : float
        Pressao de bolha
    T : float
        Temperatura
    units : list
        Lista com unidades de [pressão, temperatura], respectivamente
    """
    def __init__(self, soluratiocorrelation, oilvolcorrelation, oilviscosityselector, do=None, dg=None, P=None, Pb=None, T=None, RGL=None, BSW=None):
        
        self.do, self.dg, self.P, self.Pb, self.T = do, dg, P, Pb, T
        self.API, self.Rs, self.Bob, self.ρob = None, None, None, None
        self.soluratiocorrelation, self.oilvolcorrelation, self.oilviscosityselector = soluratiocorrelation, oilvolcorrelation, oilviscosityselector

        self.initialproperties()
        self.soluratio()
        if Pb is None:
            self.RGO = RGL/(1-BSW)
            self.standingbubblepress()
        self.oilvolforcorrelationselector()
        self.ρ_o()
        self.oilvisccorrelationselector()
        self.oilclassifier()

    def initialproperties(self):
        self.API = 141.5/self.do - 131.5
        
    def standingbubblepress(self):
        self.Pb = 18.2 * ((self.Rs/self.dg)**0.83 * 10**(0.00091 * self.T - 0.0125 * self.API) - 1.4)

    # Solubility ratio #
    def soluratio(self):
        if self.Pb is not None and self.P > self.Pb:
            self.Rs = self.dg * (((self.Pb/18.2) + 1.4) * 10**(0.0125*self.API - 0.00091*self.T))**(1/0.83) # self.Rs = self.Rsb
            
        else:
            self.solurationcorrelationselector()

    def solurationcorrelationselector(self):
        if self.soluratiocorrelation == 'Standing':
            self.standingsoluration()
    
    def standingsoluration(self):
        # if any(isinstance(self.dg, float)):
     self.Rs = self.dg * (((self.P/18.2) + 1.4) * 10**(0.0125*self.API - 0.00091*self.T))**(1/0.83)
    
    # Oil volume formation factor #
    def oilvolforcorrelationselector(self):
        if self.oilvolcorrelation == "Standing":
            self.standingdRsdP()
            self.standingdBodP()
            # condição....
            if self.P >= self.Pb: #### >= AGORA
                self.standingBob()
                self.standingCo()
                self.standingBo()
            else:
                self.standingBo()
                self.standingCo()

    def standingdRsdP(self):
        A = 0.0125*self.API - 0.00091*self.T
        B = self.P/18.2 + 1.4
        n = 1/0.83
        self.dRsdP =  self.dg * n * (10**A*B)**(n-1)*10**A/18.2
        # self.dRsdP = self.dg * (10**A*B)**n
        # self.dRsdP = self.dg * (((1/18.2) + 1.4) * 10**(0.0125*self.API - 0.00091*self.T))**(1/0.83) # OLHAR AQUI!!!!!!!


    def standingBob(self):
        self.Bob = 0.9759 + 0.00012 * (self.RGO * (self.dg / self.do) ** (0.5) + 1.25 * self.T) ** (1.2)

    def standingBo(self):
        if self.P > self.Pb + 1e-12:
            self.Bo = self.Bob * np.exp(-self.Co * (self.P - self.Pb))
        else:
            self.Bo = 0.9759 + 0.00012*(self.Rs * (self.dg/self.do)**(0.5) + 1.25*self.T)**(1.2)
            
    def standingdBodP(self):
        A = self.Rs * (self.dg / self.do) ** (0.5) + 1.25 * self.T
        self.dBodP = 0.00012 * 1.2 * A**0.2 * (self.dg / self.do) ** (0.5) * self.dRsdP
        # self.Bo = 0.9759 + 0.00012 * (self.Rs * (self.dg / self.do) ** (0.5) + 1.25 * self.T) ** (1.2)
        # self.dBodP = 0.00012*(self.dRsdP * (1.2)*(self.dg/self.do)**(0.5) + 1.25*self.T)**(0.2) # OLHAR AQUI!!!!!!!
    
    def standingCo(self):
        if self.P >= self.Pb: # dando erro no rho ob se isso e true pq n tem bo
            self.ρ_ob() # newton raphson????????????????
            num = self.ρob + 0.004347 * (self.P - self.Pb) - 79.1
            dem = 0.0007141 * (self.P - self.Pb) - 12.938
            self.Co = 1e-6 * np.exp(num/dem)
        else:
            self.Co = -1/self.Bo * self.dBodP + GasPhase("Papay", "Lee", dg=self.dg, P=self.P, T=unitsconverter.Temperature(self.T, 'F', 'R')).Bg()/self.Bo * self.dRsdP

    # Oil density #
    def ρ_ob(self):
        self.ρob = (62.4 * self.do + 0.0136 * self.RGO * self.dg) / self.Bob
    
    def ρ_o(self):
        if self.P > self.Pb:
            self.ρo = self.ρob * np.exp(self.Co*(self.P-self.Pb))

        else:
            self.ρo = (62.4 * self.do + 0.0136 * self.Rs * self.dg)/self.Bo

    # Oil viscosity #
    def oilvisccorrelationselector(self):
        if self.oilviscosityselector == "Beggs":
           self.μbeggsod() 
           self.μbeggsob()
    
    def μbeggsod(self):
        self.μod = (10 ** ((10 ** (3.0324 - 0.02023 * self.API)) * self.T ** (-1.163))) - 1
    
    def μbeggsob(self):
        coeff = 10.715*(self.Rs+100)**(-0.515)
        exp = 5.44*(self.Rs+150)**(-0.338)
        self.μob = coeff*self.μod**exp

    def oilclassifier(self):
        if self.API < 40:
            typee = "Black-Oil"
        elif 40 <= self.API <= 45:
            if self.Bo >= 2:
                typee = "Óleo volátil"
            else:
                typee = "Necessário determinar RGO para caracterizar o fluído"
        elif 45 < self.API <= 55:
            typee = "Óleo volátil"

        else:
            typee = "Modelo não aplicável - Óleo não é Black-Oil ou óleo volátil"
            raise Exception(typee)
        # print('* Tipo de fluido: ', typee, '*')
    
    def output(self):
        # print(f'Massa específica do óleo: {self.ρo:.2f} lbm/ft³')
        # print(f'Compressibilidade do óleo: {self.Co:.4e} psi⁻¹')
        # print(f'Viscosidade do óleo morto: {self.μod:.3f} cp')
        # print(f'Viscosidade do óleo saturado: {self.μob:.3f} cp')
        return self # self.ρo, self.Co, self.μod, self.μob

class WaterPhase:
    """
    Classe para calcular compressibilidade, viscosidades e massa específica da fase água.

    Parâmetros:
    ----------
    P : float
        Pressao absoluta em psia
    T : float
        Temperatura em fahrenheit
    S : float
    Salinidade da água em % de sólidos
    """
    def __init__(self, P=None, T=None, S=None):
        self.P = P
        self.T = T
        self.S = S
        self.ρw = self.ρ_wsc()
        self.Rswpure = self.soluratiopurewater()
        self.Rsw = self.soluratiobrinewater()
        self.Cw = self.compressisoterm()
        self.Bw = self.watervolfor()
        self.μwatm = self.viscosityatm()
        self.μw = self.viscosity()

    def ρ_wsc(self):
        return 62.368 + 0.438603*self.S + 1.60074e-3 * self.S**2 # water density lb/SCF

    def soluratiopurewater(self):
        # temperature in °F
        # pressure psia
        A = 8.15839 - 6.12265e-2 * self.T + 1.91663e-4 * self.T**2 - 2.1654e-7 * self.T**3
        B = 1.01021e-2 - 7.44241e-5 * self.T + 3.05553e-7 * self.T**2 - 2.94883e-10 * self.T**3
        C = (-9.02505 + 0.130237*self.T - 8.53425e-4 * self.T**2 + 2.34122e-6 * self.T**3 - 2.37049e-9 * self.T**4)*10**(-7)
        return A + B * self.P + C * self.P**2

    def soluratiobrinewater(self):
        return 10**(-0.0840655*self.S*self.T**(-0.285854))*self.Rswpure

    def compressisoterm(self):
        # Osif, T.L (1988)
        # temperature in °F
        # pressure psia
        # S mg/L
        # LIMITES 0 a 200000 mg/L de salinidade ; 1000 a 20000 psia ; temperatura 200 a 270°F
        if 200 <= self.T <= 270 and 1000 <= self.P <= 20000:
            S = unitsconverter.Density([self.S, unitsconverter.Density(self.ρw, 'lbm/ft3', 'do')], '%', 'mg/L')
            Cw = 1/(7.033*self.P + 0.5415*S - 537.0*self.T + 403300) # psia-1
        else:
            S = unitsconverter.Density([self.S,unitsconverter.Density(self.ρw,'lbm/ft3','do')],'%', 'mg/L')
            Cw = 1 / (7.033 * self.P + 0.5415 * S - 537.0 * self.T + 403300)  # psia-1
            warnings.warn('Parâmetro(s) da correlação para compressibilidade da água está(ão) fora do(s) limite(s)!')
        return Cw


    def watervolfor(self):
        # temperature in °F
        # pressure psia
        # Bw em bbl/STB
        ΔVwt = -1.0001e-2 + 1.33391e-4*self.T + 5.50654e-7*self.T**2
        ΔVwp = -1.95301e-9*self.P*self.T - 1.72834e-13*self.P**2*self.T - 3.58922e-7*self.P - 2.25341e-10*self.P**2
        return (1+ΔVwt)*(1+ΔVwp)

    def viscosityatm(self):
        # temperature in °F
        # S %
        # μwatm in cP
        A = 109.527 - 8.40564*self.S + 0.313314*self.S**2 + 8.72213*self.S**3
        B = -1.12166 + 2.63951e-2*self.S - 6.79461e-4*self.S**2 - 5.47119e-5*self.S**3 - 1.55586e-6*self.S**4
        return A * self.T**B

    def viscosity(self):
    # pressure psia
        return (0.9994 + 4.0295e-5*self.P + 3.1062e-9*self.P**2)*self.μwatm


    def output(self):
        # print(f'Massa específica da água: {self.ρw:.4f} lbm/ft³')
        # print(f'Compressibilidade da água: {self.Cw:.4e} psi⁻¹')
        # print(f'Viscosidade da água: { self.μw:.4f} cp')
        return self # self.ρw, self.Cw, self.μw