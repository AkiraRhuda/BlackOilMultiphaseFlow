import numpy as np
import select
import math
import unitsconverter
import MixtureTemperature


class MultiPhaseFlowModel:
    """
    Classe para calcular as propriedades do escoammento do fluido multifásico (fase óleo, água e gás).

    Parâmetros:
    ---------
    model : list
        Informa o tipo de modelo a ser utilizado para calcular as propriedades do fluido em cada trecho do poço.
        Suporta Homogêneo e Drift-Flux
        Se informar somente um modelo, este será utilizado em tod0 o poço.
    temperature : list
        Temperatura do fluido em K para cada trecho
    pressure : float
        Pressão do fluido em Pa
    reference_pressure : float
        Pressão de referência para calcular um erro ao final usando a queda de pressão total do sistema em Pa
    diameters : list
        Lista de diâmetros dos tubos
    rugosity : float
        Rugosidade da tubulação
    incli : float
        Inclinação do tubo em graus
    rho_l : float
        Massa específica do fluido em kg/m3
    rho_o : float
        Massa específica do óleo em kg/m3
    rho_g : float
        Massa específica do gás em kg/m3
    rho_w : float
        Massa específica da água em kg/m3
    do : float
        Densidade do óleo em relação à água
    dg : float
        Densidade do gás em relação ao ar
    salinity : float
        Salinidade da água em % (teor de sólidos)
    QLsc : float
        Vazão de líquidos na condição standard
    QOsc : float
        Vazão de óleo na condição standard
    QWsc : float
        Vazão de água na condição standard
    QL : float
        Vazão de líquidos in-situ (condição de reservatório)
    QO : float
        Vazão de óleo in-situ (condição de reservatório)
    QW : float
        Vazão de água in-situ (condição de reservatório)
    Bo : float
        Fator volume-formação de óleo em sm3/sm3
    Bw : float
        Fator volume-formação de água em sm3/sm3
    Bg : float
        Fator volume-formação de gás em sm3/sm3
    Bsw : float
        Teor de água e sedimentos no líquido em %
    RGL : float
        Razão gás-líquido em sm3/sm3
    RGO : float
        Razão gás-óleo em sm3/sm3
    Rs : float
        Razão de solubilidade do gás no óleo em sm3sm3
    Rsw : float
        Razão de solubilidade do gás na água em sm3sm3
    mu_o : float
        Viscosidade do óleo em Pa.s
    mu_w : float
        Viscosidade da água em Pa.s
    mu_g : float
        Viscosidade do gás em Pa.s
    mu_l : float
        Viscosidade do líquido em Pa.s
    sig_og : float
        Tensão óleo-gás em N/m
    sig_og : float
        Tensão água-gás em N/m
    sig_lg : float
        Tensão líquido-gás em N/m
    Z : float
        Fator Z para correção do modelo de gases ideais para gases reais
    length : list
        Comprimentos das tubulações de toda a arquitetura do sistema em m
    """

    def __init__(self, model : list, temperature, pressure, diameters : list, rugosity, incli : list, rho_l = None, rho_o = None, rho_g = None, rho_w = None,
            do = None, dg=None, salinity = None,QLsc = None, QOsc = None, QWsc = None, QL = None, QO = None, QW = None, QG=None, Bo = None, Bw = None, Bg = None,
            Bsw = None,RGL = None, RGO = None, Rs = None, Rsw = None, mu_o = None, mu_w = None, mu_g = None, mu_l = None, sig_og = None,
            sig_wg = None, sig_lg = None, Z = None, length : list = None, N_division=None, reference_pressure = None, postprocess=False):
        self.model = model
        self.temperature = [temperature]
        self.initial_pressure = pressure
        self.Pressure = [pressure]
        self.rho_o = rho_o
        self.rho_g = rho_g
        self.rho_w = rho_w
        self.rho_l = rho_l # caso ele forneca direto esse valor
        self.do = do
        self.dg = dg
        self.salinity = salinity
        self.diameters = diameters
        self.rugos = rugosity
        self.incli = incli
        self.QLsc = QLsc
        self.QOsc = QOsc
        self.QWsc = QWsc
        self.QL = QL
        self.QO = QO
        self.QW = QW
        self.QG = QG
        self.Bo = Bo
        self.Bw = Bw
        self.Bg = Bg
        self.Bsw = Bsw/100 if Bsw is not None else 0
        self.RGL = RGL
        if RGO is not None and self.RGL is None and self.Bsw != 0:
            self.RGL = RGO*(1-self.Bsw)
        self.RGO = RGO
        self.Rs = Rs
        self.Rsw = Rsw
        self.mu_o = mu_o
        self.mu_w = mu_w
        self.mu_g = mu_g
        self.mu_l = mu_l # caso ele forneca direto esse valor
        self.sig_og = sig_og
        self.sig_wg = sig_wg
        self.sig_lg = sig_lg # caso ele forneca direto esse valor
        self.Z = Z
        self.length = length
        self.N_division = N_division
        self.Pref = reference_pressure
        self.postprocess = postprocess
        self.exception = False
        self.Qmassica = []

        if (self.Bo is None and self.Bw is None or self.Bg is None or self.RGL is None  and self.RGO is None or self.Rs is None  and self.Rsw
                is None or self.mu_o is None and self.mu_w is None and self.mu_l is None or self.mu_g is None or self.sig_og is None  and
                self.sig_wg is None and self.sig_lg is None) and self.QL is None:
            if self.do is None:
                self.do = unitsconverter.Density(self.rho_o, 'kg/m3', 'do')
            if self.dg is None:
                self.dg = unitsconverter.Density(self.rho_g, 'kg/m3', 'dg')
            self.update_PVT_properties(unitsconverter.Temperature(self.temperature[-1],'K','C'), self.Pressure[-1], self.do, self.dg, self.RGL,self.Bsw, self.salinity)

        if self.Bw is not None and self.RGL is None or self.RGL is not None and self.Bw is None and self.QLsc is None:
            raise Exception('Verifique se está sendo fornecido RGL e Bw!')

        # if self.RGL is not None and self.RGO is not None and self.QLsc is None:
        #     raise Exception('Não é possível utilizar RGO e RGL simultaneamente!')

        if self.QLsc is not None and self.QOsc is not None and self.QWsc is not None:
            if self.QLsc != self.QOsc + self.QWsc:
                raise Exception('A vazão total de líquidos não corresponde à soma das vazões de água e óleo em condição standard')

        self.pipe_properties()
        self.architecture()
        # self.run()


    def update_PVT_properties(self, T, P, do, dg, RGL, BSW, salinidade):
        """
        T : °C
        P : Pa
        do : d
        dg : d
        RGL : sm3/sm3
        BSW : % em decimal
        salinidade : % em decimal
        """
        from BlackOilModel import GasPhase
        from BlackOilModel import OilPhase
        from BlackOilModel import WaterPhase
        Gas = GasPhase(zcorrelation="Papay", μcorrelation='Lee', dg=dg, P=unitsconverter.Pressure(P, 'Pa', 'psi'),
                       T=unitsconverter.Temperature(T, 'C', 'R')).output()
        Oil = OilPhase("Standing", "Standing", "Beggs", do=do, dg=dg, P=unitsconverter.Pressure(P, 'Pa', 'psi'),
                       T=unitsconverter.Temperature(T, 'C', 'F'), RGL=RGL, BSW=BSW).output()
        Water = WaterPhase(P=unitsconverter.Pressure(P, 'Pa', 'psi'),
                           T=unitsconverter.Temperature(T, 'C', 'F'), S=salinidade).output()
        self.rho_o = unitsconverter.Density(Oil.ρo, 'lbm/ft3','kg/m3')
        self.rho_g = unitsconverter.Density(Gas.ρ, 'lbm/ft3','kg/m3')
        self.rho_w = unitsconverter.Density(Water.ρw, 'lbm/ft3','kg/m3')

        self.Bo = unitsconverter.Standard(Oil.Bo, 'bbl/stb', 'sm3/sm3')
        self.Bw = unitsconverter.Standard(Water.Bw, 'bbl/stb', 'sm3/sm3')
        self.Bg =  unitsconverter.Standard(Gas.Bg(), 'scfstb', 'sm3sm3')


        #
        # self.Bsw = Bsw
        # self.RGL = RGL
        # self.RGO = RGO

        # self.Rs = Oil.Rs
        # self.Rsw = Water.Rsw
        self.Rs = unitsconverter.Standard(Oil.Rs, 'scfstb', 'sm3sm3') if Oil.Rs else 0
        self.Rsw = unitsconverter.Standard(Water.Rsw, 'scfstb', 'sm3sm3') if Water.Rsw else 0
        self.mu_o = unitsconverter.Viscosity(Oil.μob,'cp','Pa.s')
        self.mu_w = unitsconverter.Viscosity(Water.μw,'cp','Pa.s')
        self.mu_g = unitsconverter.Viscosity(Gas.μ,'cp','Pa.s')
        # self.run()
        self.flow(force_update=True)
        self.fractions()
        self.mixtureproperties(force_update=True)


    def run(self):
        self.flow()
        self.fractions()
        self.mixtureproperties()

        if isinstance(self.model, str):
            self.model = [self.model for _ in range(len(self.incli))]
        elif isinstance(self.model, list):
            if len(self.model) != len(self.incli):
                raise Exception("Verifique se há um modelo multifásico para cada trecho do poço!")
        else:
            raise Exception(f'{self.model} não é suportado!\nUse o modelo Homogêneo ou Drift-Flux')

        if self.length is not None:
            self.Pos = [0]
            self.VM = [0] # supondo velocidade inicial igual a 0 no início do poço
            actual_pos = 0
            actual_depth = 0
            for L, theta in zip(self.length, self.incli):
                actual_depth += L*np.sin(theta*np.pi/180)
            for L, theta, model in zip(self.length, self.incli, self.model):
                self.select_model(model)
                accumulated_L = 0
                while accumulated_L < L - 1e-4:
                    """if 1000 < actual_pos < 1099:
                        print('aqui')
                        print(actual_pos)
                        print(self.vm)
                        print(self.QG)
                        print(self.QL)
                        print(self.Bg)
                        print(self.Bo)
                        print(self.Pressure[-1])
                        print(self.Ap)
                        print(theta)
                        print("Bo", self.Bo)
                        print("QOsc", self.QOsc)
                        print("QO", self.QO)
                        print("Bw", self.Bw)
                        print("QW", self.QW)
                        print("QL", self.QL)"""
                    # O dl pode ser alterado caso o poço não seja totalmente percorrido pelo algoritmo, gerando um dl que completa o poço!!!
                    current_dl = self.dl
                    if accumulated_L + current_dl > L:
                        current_dl = L - accumulated_L
                    actual_depth -= current_dl * np.sin(theta * np.pi / 180)

                    PerdaCarga = self.model_instance.output(dl=current_dl, incli=theta)
                    Temp_Pressure = self.Pressure[-1] + float(PerdaCarga)

                    try:
                        Temp_Temperature = MixtureTemperature.Temperature(T=self.temperature[-1], P=Temp_Pressure, rho_o=self.rho_o,
                                            QO=self.QO, QW=self.QW, QL=self.QL, Mg=self.model_instance.Mg, HL=self.model_instance.H_L).Calculate_Temperature(
                            dl=current_dl, depth=actual_depth, m_dot_m=self.model_instance.mass_flow, theta=theta)
                        if Temp_Temperature!=Temp_Temperature:
                            Temp_Temperature = self.temperature[-1]
                            print(f'Erro no cálculo da temperatura na profundidade {actual_depth} m!\nReplicando a temperatura anterior para não quebrar o códgio')

                    except:
                        Temp_Temperature = self.temperature[-1]
                        print(f'Erro no cálculo da temperatura na profundidade {actual_depth} m!\nReplicando a temperatura anterior para não quebrar o códgio')
                    print(Temp_Temperature)
                    if Temp_Pressure <= 0 or Temp_Pressure != Temp_Pressure:
                        print(
                            f"\nPressão chegou a zero na posição L={actual_pos:.2f} m. O poço não tem energia para fluir até a superfície.")
                        self.exception = True
                        break

                    if isinstance(Temp_Pressure, float) is not True and isinstance(Temp_Pressure, int) is not True:
                        raise ValueError('Pressão assumiu not a number! Verifique a entrada de dados!')
                    self.update_PVT_properties(unitsconverter.Temperature(self.temperature[-1], 'K', 'C'),
                                                                  Temp_Pressure, self.do,
                                                                  self.dg,
                                                                  self.RGL,
                                                                  self.Bsw,
                                                                  self.salinity)
                    actual_pos += current_dl
                    accumulated_L += current_dl
                    self.temperature.append(Temp_Temperature)
                    self.Pressure.append(Temp_Pressure) # 44866301.289638236
                    self.VM.append(self.vm)
                    self.Pos.append(actual_pos)



                if self.exception is True:
                    break

            self.finalPressure = self.Pressure[-1]
            print(f'Pressão Calculada no final do poço: {self.finalPressure:.2f} Pa')

            if self.Pref is not None:
                if self.exception is False:
                    self.abserro = self.finalPressure - self.Pref
                    self.relerro = abs((self.abserro / self.Pref)) * 100

                    if abs(self.abserro) > 0:
                        print(f'Erro Absoluto: {self.abserro:.2f} Pa')
                        print(f'Erro Relativo: {self.relerro:.2f} %')
                    else:
                        raise ValueError('Erro assumiu not a number! Verifique a entrada de dados!')
                elif self.exception is True:
                    self.abserro = np.inf
                    self.relerro = np.inf

            self.PostProcess(self.Pressure, self.Pos,self.VM, self.temperature) # arrumar essas variaveis de entrada...

        else:
            for theta, model in zip(self.incli, self.model):
                self.select_model(model)
                self.model_instance.output(incli=theta)

    def select_model(self, model):
        if model.lower() == 'homogeneo' or model.lower() == 'homogêneo':
            self.model_instance = HomogeneousFlowModel(self)
        elif model.lower() == 'drift-flux' or model.lower() == 'driftflux':
            self.model_instance = DriftFluxModel(self)
            pass
        elif model.lower() == 'beggs' or model.lower() == 'beggsandbrill' or model.lower() == 'beggs&brill':
            self.model_instance = BeggsandBrillModel(self)
        else:
            raise Exception(f'{self.model} não é suportado!\nUse o modelo Homogêneo ou Drift-Flux')


    def prints(self):
        print('Densidade do óleo: ',self.rho_o)
        print('Densidade da água: ', self.rho_w)
        print('Densidade do gás: ', self.rho_g)
        print('Densidade da mistura: ', self.rho_m_ns)
        print('Viscosidade da mistura: ', self.mu_m_ns)
        print('Lambda L: ', self.lambda_L)
        print('Lambda G: ', self.lambda_G)


    def pipe_properties(self):
        if len(self.diameters) > 1:
            self.Dh = max(self.diameters) - min(self.diameters)
            self.Ap = np.pi/4*(max(self.diameters)**2 - min(self.diameters)**2)
        else:
            self.Dh = self.diameters[0]
            self.Ap = np.pi/4*self.diameters[0]**2

    def flow(self, force_update=False):
        # using SI units
        # for _ in range(2):
        if self.QL is not None and self.Bo is not None and self.Bw is not None and self.QOsc is not None:
            if self.QO is None:
                self.QO = self.QOsc * self.Bo
            if self.QW is None:
                self.QW = self.QWsc * self.Bw
        else:
            if self.QO is not None and self.QW is not None:
                self.QL = self.QO + self.QW
            elif self.QL is not None and self.QG is not None:
                pass
            else:
                # complicado, assumindo que deve ser calculado as vazoes na superficie
                if any([self.Bsw is None, self.Bo is None, self.Bw is None]) is True:
                    raise ValueError("Verifique se está sendo fornecido Bsw, Bo e Bw!")

                else:
                    if self.QLsc is not None:
                        if self.QOsc is None:
                            self.QOsc =self.QLsc * (1-self.Bsw)
                        if self.QWsc is None:
                            self.QWsc = self.QLsc * self.Bsw
                    else:
                        if self.QOsc is not None and self.QWsc is not None:
                            self.QLsc = self.QOsc + self.QWsc
                        elif self.QOsc is not None:
                            self.QLsc = self.QOsc/(1-self.Bsw)
                            self.QWsc = self.QLsc * self.Bsw
                        else:
                            if self.QWsc is not None:
                                self.QLsc = self.QWsc / self.Bsw
                                self.QOsc = self.QLsc * (1-self.Bsw)
                    if self.QOsc is not None and self.QWsc is not None:
                        self.QO = self.QOsc * self.Bo
                        self.QW = self.QWsc * self.Bw
                        self.QL = self.QO + self.QW

        if force_update is True:
            self.QO = self.QOsc * self.Bo
            self.QW = self.QWsc * self.Bw
            self.QL = self.QO + self.QW

        if self.Rs == 0:
            if self.Bsw == 0:
                self.QGsc = self.RGO * self.QOsc
                self.QG = self.Bg * self.QGsc

            elif self.Bsw > 0:
                self.QGsc = self.RGL * self.QLsc
                self.QG = (self.QGsc - self.QWsc * self.Rsw) * self.Bg

        else:
            if self.RGO is not None or self.RGL is not None and self.Rs is not None and self.Bg is not None:
                if self.Bsw == 0:
                    self.QGsc = self.RGO * self.QOsc
                    self.QG = (self.RGO - self.Rs) * self.Bg * self.QOsc

                else:
                    self.QGsc = self.RGL * self.QLsc
                    self.QG = (self.QGsc - self.QOsc * self.Rs - self.QWsc * self.Rsw) * self.Bg

        self.vsl = self.QL / self.Ap
        self.vsg = self.QG / self.Ap
        print('ql',self.QL)
        print('qg',self.QG)
        self.vm = self.vsl + self.vsg
        if self.Bo is not None and self.Bg is not None and self.rho_o is not None and self.rho_w is not None:
            print('bo',self.Bo)
            print('Bg',self.Bg)
            print('massica_o: ',self.rho_o * self.QO)
            print('massica_w: ',self.rho_w * self.QW)
            print('massica_g: ',self.rho_g * self.QG)
            print('massica_total: ',self.rho_o * self.QO+self.rho_w * self.QW+self.rho_g * self.QG)
            self.Qmassica.append(self.rho_o * self.QO+self.rho_w * self.QW+self.rho_g * self.QG)
            print('\n\n\n')

    def fractions(self):
        # Drift Flux nõa usa isso...
        self.lambda_L = self.QL / (self.QL + self.QG)
        self.lambda_G = self.QG / (self.QL + self.QG)

    def mixtureproperties(self, force_update=False):
        if force_update is False:
            if self.rho_l is None:
                self.rho_l = self.rho_o * self.QO/self.QL + self.rho_w * self.QW/self.QL
            if self.mu_l is None:
                self.mu_l = self.mu_o * self.QO/self.QL  + self.mu_w * self.QW/self.QL
            if self.sig_lg is None:
                self.sig_lg = self.sig_og * self.QO/self.QL + self.sig_wg * self.QW/self.QL
            self.rho_m_ns = self.lambda_L * self.rho_l + self.lambda_G * self.rho_g
            self.mu_m_ns = self.lambda_L * self.mu_l + self.lambda_G * self.mu_g
        elif force_update is True:
            self.rho_l = self.rho_o * self.QO / self.QL + self.rho_w * self.QW / self.QL
            self.mu_l = self.mu_o * self.QO / self.QL + self.mu_w * self.QW / self.QL
            self.sig_lg = self.sig_og * self.QO / self.QL + self.sig_wg * self.QW / self.QL
            self.rho_m_ns = self.lambda_L * self.rho_l + self.lambda_G * self.rho_g
            self.mu_m_ns = self.lambda_L * self.mu_l + self.lambda_G * self.mu_g

    def architecture(self):
        if self.length is not None:
            self.L_total = 0
            for L in self.length:
                self.L_total += L
            if self.N_division is None:
                self.dl = self.L_total/100
            else:
                self.dl = self.L_total/self.N_division
        else:
            self.dl = None


    def PostProcess(self, pressure, position, vm, temperature):
        if self.postprocess is True:
            print(self.Qmassica[0]- self.Qmassica[-1])
            print(100*(self.Qmassica[0] - self.Qmassica[-1])/self.Qmassica[-1])
            print(self.Qmassica[0])
            print(self.Qmassica[-1])
            import PostProcess

            PostProcess.run(pressure, position, vm, temperature)



class HomogeneousFlowModel:
    def __init__(self, initial_properties):
        self.initial_properties = initial_properties
        self.TPD = 0
        self.PDG = 0
        self.PDA = 0
        self.PDF = 0

    # def mixture_velocity(self):
    #     self.vm = (self.initial_properties.QL + self.initial_properties.QG) / self.initial_properties.Ap
    #     print(self.initial_properties.QL)
    #     print(self.initial_properties.QG)

    def HoldUp(self):
        self.H_L = self.initial_properties.lambda_L
        self.H_G = self.initial_properties.lambda_G

    def calculate_reynolds(self):
        self.Re = self.initial_properties.rho_m_ns * self.initial_properties.vm * self.initial_properties.Dh/self.initial_properties.mu_m_ns

    def friction_coefficient(self):
        self.f = 0.0055 * (1+(2*10**4*self.initial_properties.rugos/self.initial_properties.Dh + 1e6/self.Re)**(1/3))

    def pressure_drop_friction(self):
        self.dPDFdl = self.f * self.initial_properties.rho_m_ns * self.initial_properties.vm ** 2 / (2 * self.initial_properties.Dh)

    def pressure_drop_gravity(self, theta=0):
        self.dPDGdl = self.initial_properties.rho_m_ns * 9.81 * np.sin(theta*np.pi/180)

    def pressure_drop_acceleration_coefficent(self):
        self.mass_flow = self.initial_properties.rho_m_ns * self.initial_properties.vm * self.initial_properties.Ap
        self.mass_flow_g = self.initial_properties.rho_g * self.initial_properties.QG / self.initial_properties.Ap * self.initial_properties.Ap
        self.X = self.mass_flow_g/self.mass_flow
        if self.initial_properties.Z is None and self.initial_properties.Bg is not None:
            self.z_factor = self.initial_properties.Bg * 288.15/101325 * (self.initial_properties.Pressure[-1]/self.initial_properties.temperature[-1])
        else:
            self.z_factor = self.initial_properties.Z

        if self.z_factor is not None:
            self.Mg = self.z_factor * 8314 * self.initial_properties.temperature[-1] * self.initial_properties.rho_g / self.initial_properties.Pressure[-1]
        else:
            self.Mg = 28.97 * self.initial_properties.dg

        self.K = - (self.mass_flow/self.initial_properties.Ap)**2 * self.Mg * self.X / (8314 * self.initial_properties.temperature[-1] * self.initial_properties.rho_g**2)

    def total_pressure_drop(self, dl=None):
        self.dTPDdl = -(self.dPDFdl + self.dPDGdl)/(1+self.K)
        self.dPDAdl = self.dTPDdl * self.K
        if dl is not None:
             # ACUMULAVA ANTES MAS AGORA N MAISS!!!!!!!!!!!!!!
            self.TPD = -(self.dPDFdl + self.dPDGdl) / (1 + self.K) * dl
            self.PDA = self.dPDAdl * dl
            self.PDG = self.dPDGdl * dl
            self.PDF = self.dPDFdl * dl

    def output(self, incli=0, dl=1):
        self.run(incli=incli, dl=dl)
        if dl is None:
            print('Perda de carga por gravidade: ',self.dPDGdl,'Pa/m')
            print('Perda de carga por fricção: ', self.dPDFdl, 'Pa/m')
            print('Perda de carga por aceleração: ', self.dPDAdl, 'Pa/m')
            print('Perda de carga total: ', self.dTPDdl, 'Pa/m')
        else:
            return self.TPD

    def run(self,incli=0, dl=None):
        # self.mixture_velocity()
        self.HoldUp()
        self.calculate_reynolds()
        self.friction_coefficient()
        self.pressure_drop_gravity(incli)
        self.pressure_drop_friction()
        self.pressure_drop_acceleration_coefficent()
        self.total_pressure_drop(dl)

class DriftFluxModel:
    """
    Modelo Drift Flux de Zuber & Findlay (1965)
    Todas as velocidades são em m/s
    Ângulo em graus.
    """

    def __init__( self, initial_properties):
        self.initial_properties = initial_properties
        self.TPD = 0
        self.PDG = 0
        self.PDA = 0
        self.PDF = 0

    def froude(self):
        return abs(self.initial_properties.vm) / np.sqrt(9.81 * self.initial_properties.Dh)

    def C0(self, Fr, theta):
        if Fr < 3.5:
            return 1.05 + 0.15 * np.sin(theta)
        else:
            return 1.20

    def v_drift(self,theta, Fr):
        if Fr < 3.5:
            return np.sqrt(9.81*self.initial_properties.Dh)*(0.35*np.sin(theta)+(0.54*np.cos(theta)))

        else:
            return 0.35*np.sqrt(9.81*self.initial_properties.Dh)*np.sin(theta)

    def v_G(self,theta, C0, vd):
        return C0*self.initial_properties.vm + vd

    def HoldUp(self,theta, vg):
        self.H_G = self.initial_properties.vsg / vg
        self.H_L = 1-self.H_G

    def mixtureproperties(self):
        rho_m = self.H_G*self.initial_properties.rho_g + self.H_L*self.initial_properties.rho_l
        mu_m = self.H_G*self.initial_properties.mu_g + self.H_L*self.initial_properties.mu_l
        return rho_m, mu_m

    @staticmethod
    def HL(HG):
        return 1 - HG

    def dPdL_total(self, vm, HL, f, dpdl_a):
        rho_m = self.rho_m(HL)
        dpdl = ((-f *rho_m*(vm**2))/((2*self.Dh)) - (rho_m*self.g*np.sin(self.theta)) - (dpdl_a))

        return dpdl

    def calculate_reynolds(self, rho_m, mu_m):
        self.Re = rho_m * self.initial_properties.vm * self.initial_properties.Dh/mu_m


    def friction_coefficient(self):
        self.f = 0.0055 * (1+(2*10**4*self.initial_properties.rugos/self.initial_properties.Dh + 1e6/self.Re)**(1/3))

    def pressure_drop_friction(self, rho_m):
        self.dPDFdl = self.f * rho_m * self.initial_properties.vm ** 2 / (2 * self.initial_properties.Dh)

    def pressure_drop_gravity(self, rho_m, theta=0):
        self.dPDGdl = rho_m * 9.81 * np.sin(theta)

    def pressure_drop_acceleration_coefficent(self):
        self.mass_flow = self.initial_properties.rho_m_ns * self.initial_properties.vm * self.initial_properties.Ap
        self.mass_flow_g = self.initial_properties.rho_g * self.initial_properties.QG / self.initial_properties.Ap * self.initial_properties.Ap
        self.X = self.mass_flow_g/self.mass_flow
        if self.initial_properties.Z is None and self.initial_properties.Bg is not None:
            self.z_factor = self.initial_properties.Bg * 288.15/101325 * (self.initial_properties.Pressure[-1]/self.initial_properties.temperature[-1])
        else:
            self.z_factor = self.initial_properties.Z

        if self.z_factor is not None:
            self.Mg = self.z_factor * 8314 * self.initial_properties.temperature[-1] * self.initial_properties.rho_g / self.initial_properties.Pressure[-1]

        else:
            self.Mg = 28.97 * self.initial_properties.dg

        self.K = - (self.mass_flow/self.initial_properties.Ap)**2 * self.Mg * self.X / (8314 * self.initial_properties.temperature[-1] * self.initial_properties.rho_g**2)

    def total_pressure_drop(self, dl=None):
        self.dTPDdl = -(self.dPDFdl + self.dPDGdl)/(1+self.K)
        self.dPDAdl = self.dTPDdl * self.K
        if dl is not None:
             # ACUMULAVA ANTES MAS AGORA N MAISS!!!!!!!!!!!!!!
            self.TPD = -(self.dPDFdl + self.dPDGdl) / (1 + self.K) * dl
            self.PDA = self.dPDAdl * dl
            self.PDG = self.dPDGdl * dl
            self.PDF = self.dPDFdl * dl

    def output(self, incli, dl=None):
        theta = incli * np.pi / 180
        Fr = self.froude()
        C0 = self.C0(Fr, theta)
        vd = self.v_drift(theta, Fr)
        vg = self.v_G(theta, C0, vd)
        self.HoldUp(theta, vg)
        rho_m, mu_m = self.mixtureproperties()
        self.calculate_reynolds(rho_m, mu_m)
        self.friction_coefficient()
        self.pressure_drop_gravity(rho_m, theta)
        self.pressure_drop_friction(rho_m)
        self.pressure_drop_acceleration_coefficent()
        self.total_pressure_drop(dl)
        if dl is None:
            print('Perda de carga por gravidade: ',self.dPDGdl,'Pa/m')
            print('Perda de carga por fricção: ', self.dPDFdl, 'Pa/m')
            print('Perda de carga por aceleração: ', self.dPDAdl, 'Pa/m')
            print('Perda de carga total: ', self.dTPDdl, 'Pa/m')
        else:
            return self.TPD
        return {
            "Fr": Fr,
            "C0": C0,
            "vd": vd,
            "vg": vg,
            "HG": self.H_G,
            "rho_m": rho_m,
        }


class BeggsandBrillModel:
    """
    Modelo Beggs and Brill (1973)
    Todas as velocidades são em m/s
    Ângulo em graus.
    """
    def __init__(self, initial_properties):
        self.initial_properties = initial_properties
        self.TPD = 0
        self.PDG = 0
        self.PDA = 0
        self.PDF = 0

    def Froude2(self):
        self.Fr2 = self.initial_properties.vm**2 / (9.81 * self.initial_properties.Dh)

    def HoldUp(self):
        N_LV = self.initial_properties.vsl * (self.initial_properties.rho_l/(9.81 * self.initial_properties.sig_gl))
        if self.pattern == 'Distribuido':
            a, b, c = 1.065, 0.5824, 0.0609
        elif self.pattern == 'Segregado':
            a, b, c = 0.98, 0.4846, 0.0868
            # d, e, f, g =
        elif self.pattern == 'Intermitente':
            a, b, c = 0.845, 0.5351, 0.0173

        if self.pattern == 'Distribuido':
            H_LO = a*self.initial_properties.lambda_L**b/(self.Fr2**c)
            psi = 1
        # hold up ................ FAZERRRRRRRRRRRRR
        self.H_L = H_LO * psi

    def mixtureproperties(self):
        self.rho_m = self.H_G*self.initial_properties.rho_g + self.H_L*self.initial_properties.rho_l

    def calculate_reynolds(self):
        self.Re = self.initial_properties.rho_m_ns * self.initial_properties.vm * self.initial_properties.Dh/self.initial_properties.mu_m_ns

    def friction_coefficient(self):
        fn = 0.0055 * (1+(2*10**4*self.initial_properties.rugos/self.initial_properties.Dh + 1e6/self.Re)**(1/3))
        y = self.initial_properties.lambda_L/self.H_L**2
        if 1.0 < y < 1.2:
            s = math.log(2.2*y -1.2)
        else:
            a = -0.0523 + 3.182*math.log(y) - 0.8725*math.log(y)**2 + 0.001855* math.log(y)**4
            s= math.log(y)/a
        self.f =fn * np.e**s

    def pressure_drop_friction(self):
        self.dPDFdl = self.f * self.initial_properties.rho_m_ns * self.initial_properties.vm ** 2 / (2 * self.initial_properties.Dh)

    def pressure_drop_gravity(self, theta=0):
        self.dPDGdl = self.rho_m * 9.81 * np.sin(theta*np.pi/180)

    def pressure_drop_acceleration_coefficent(self):
        self.K = -(self.rho_m * self.initial_properties.vm * self.initial_properties.vsg)/self.initial_properties.Pressure[-1]

    def total_pressure_drop(self, dl=None):
        self.dTPDdl = -(self.dPDFdl + self.dPDGdl)/(1+self.K)
        self.dPDAdl = self.dTPDdl * self.K
        if dl is not None:
            self.TPD = -(self.dPDFdl + self.dPDGdl) / (1 + self.K) * dl
            self.PDA = self.dPDAdl * dl
            self.PDG = self.dPDGdl * dl
            self.PDF = self.dPDFdl * dl

    def output(self, incli=0, dl=1):
        self.run(incli=incli, dl=dl)
        if dl is None:
            print('Perda de carga por gravidade: ',self.dPDGdl,'Pa/m')
            print('Perda de carga por fricção: ', self.dPDFdl, 'Pa/m')
            print('Perda de carga por aceleração: ', self.dPDAdl, 'Pa/m')
            print('Perda de carga total: ', self.dTPDdl, 'Pa/m')
        else:
            return self.TPD

    def FlowPattern(self):
        L1 = 316 * self.initial_properties.lambda_L ** 0.302
        L2 = 0.0009252 * self.initial_properties.lambda_L ** -2.4684
        L3 = 0.1 * self.initial_properties.lambda_L ** -1.4516
        L4 = 0.5 * self.initial_properties.lambda_L ** -6.738
        if self.initial_properties.lambda_L < 0.4 and self.Fr2 >= L1 or self.initial_properties.lambda_L >= 0.4 and self.Fr2 > L4:
            self.pattern = 'Distribuido'
        elif self.initial_properties.lambda_L < 0.001 and self.Fr2 < L1 or self.initial_properties.lambda_L >= 0.001 and self.Fr2 < L2:
            self.pattern = 'Segregado'
        elif self.initial_properties.lambda_L >= 0.01 and L2 <= self.Fr2 <= L3:
            self.pattern = 'Transição'
        elif 0.01 <= self.initial_properties.lambda_L <= 0.4 and L3 <= self.Fr2 <= L1 or self.initial_properties.lambda_L >= 0.4 and L3 <= self.Fr2 <= L4:
            self.pattern = 'Intermitente'

    def run(self,incli=0, dl=None):
        self.Froude2()
        self.FlowPattern()
        self.HoldUp()
        self.mixtureproperties()
        self.calculate_reynolds()
        self.friction_coefficient()
        self.pressure_drop_gravity(incli)
        self.pressure_drop_friction()
        self.pressure_drop_acceleration_coefficent()
        self.total_pressure_drop(dl)

