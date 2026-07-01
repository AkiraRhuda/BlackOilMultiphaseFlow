from ctypes.wintypes import PBOOL

from BlackOilModel import GasPhase
from BlackOilModel import OilPhase
from BlackOilModel import WaterPhase
# from OLD import MultiPhaseFlowModel
from MultiPhaseFlowModels import MultiPhaseFlowModel
import unitsconverter
import math
import numpy as np

# Dados iniciais
P_sep = 10 # 10 # bar
P_res = 450 # bar
L_riser = 1050 # m
L_solomarinho = 700 # m
L_poco = 850/math.sin(60*math.pi/180)
inclinacoes = [60,0,90]
# inclinacoes = [60*math.pi/180,0,90*math.pi/180]
T_poco_inicial = 80 # °C
T_solomarinho = 4 # °C
T_riser_final = 15 # °C

BSW = 25 # %
RGL = 800 # sm3/sm3
API = 30 # °
dg = 0.72
sig_og = 0.00841 # N/m
sig_wg = 0.030 # N/m
TEC_poco = 2 # W/m.K
TEC_marinho = 1 # W/m.K

Q = 3000 # sm3/d

# Calculos
do = unitsconverter.Density(API, 'API', 'do')

# Dados estimados
salinidade = 0.15

Gas = GasPhase(zcorrelation="Papay", μcorrelation='Lee', dg=dg, P=unitsconverter.Pressure(P_res,'bar','psi'), T=unitsconverter.Temperature(T_riser_final, 'C','R')).output()
Oil = OilPhase("Standing", "Standing", "Beggs",do=do, dg=dg, P=unitsconverter.Pressure(P_res,'bar','psi'), T=unitsconverter.Temperature(T_riser_final,'C','F'), RGL=RGL, BSW=BSW).output()
Water = WaterPhase(P=unitsconverter.Pressure(P_res,'bar','psi'), T=unitsconverter.Temperature(T_riser_final,'C','F'), S=salinidade).output()
Dh_ideal = 3.376504 # gabarito
diam_temp = [unitsconverter.Length(Dh_ideal, 'in', 'm')]
rugos_temp = unitsconverter.Length(0.0075, 'in', 'm')

# conversoes
T = T_poco_inicial+273.15
P = unitsconverter.Pressure(P_res,'bar','Pa')
P_ref = unitsconverter.Pressure(P_sep, 'bar','Pa')
rho_o = unitsconverter.Density(Oil.ρo, 'lbm/ft3', 'kg/m3')
rho_g = unitsconverter.Density(Gas.ρ, 'lbm/ft3', 'kg/m3')
rho_w = unitsconverter.Density(Water.ρw, 'lbm/ft3', 'kg/m3')
mu_o = unitsconverter.Viscosity(Oil.μob, 'cp', 'Pa.s')
mu_g = unitsconverter.Viscosity(Gas.μ, 'cp', 'Pa.s')
mu_w = unitsconverter.Viscosity(Water.μw, 'cp', 'Pa.s')
Q = Q/86400
Bo = unitsconverter.Standard(Oil.Bo, 'scfstb', 'sm3sm3')
Bw = unitsconverter.Standard(Water.Bw, 'scfstb', 'sm3sm3') # precisa ser sm3/sm3 .;........
Bg = unitsconverter.Standard(Gas.Bg(), 'scfstb', 'sm3sm3')
Rs = unitsconverter.Standard(Oil.Rs, 'scfstb', 'sm3sm3')
Rsw = unitsconverter.Standard(Water.Rsw, 'scfstb', 'sm3sm3')
Pb = unitsconverter.Pressure(Oil.Pb, 'psi', 'Pa')
# Se P > Pb é monofásico.

def MultiPhaseFlowModeloptimize(Dh):
    try:
        Dh = Dh[0]
        Dh_list = np.ones(3) * unitsconverter.Length(Dh, 'in', 'm')
        Modelo = MultiPhaseFlowModel(model='beggs', temperature=T, pressure=P, bubblepressure=Pb, dg=dg, salinity=salinidade, do=do,
                            diameters=Dh_list, rugosity=rugos_temp, incli=inclinacoes, rho_o=rho_o,
                            rho_g=rho_g, rho_w=rho_w, QLsc=Q, Bo=Bo, Bw=Bw, Bg=Bg, Bsw=BSW, RGL=RGL, Rs=Rs, Rsw=Rsw,
                            mu_o=mu_o, mu_w=mu_w, mu_g=mu_g, sig_og=sig_og, sig_wg=sig_wg, Z=Gas.Z,
                            length=[L_poco, L_solomarinho, L_riser], reference_pressure=P_ref)
        Modelo.run()
        print(Dh)
        if isinstance(Modelo.abserro, float) is True:
            return abs(Modelo.abserro)
        else:
            print(abs(Modelo.abserro))
            return 1e5
    except:
        return 1e10

#### DESCOMENTAR SE QUISER ACHAR O DIAMETRO IDEAL!!!!

# from scipy.optimize import minimize
# D_ideal = minimize(MultiPhaseFlowModeloptimize,x0=np.array(4).astype(float), bounds=[[0.5,7]],method='Nelder-Mead')

#################################################################################
################################### Homogeneo ###################################
#################################################################################

# D bizarro 1.700105911810369
# novo D ideal?? 2.1936011396849584
# D ideal para homogeneo = 2.196582005915114 pol
D_ideal = [unitsconverter.Length(1.700105911810369, 'in', 'm')] # pol
D_ideal = np.ones(3) * D_ideal
Modelo = MultiPhaseFlowModel(model='homogeneo', temperature=T, pressure=P, bubblepressure=Pb, dg=dg, salinity=salinidade, do=do, diameters=D_ideal,rugosity=rugos_temp, incli=inclinacoes, rho_o=rho_o,
                             rho_g=rho_g, rho_w=rho_w, QLsc=Q,Bo=Bo, Bw=Bw, Bg=Bg, Bsw=BSW, RGL=RGL, Rs=Rs, Rsw=Rsw,
                             mu_o=mu_o, mu_w=mu_w, mu_g=mu_g, sig_og=sig_og, sig_wg=sig_wg, Z=Gas.Z, length=[L_poco,L_solomarinho,L_riser], N_division= 100, reference_pressure=P_ref, postprocess=True)
Modelo.run()


##################################################################################
################################### Drift-flux ###################################
##################################################################################

# novo D ideal?? 2.313276380648494
# D ideal para drift-flux = 2.3038571715354914 pol
D_ideal = [unitsconverter.Length(2.313276380648494, 'in', 'm')]
D_ideal = np.ones(3) * D_ideal
Modelo = MultiPhaseFlowModel(model='drift-flux', temperature=T, pressure=P, bubblepressure=Pb, dg=dg, salinity=salinidade, do=do, diameters=D_ideal,rugosity=rugos_temp, incli=inclinacoes, rho_o=rho_o,
                             rho_g=rho_g, rho_w=rho_w, QLsc=Q,Bo=Bo, Bw=Bw, Bg=Bg, Bsw=BSW, RGL=RGL, Rs=Rs, Rsw=Rsw,
                             mu_o=mu_o, mu_w=mu_w, mu_g=mu_g, sig_og=sig_og, sig_wg=sig_wg, Z=Gas.Z, length=[L_poco,L_solomarinho,L_riser], N_division= 100, reference_pressure=P_ref, postprocess=True)
# Modelo.run()



#######################################################################################
################################### Beggs and Brill ###################################
#######################################################################################

# 2.3865680840648547 novo D ideal??
# diametro ótimo beggs 2.382999665 pol
D_nao_ideal = [unitsconverter.Length(3.0, 'in', 'm')]
D_nao_ideal = np.ones(3) * D_nao_ideal
Modelo1 = MultiPhaseFlowModel(model='beggs and brill', temperature=T, pressure=P, bubblepressure=Pb, dg=dg, salinity=salinidade, do=do, diameters=D_nao_ideal,rugosity=rugos_temp, incli=inclinacoes, rho_o=rho_o,
                             rho_g=rho_g, rho_w=rho_w, QLsc=Q,Bo=Bo, Bw=Bw, Bg=Bg, Bsw=BSW, RGL=RGL, Rs=Rs, Rsw=Rsw,
                             mu_o=mu_o, mu_w=mu_w, mu_g=mu_g, sig_og=sig_og, sig_wg=sig_wg, Z=Gas.Z, length=[L_poco,L_solomarinho,L_riser], N_division= 100, reference_pressure=P_ref, postprocess=True, direction='ascendente')
Modelo1.run()
 # 1.9271605419544349
D_ideal = [unitsconverter.Length(2.3865680840648547, 'in', 'm')]
D_ideal = np.ones(3) * D_ideal
Modelo2 = MultiPhaseFlowModel(model='beggs and brill', temperature=T, pressure=P, bubblepressure=Pb, dg=dg, salinity=salinidade, do=do, diameters=D_ideal,rugosity=rugos_temp, incli=inclinacoes, rho_o=rho_o,
                             rho_g=rho_g, rho_w=rho_w, QLsc=Q,Bo=Bo, Bw=Bw, Bg=Bg, Bsw=BSW, RGL=RGL, Rs=Rs, Rsw=Rsw,
                             mu_o=mu_o, mu_w=mu_w, mu_g=mu_g, sig_og=sig_og, sig_wg=sig_wg, Z=Gas.Z, length=[L_poco,L_solomarinho,L_riser], N_division= 100, reference_pressure=P_ref, postprocess=True, direction='ascendente')
Modelo2.run()

import PostProcess
PostProcess.run([Modelo1.Pressure,Modelo2.Pressure], Modelo2.Pos)