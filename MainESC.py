from BlackOilModel import GasPhase
from BlackOilModel import OilPhase
from BlackOilModel import WaterPhase
# from backup.MultiPhaseFlowModels import MultiPhaseFlowModel
from MultiPhaseFlowModels import MultiPhaseFlowModel
import unitsconverter
import math
# Dados iniciais
P_sep = 10 # 10 # bar
P_res = 450 # bar
L_riser = 1050 # m
L_solomarinho = 700 # m
L_poco = 850/math.sin(60*math.pi/180)
inclinacoes = [60,0,90]
# inclinacoes = [60*math.pi/180,0,90*math.pi/180]
T_poco_inicial = 55 # 80 # °C
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
T = unitsconverter.Temperature(T_riser_final,'C','K')
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

def MultiPhaseFlowModeloptimize(Dh):
    try:
        Modelo = MultiPhaseFlowModel(model='homogeneo', temperature=T, pressure=P, dg=dg, salinity=salinidade, do=do,
                            diameters=[unitsconverter.Length(Dh, 'in', 'm')], rugosity=rugos_temp, incli=inclinacoes, rho_o=rho_o,
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


from scipy.optimize import minimize
# from scipy.optimize import root
# from scipy.optimize import root_scalar
# D_ideal  = root_scalar(
#     MultiPhaseFlowModeloptimize,
#     bracket=[3.0,3.5],
#     method='brentq',
#     xtol=1e-4
# )

D_ideal = minimize(MultiPhaseFlowModeloptimize,x0=[4], bounds=[[0.5,7]],method='Nelder-Mead')
# 2.2341
# print(D_ideal.x)
# D_ideal = root(MultiPhaseFlowModeloptimize,x0=7,method='hybr')
# D_ideal_novo = 3.9 # pol
D_ideal = [unitsconverter.Length(D_ideal.x , 'in', 'm')]
Modelo = MultiPhaseFlowModel(model='homogeneo', temperature=T, pressure=P, dg=dg, salinity=salinidade, do=do, diameters=D_ideal,rugosity=rugos_temp, incli=inclinacoes, rho_o=rho_o,
                             rho_g=rho_g, rho_w=rho_w, QLsc=Q,Bo=Bo, Bw=Bw, Bg=Bg, Bsw=BSW, RGL=RGL, Rs=Rs, Rsw=Rsw,
                             mu_o=mu_o, mu_w=mu_w, mu_g=mu_g, sig_og=sig_og, sig_wg=sig_wg, Z=Gas.Z, length=[L_poco,L_solomarinho,L_riser], N_division= 100, reference_pressure=P_ref, postprocess=True)
Modelo.run()

"""
Massa específica da água: 62.3899 lbm/ft³
Compressibilidade da água: 2.6816e-06 psi⁻¹
Viscosidade da água: 1.1381 cp
"""