from BlackOilModel import GasPhase
from BlackOilModel import OilPhase
from BlackOilModel import WaterPhase
import unitsconverter
import math
# Dados iniciais
P_sep = 10 # bar
L_riser = 1050 # m
L_solomarinho = 700 # m
L_poco = 850/math.sin(60*math.pi/180)

T_poco_inicial = 80 # °C
T_solomarinho = 4 # °C
T_riser_final = 15 # °C

BSW = 0.25 # %
RGL = 800 # sm3/sm3
API = 30 # °
dg = 0.72
P_res = 450 # bar
sig_og = 0.00841 # N/m
sig_wg = 0.030 # N/m
TEC_poco = 2 # W/m.K
TEC_marinho = 1 # W/m.K

Q = 3000 # sm3/d
# Calculos
do = unitsconverter.Density(API, 'API', 'd')

# Dados estimados
salinidade = 0.15

GasPhase(zcorrelation="Papay", μcorrelation='Lee', dg=dg, P=unitsconverter.Pressure(P_sep,'bar','psi'), T=unitsconverter.Temperature(T_riser_final, 'C','R')).output()
OilPhase("Standing", "Standing", "Beggs",do=do, dg=dg, P=unitsconverter.Pressure(P_sep,'bar','psi'), T=unitsconverter.Temperature(T_riser_final,'C','F'), RGL=RGL, BSW=BSW).output()
WaterPhase(P=unitsconverter.Pressure(P_sep,'bar','psi'), T=unitsconverter.Temperature(T_riser_final,'C','F'), S=salinidade).output()

"""
Massa específica da água: 62.3899 lbm/ft³
Compressibilidade da água: 2.6816e-06 psi⁻¹
Viscosidade da água: 1.1381 cp
"""