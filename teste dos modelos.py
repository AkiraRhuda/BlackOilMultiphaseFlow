from MultiPhaseFlowModels import MultiPhaseFlowModel
import unitsconverter

#################################################################################
################################### Homogeneo ###################################
#################################################################################

# Dados iniciais
P_sep = unitsconverter.Pressure(5, 'bar', 'Pa') # bar
T = unitsconverter.Temperature(55, 'C', 'K') # ° C
P_b = unitsconverter.Pressure(250, 'bar', 'Pa')  # bar
L = 1200 # m
diam = [unitsconverter.Length(7, 'in', 'm'), unitsconverter.Length(3, 'in', 'm')]
rugos = unitsconverter.Length(0.0075, 'in', 'm')

QOsc = 0.05787 # m3/d
Bsw = 30 # %
Rgl = 150 # sm3 / sm3
Bw = 1.05 # m3 / sm3
Bo = 1.18 # m3 / sm3
Bg = 0.008 # m3 / sm3
Rs = 50 # sm3 / sm3
Rsw = 10 # sm3 / s
sig_og = 0.00841# N/m
sig_wg = 0.004 # N/m
rho_w = 1032 # km/m3
rho_o = 800 # kg/m3
rho_g = 85 # km/m3
mu_w = unitsconverter.Viscosity(1.2, 'cp', 'Pa.s') # cp
mu_o = unitsconverter.Viscosity(0.8, 'cp', 'Pa.s') # cp
mu_g = unitsconverter.Viscosity(0.016, 'cp', 'Pa.s') # cp


Modelo = MultiPhaseFlowModel(model='homogeneo', temperature=T, pressure=P_b, diameters=diam,rugosity=rugos, incli=[90], rho_o=rho_o, rho_g=rho_g, rho_w=rho_w,
                    do=unitsconverter.Density(rho_o, 'kg/m3', 'do'), dg=unitsconverter.Density(rho_g, 'kg/m3', 'dg') ,  salinity=0, QOsc=QOsc,Bo=Bo, Bw=Bw, Bg=Bg, Bsw=Bsw, RGL=Rgl, Rs=Rs, Rsw=Rsw, mu_o=mu_o, mu_w=mu_w, mu_g=mu_g, sig_og=sig_og, sig_wg=sig_wg)
Modelo.run()
# print(Modelo.model_instance.dTPDdl)

#################################################################################
#################################################################################
#################################################################################

"""from Tools.scripts.fixdiv import multi_ok

from BlackOilModel import GasPhase
from BlackOilModel import OilPhase
from BlackOilModel import WaterPhase
from MultiPhaseFlowModels import MultiPhaseFlowModel
import unitsconverter
import math
# Dados iniciais
P_sep = 250 # 10 # bar
L_riser = 1050 # m
L_solomarinho = 700 # m
L_poco = 850/math.sin(60*math.pi/180)

T_poco_inicial = 55 # 80 # °C
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

Gas = GasPhase(zcorrelation="Papay", μcorrelation='Lee', dg=dg, P=unitsconverter.Pressure(P_sep,'bar','psi'), T=unitsconverter.Temperature(T_riser_final, 'C','R')).output()
Oil = OilPhase("Standing", "Standing", "Beggs",do=do, dg=dg, P=unitsconverter.Pressure(P_sep,'bar','psi'), T=unitsconverter.Temperature(T_riser_final,'C','F'), RGL=RGL, BSW=BSW).output()
Water = WaterPhase(P=unitsconverter.Pressure(P_sep,'bar','psi'), T=unitsconverter.Temperature(T_riser_final,'C','F'), S=salinidade).output()

diam_temp = [unitsconverter.Length(7, 'in', 'm'), unitsconverter.Length(3, 'in', 'm')]
rugos_temp = unitsconverter.Length(0.0075, 'in', 'm')
# conversoes
T = unitsconverter.Temperature(T_riser_final,'C','K')
P = unitsconverter.Pressure(P_sep,'bar','Pa')
rho_o = unitsconverter.Density(Oil.ρo, 'lbm/ft3', 'kg/m3')
rho_g = unitsconverter.Density(Gas.ρ, 'lbm/ft3', 'kg/m3')
rho_w = unitsconverter.Density(Water.ρw, 'lbm/ft3', 'kg/m3')
mu_o = unitsconverter.Viscosity(Oil.μob, 'cp', 'Pa.s')
mu_g = unitsconverter.Viscosity(Gas.μ, 'cp', 'Pa.s')
mu_w = unitsconverter.Viscosity(Water.μw, 'cp', 'Pa.s')
QoSc = 0.05787
Bo = unitsconverter.Density(Oil.Bo, 'scfstb', 'sm3sm3')
Bw = unitsconverter.Density(Water.Bw, 'scfstb', 'sm3sm3')
Bg = unitsconverter.Density(Gas.Bg(), 'scfstb', 'sm3sm3')
Rs = unitsconverter.Density(Oil.Rs, 'scfstb', 'sm3sm3')
Rsw = unitsconverter.Density(Water.Rsw, 'scfstb', 'sm3sm3')

Modelo = MultiPhaseFlowModel(model='homogeneo', temperature=T, pressure=P, diameters=diam_temp,rugosity=rugos_temp, incli=90, rho_o=rho_o,
                            rho_g=rho_g, rho_w=rho_w, QOsc=QoSc,Bo=Bo, Bw=Bw, Bg=Bg, Bsw=BSW, RGL=RGL, Rs=Rs, Rsw=Rsw,
                            mu_o=mu_o, mu_w=mu_w, mu_g=mu_g, sig_og=sig_og, sig_wg=sig_wg, Z=Gas.Z)

"""



"""
Massa específica da água: 62.3899 lbm/ft³
Compressibilidade da água: 2.6816e-06 psi⁻¹
Viscosidade da água: 1.1381 cp
"""

##################################################################################
################################### DRIFT FLUX ###################################
##################################################################################

# Dados iniciais
P_sep = unitsconverter.Pressure(5, 'bar', 'Pa') # bar
T = unitsconverter.Temperature(25, 'C', 'K') # ° C
P_b = unitsconverter.Pressure(250, 'bar', 'Pa')  # bar
L = 1200 # m
diam = [5.1e-2]
rugos = 0.045e-3
m_l = 1.5
import numpy as np
 # m3/d
Bsw = 30 # %
Rgl = 150 # sm3 / sm3
Bw = 1.05 # m3 / sm3
Bo = 1.18 # m3 / sm3
Bg = 0.008 # m3 / sm3
Rs = 50 # sm3 / sm3
Rsw = 10 # sm3 / s
sig_lg = 0.00841# N/m
rho_l = 850
rho_g = 1.5 # km/m3
dg = 0.7
QLsc = m_l/rho_l
QG = 0.015/rho_g
mu_l = unitsconverter.Viscosity(2, 'cp', 'Pa.s') # cp
mu_g = unitsconverter.Viscosity(0.02, 'cp', 'Pa.s') # cp


Modelo = MultiPhaseFlowModel(model='driftflux', temperature=T, pressure=P_b, diameters=diam,rugosity=rugos, incli=[10], rho_l=rho_l, rho_g=rho_g,
                    dg=dg, salinity=0, QL=QLsc, QG=QG,mu_l=mu_l, mu_g=mu_g, sig_lg=sig_lg)
# Modelo = MultiPhaseFlowModel(model='driftflux', temperature=T, pressure=P_b, diameters=diam,rugosity=rugos, incli=[10], rho_l=rho_l, rho_g=rho_g,
#                     do=unitsconverter.Density(rho_o, 'kg/m3', 'do'), dg=unitsconverter.Density(rho_g, 'kg/m3', 'dg') ,  salinity=0, QOsc=QOsc,Bo=Bo, Bw=Bw, Bg=Bg, Bsw=Bsw, RGL=Rgl, Rs=Rs, Rsw=Rsw, mu_o=mu_o, mu_w=mu_w, mu_g=mu_g, sig_og=sig_og, sig_wg=sig_wg)
#
#Modelo.run()
# print(Modelo.model_instance.dTPDdl)

#######################################################################################
################################### Beggs and Brill ###################################
#######################################################################################




# Dados iniciais
P_sep = unitsconverter.Pressure(5, 'bar', 'Pa') # bar
T = unitsconverter.Temperature(75, 'C', 'K') # ° C
P_b = unitsconverter.Pressure(150, 'bar', 'Pa')  # bar
L = 1500 # m
diam = [0.1016] # m
rugos = 0.005 * 0.0256
QLsc = 4000 # m3/d
Bsw = 20 # %
RGO = 140 # sm3 / sm3
Bw = 1.55 # m3 / sm3
Bo = 1.15 # m3 / sm3
Bg = 0.0095 # m3 / sm3
Rs = 50 # sm3 / sm3
Rsw = 10 # sm3 / s
sig_og = 0.00841 # N/m
sig_wg = 0.004 # N/m
rho_o = 920 # km/m3
rho_w = 981.4 # km/m3
rho_g = 85 # km/m3
QLsc = QLsc/86400
mu_o = unitsconverter.Viscosity(0.8, 'cp', 'Pa.s') # cp
mu_w = unitsconverter.Viscosity(0.381, 'cp', 'Pa.s') # cp
mu_g = unitsconverter.Viscosity(0.016, 'cp', 'Pa.s') # cp


Modelo = MultiPhaseFlowModel(model='beggsandbrill', temperature=T, pressure=P_b, diameters=diam,rugosity=rugos, incli=[45], rho_o=rho_o, rho_w=rho_w,
                             rho_g=rho_g, salinity=0, QLsc=QLsc,mu_o=mu_o, mu_w=mu_w, mu_g=mu_g, sig_og=sig_og, sig_wg=sig_wg,
                             Rs=Rs, Rsw=Rsw, Bg=Bg, Bo=Bo, Bw=Bw, RGO=RGO, Bsw=Bsw)
# Modelo = MultiPhaseFlowModel(model='driftflux', temperature=T, pressure=P_b, diameters=diam,rugosity=rugos, incli=[10], rho_l=rho_l, rho_g=rho_g,
#                     do=unitsconverter.Density(rho_o, 'kg/m3', 'do'), dg=unitsconverter.Density(rho_g, 'kg/m3', 'dg') ,  salinity=0, QOsc=QOsc,Bo=Bo, Bw=Bw, Bg=Bg, Bsw=Bsw, RGL=Rgl, Rs=Rs, Rsw=Rsw, mu_o=mu_o, mu_w=mu_w, mu_g=mu_g, sig_og=sig_og, sig_wg=sig_wg)
#
Modelo.run()
print(Modelo.model_instance.dTPDdl)