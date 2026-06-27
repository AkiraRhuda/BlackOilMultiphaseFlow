import numpy as np
import matplotlib.pyplot as plt
#Eng tools - valores experimentais
T_CH4 = np.array([200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850,900, 950, 1000, 1050, 1100 ]) # K
Cp_CH4 = np.array([2.087, 2.121, 2.156,2.191, 2.226, 2.293, 2.365, 2.442, 2.525, 2.703, 2.889, 3.074, 3.256, 3.432, 3.602, 3.766, 3.923, 4.072, 4.214, 4.348, 4.475, 4.595, 4.708]) # kJ/(kg·K)
T_C2H6 = np.array([250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900 ])
Cp_C2H6 = np.array([1.535, 1.651, 1.766, 1.878, 1.987, 2.095, 2.199, 2.402, 2.596, 2.782, 2.958, 3.126, 3.286, 3.438, 3.581, 3.717, 3.846])

curva_CH4 = np.poly1d(np.polyfit(T_CH4, Cp_CH4, 3))
curva_C2H6 = np.poly1d(np.polyfit(T_C2H6, Cp_C2H6, 3))

# print(curva_CH4)
# print(curva_C2H6)

T_fit_CH4 = np.linspace(200, 1100, 400)
T_fit_C2H6 = np.linspace(250, 900, 400)

"""plt.figure(figsize=(10,6))
plt.scatter(T_CH4, Cp_CH4, label = f'Metano ($CH_4$)')
plt.scatter(T_C2H6, Cp_C2H6, label = f'Etano ($C_2H_6$)')
plt.plot(T_fit_CH4, curva_CH4(T_fit_CH4), color='blue', linestyle='-', label='CH4 (ajuste)')
plt.plot(T_fit_C2H6, curva_C2H6(T_fit_C2H6), color='green', linestyle='-', label='C2H6 (ajuste)')
plt.title(f'$C_p(T)$')
plt.xlabel(f'Temperatura [K]')
plt.ylabel(r'$C_p\,[\mathrm{kJ\,kg^{-1}\,K^{-1}}]$')
plt.grid(True)
plt.legend()
plt.show()"""

from iapws import IAPWS97

class Temperature:
      def __init__(self, T, P, rho_o, QO, QW, QL, Mg, HL):
          # T em K, P em Pa, rho em kg/m3, Mg em Kg/kmol
          Mg = Mg/1000 # Kg/mol
          self.T = T
          self.P = P /1e6 # MPa
          CpW = self.Calculate_CpW() # KJ/(kg·K)
          CpW *= 1000 # J/(kg·K)
          CpO = self.Calculate_CpO(T, rho_o) # KJ/(kg·K)
          CpO *= 1000 # J/(kg·K)
          CpL = self.Calculate_CpL(CpO, CpW, QO, QW, QL) # KJ/(kg·K)
          CpG = self.Calculate_CpG(Mg, T) # KJ/(kg·K)
          CpG *= 1000 # J/(kg·K)
          self.CpM = self.Calculate_CpM(CpL, CpG, HL)

          # T = self.Calculate_Temperature(dl, depth, T, m_dot_m, CpM, theta)
          print(f"Cp água   = {CpW:.3f} [J/(kg·K)]")
          print(f"Cp óleo   = {CpO:.3f} [J/(kg·K)]")
          print(f"Cp líquido= {CpL:.3f} [J/(kg·K)]")
          print(f"Cp gás    = {CpG:.3f} [J/(kg·K)]")
          print(f"Cp mistura= {self.CpM:.3f} [J/(kg·K)]")
          # print(f"Temp da mistura= {T:.3f} [K]")
      def Calculate_CpG(self, Mg, T):
          # Valores experimentais
          M_CH4 = 0.01604  # kg/mol
          M_C2H6 = 0.03007 # kg/mol
          T_CH4 = np.array([200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850,900, 950, 1000, 1050, 1100 ]) # K
          Cp_CH4 = np.array([2.087, 2.121, 2.156,2.191, 2.226, 2.293, 2.365, 2.442, 2.525, 2.703, 2.889, 3.074, 3.256, 3.432, 3.602, 3.766, 3.923, 4.072, 4.214, 4.348, 4.475, 4.595, 4.708]) # kJ/(kg·K)
          T_C2H6 = np.array([250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900 ])
          Cp_C2H6 = np.array([1.535, 1.651, 1.766, 1.878, 1.987, 2.095, 2.199, 2.402, 2.596, 2.782, 2.958, 3.126, 3.286, 3.438, 3.581, 3.717, 3.846])

          curva_CH4 = np.poly1d(np.polyfit(T_CH4, Cp_CH4, 3))
          curva_C2H6 = np.poly1d(np.polyfit(T_C2H6, Cp_C2H6, 3))

          y_CH4 = (M_C2H6 - Mg) / (M_C2H6 - M_CH4)
          y_C2H6 = 1 - y_CH4

          Cp_CH4 = curva_CH4(T)  # kJ/(kg·K)
          Cp_C2H6 = curva_C2H6(T) # kJ/(kg·K)

          CpG = y_CH4 * Cp_CH4 + y_C2H6 * Cp_C2H6 # kJ/(kg·K)

          return CpG

      def Calculate_CpW(self):
          state = IAPWS97(T=self.T, P=self.P)
          CpW = state.cp
          return CpW # kJ/(kg·K)

      def Calculate_CpO(self,T, rho_o): # Correlação
          rho_w = 1000
          do = rho_o/rho_w
          CpO = (2e-3 * (T - 273.15) - 1.429) * do + (2.67e-3 * (T - 273.15)) + 3.049 # T em K
          return CpO  # kJ/(kg·K)

      def Calculate_CpL(self, CpO, CpW, QO, QW, QL):
          CpL = CpO * (QO/QL)+ CpW * (QW/QL)
          return CpL

      def Calculate_CpM(self,CpL, CpG, Holdup_Liquido):
          CpM = Holdup_Liquido * CpL + (1 - Holdup_Liquido) * CpG
          return CpM

      def Calculate_Temperature(self, dl, depth, m_dot_m, theta):
          theta = theta * np.pi / 180
          # T inf e TEC dependem de qual trecho esta, entao tem que saber a profundidade
          depth = np.asarray(depth, dtype=float)
          g = 9.81
          Tinf = np.zeros_like(depth)
          TEC = np.zeros_like(depth)

          trecho_riser = (depth >= 0) & (depth < 1050)
          Tinf[trecho_riser] = (15+273.15) - (11 / 1050) * depth[trecho_riser]
          TEC[trecho_riser] = 1

          trecho_leito = depth == 1050
          Tinf[trecho_leito] = 4+273.15
          TEC[trecho_leito] = 1

          trecho_poco = (depth > 1050) & (depth <= 1900)
          Tinf[trecho_poco] = 4 + 273.15 + (76 / 850) * (depth[trecho_poco] - 1050)
          TEC[trecho_poco] = 2 #

          T = Tinf - (m_dot_m * g * np.sin(theta)) / TEC - ((Tinf - (m_dot_m * g * np.sin(theta)) / TEC - self.T) * np.exp(-(TEC * dl) / (m_dot_m * self.CpM)))

          return T



"""T = 1000         # K
P = 0.101325       # MPa
rho_o = 850 #kg/m3
Mg = 0.0200        # kg/mol
QO = 0.020
QW = 0.010
QL = QO + QW
Holdup_Liquido = 0.70
m_dot_m = 5"""

# T = 353.15        # K
# P = 44597790      # Pa
# rho_o = 584 #kg/m3
# Mg = 25.55        # kg/kmol
# QO = 0.009576 # m3/s
# QW = 0.001525 # m3/s
# QL = QO + QW # m3/s
# Holdup_Liquido = 0.5169
# m_dot_m = 9.935593 # kg/s
#
# theta = 60 # graus
#
# dl = 12.0317 # m
# depth = -2.9700686354772188e-12 # m
#
# # Temperature(T, P, rho_o, QO, QW, QL, Mg, Holdup_Liquido)
# T = Temperature(T, P, rho_o, QO, QW, QL, Mg, Holdup_Liquido).Calculate_Temperature(dl, depth, m_dot_m, theta)
#
# print(T)
# # Calculate_Temperature(self, dl, depth, m_dot_m, theta)"""
