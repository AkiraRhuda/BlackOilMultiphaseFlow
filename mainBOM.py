from BlackOilModel import GasPhase
from BlackOilModel import OilPhase
import unitsconverter

# Dados iniciais
densidaderelativagas = 0.84
densidaderelativaoleo = 0.86
tempseparador = 80 # °F
Pressaoseparador = 100 # psia
Pressaodebolha = 5000 # psia
pressaolinha = 3626 # psia
temperaturalinha = 122 # °F

GasPhase(zcorrelation="Papay", μcorrelation='Lee', dg=densidaderelativagas, P=pressaolinha, T=unitsconverter.Temperature(temperaturalinha, 'F','R')).output()
OilPhase("Standing", "Standing", "Beggs",densidaderelativaoleo, densidaderelativagas, pressaolinha, Pressaodebolha, temperaturalinha).output()