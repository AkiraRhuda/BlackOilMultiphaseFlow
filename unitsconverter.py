"""
Módulo para conversão de unidades de temperatura e pressão.

Métodos:
----------
    Temperature(value, from_unity, to_unity): Converte valores de temperatura
    Pressure(value, from_unity, to_unity): Converte valores de pressão
    Density(value, from_unity, to_unity): Converte valores de massa específica

Parâmetros:
----------
        value (float): Valor numérico da temperatura a ser convertida
        from_unity (str): Unidade de origem ('C', 'F', 'K' ou 'R')
        to_unity (str): Unidade de destino ('C', 'F', 'K' ou 'R')

Retorna:
--------
        float: Valor convertido na unidade desejada

    Escalas suportadas:
    - Temperatura:
        - Celsius (°C)
        - Fahrenheit (°F)
        - Kelvin (K)
        - Rankine (°R)
    - Pressão:
        - psi (libra-força por polegada quadrada)
        - psia (psi absoluto - mesmo valor numérico que psi)
        - Pa (Pascal)
        - bar
        - atm (atmosfera padrão)
    - Massa especifífica:
        - °API (grau API)
        - d (densidade relativa)
        - % (porcentagem de teor de sólidos)
        - mg/L (miligramas por litro)
    """


def Temperature(value, from_unity, to_unity):

    if from_unity == 'C':
        if to_unity == 'F':
            return (value * 9 / 5) + 32
        elif to_unity == 'K':
            return value + 273.15
        elif to_unity == 'R':
            return (value + 273.15) * 9 / 5


    elif from_unity == 'F':
        if to_unity == 'C':
            return (value - 32) * 5 / 9
        elif to_unity == 'K':
            return (value - 32) * 5 / 9 + 273.15
        elif to_unity == 'R':
            return value + 459.67


    elif from_unity == 'K':
        if to_unity == 'C':
            return value - 273.15
        elif to_unity == 'F':
            return (value - 273.15) * 9 / 5 + 32
        elif to_unity == 'R':
            return value * 9 / 5


    elif from_unity == 'R':
        if to_unity == 'C':
            return (value - 491.67) * 5 / 9
        elif to_unity == 'F':
            return value - 459.67
        elif to_unity == 'K':
            return value * 5 / 9

    if from_unity == to_unity:
        return value

    raise ValueError(f"Conversão de {from_unity} para {to_unity} não suportada")


def Pressure(value, from_unity, to_unity):

    if to_unity == 'psi':
        if from_unity == 'Pa':
            return value / 6894.76
        elif from_unity == 'bar':
            return value * 14.5038
        elif from_unity == 'atm':
            return value * 14.6959
        elif from_unity == 'psi':
            return value


    elif from_unity == 'psi':
        if to_unity == 'Pa':
            return value * 6894.76
        elif to_unity == 'bar':
            return value / 14.5038
        elif to_unity == 'atm':
            return value / 14.6959


    elif from_unity == 'bar':
        if to_unity == 'atm':
            return value * 0.986923
        elif to_unity == 'Pa':
            return value * 100000

    if from_unity == to_unity:
        return value

    raise ValueError(f"Conversão de {from_unity} para {to_unity} não suportada")

def Density(value, from_unity, to_unity):

    if from_unity == 'API':
        if to_unity == 'do':
            return 141.5/(value+131.5) # do=141.5/(api+131) ----- do((api+131)) = 141.5


    elif from_unity == 'do':
        if to_unity == 'API':
            return 141.5/value - 131.5

    elif from_unity =='%':
        if to_unity == 'mg/L':
            return value[0] * 10000 * value[1]

    elif from_unity == 'lbm/ft3':
        if to_unity == 'do':
            # Valor de referência é água doce a 60°F
            return value/62.428
        if to_unity == 'kg/m3':
            return value * 16.01846337

    elif from_unity == 'kg/m3':
        if to_unity == 'do':
            return value/1000
        if to_unity == 'dg':
            return value / 1.2754
    if from_unity == to_unity:
        return value

    raise ValueError(f"Conversão de {from_unity} para {to_unity} não suportada")

def Length(value, from_unity, to_unity):
    if from_unity == 'in':
        if to_unity == 'ft':
            return value/12
    if from_unity == 'ft':
        if to_unity == 'in':
            return value*12
    if from_unity == 'ft':
        if to_unity == 'm':
            return value * 0.3048
    if from_unity == 'in':
        if to_unity == 'm':
            return value * 0.0254
    if from_unity == 'm':
        if to_unity == 'ft':
            return value / 0.3048
    if from_unity == 'm':
        if to_unity == 'in':
            return value / 0.0254

    raise ValueError(f"Conversão de {from_unity} para {to_unity} não suportada")

def Viscosity(value, from_unity, to_unity):
    if from_unity == 'cp':
        if to_unity == 'Pa.s':
            return value * 10**-3

    raise ValueError(f"Conversão de {from_unity} para {to_unity} não suportada")

def Standard(value, from_unity, to_unity):

    if from_unity == 'sm3sm3':
        if to_unity == 'scfstb':
            return value / 0.1781076067
        if to_unity == 'bbl/stb':
            return value * 0.1781076067

    elif from_unity == 'scfstb':
        if to_unity == 'sm3sm3':
            return value * 0.1781076067
    elif from_unity == 'bbl/stb':
        if to_unity == 'sm3/sm3':
            return value * 0.1781076067

    raise ValueError(f"Conversão de {from_unity} para {to_unity} não suportada")