from matplotlib import pyplot as plt

def run(pressure, position):
    plt.plot(position,pressure)
    plt.grid(True)
    plt.xlabel('Position')
    plt.ylabel('Pressure')
    plt.show()