from matplotlib import pyplot as plt

def run(pressure, position, vm):
    plt.plot(position,pressure)
    plt.grid(True)
    plt.xlabel('Position')
    plt.ylabel('Pressure')
    plt.show()

    plt.plot(position,vm)
    plt.grid(True)
    plt.xlabel('Position')
    plt.ylabel('Velocity')
    plt.show()