from matplotlib import pyplot as plt

def run(pressure, position, vm,temperature):
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

    plt.plot(position, temperature)
    plt.grid(True)
    plt.xlabel('Position')
    plt.ylabel('Temperature')
    plt.show()