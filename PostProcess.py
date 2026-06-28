from matplotlib import pyplot as plt

def run(pressure=None, position=None, vm=None,temperature=None):
    if pressure is not None:
        if isinstance(pressure[0], list):
            for i in range(len(pressure)):
                plt.plot(position,pressure[i])
        else:
            plt.plot(position, pressure)
        plt.grid(True)
        plt.xlabel('Position')
        plt.ylabel('Pressure')
        plt.show()

    if vm is not None:
        plt.plot(position,vm)
        plt.grid(True)
        plt.xlabel('Position')
        plt.ylabel('Velocity')
        plt.show()
    if temperature is not None:
        plt.plot(position, temperature)
        plt.grid(True)
        plt.xlabel('Position')
        plt.ylabel('Temperature')
        plt.show()