# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from matplotlib.pyplot import xlim, ylim, figure, ion, show
from numpy import zeros
from auxiliarFunctions import interactionEnergy, potential

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.

from numpy import mean
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    position = zeros([2, 100])
    energies = zeros([5000001])
    with open("prueba2.txt", "r") as file:
        a = file.readline()
        b = file.readline()
        c = file.readline()

        a = a.split(',')
        a[0] = a[0].split('[')[1]
        a[-1] = a[-1].split(']')[0]

        b = b.split(',')
        b[0] = b[0].split('[')[1]
        b[-1] = b[-1].split(']')[0]

        c = c.split(',')
        c[0] = c[0].split('[')[1]
        c[-1] = c[-1].split(']')[0]

        position[0][:] = a
        position[1][:] = b
        energies[:] = c

    X = position[0]
    Y = position[1]

    fig1 = figure()
    ax1 = fig1.add_subplot(111)
    system, = ax1.plot(X, Y, ".")
    xlim([-3, 5])
    ylim([-0.1, 5])

    fig2 = figure()
    ax2 = fig2.add_subplot(111)
    system, = ax2.plot(energies, ".")

    fig3 = figure()
    ax3 = fig3.add_subplot(111)
    system, = ax3.plot(energies[4000000:], ".")
    show()

    print(mean(energies[4000000:]))

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
