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
    e = zeros([150])
    p = zeros([150])

    with open("datosEq2.txt", "r") as file:
        a = file.readline()
        b = file.readline()

        a = a.split(',')
        a[0] = a[0].split('[')[1]
        a[-1] = a[-1].split(']')[0]

        b = b.split(',')
        b[0] = b[0].split('[')[1]
        b[-1] = b[-1].split(']')[0]

        e[:] = a.copy()[:150]
        p[:] = b.copy()[:150]


    fig3 = figure()
    ax3 = fig3.add_subplot(111)
    system, = ax3.plot(e)

    fig3 = figure()
    ax3 = fig3.add_subplot(111)
    system, = ax3.plot(p)

    show()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
