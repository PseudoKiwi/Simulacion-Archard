# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from matplotlib.pyplot import xlim, ylim, figure, ion, show
from numpy import zeros

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    position = zeros([2, 100])
    energies = zeros([1000001])
    with open("prueba.txt", "r") as file:
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
    show()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
