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
    #e = zeros([30])
    #p = zeros([30])

    #with open("resultadosEq3.txt", "r") as file:
    #    a = file.readline()
    #    b = file.readline()

    #    a = a.split(',')
    #    a[0] = a[0].split('[')[1]
    #    a[-1] = a[-1].split(']')[0]

    #    b = b.split(',')
    #    b[0] = b[0].split('[')[1]
    #    b[-1] = b[-1].split(']')[0]

    #    e[:] = a.copy()[:30]
    #    p[:] = b.copy()[:30]


    fig3 = figure()
    ax3 = fig3.add_subplot(111)
    system, = ax3.plot([2,4,4,4,6,7,6,6,6,16,14,19],[8.196310363416938,12.42323983598615,13.150891507499873,13.837511387561603,19.735858338557033,21.333120170063086,23.09459734070445,32.17306003810418,68.3456504108901,64.69017723791073,61.013339448827146,66.2441167672423])

    #fig3 = figure()
    #ax3 = fig3.add_subplot(111)
    #system, = ax3.plot(p)

    show()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
