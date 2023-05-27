# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from matplotlib.pyplot import figure, show
from numpy import zeros, polyfit, linspace

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

    #particulas = []
    #fuerza = []
    #with open("force-loss.txt", "r") as file:
    #    for linea in file:
    #        aux = linea.split(', ')
    #        particulas.append(float(aux[1]))
    #        fuerza.append(float(aux[0]))
    #m, b = polyfit(particulas, fuerza, 1)
    #p2 = linspace(min(particulas), max(particulas), 2)
    #f2 = b + m*p2

    cantidad = zeros([50])
    fuerzas = zeros([50], float)
    particulas = []

    with open("force-loss.txt", "r") as file:
        for linea in file:
            aux = linea.split(', ')
            p = int(aux[1])
            fuerzas[p] += float(aux[0])
            cantidad[p] += 1

    fuerza = []
    for i, f in enumerate(fuerzas):
        if f != 0:
            c = cantidad[i]
            fuerza.append(f/c)
            particulas.append(i)

    m, b = polyfit(particulas, fuerza, 1)
    p2 = linspace(min(particulas), max(particulas), 2)
    f2 = b + m * p2

    fig3 = figure()
    ax3 = fig3.add_subplot(111)
    system, = ax3.plot(particulas, fuerza, ".")
    system, = ax3.plot(p2, f2)

    #fig3 = figure()
    #ax3 = fig3.add_subplot(111)
    #system, = ax3.plot(p)

    show()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
