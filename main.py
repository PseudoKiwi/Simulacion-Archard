# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import matplotlib.pyplot as plt
from numpy import zeros, polyfit, linspace, genfromtxt, unique, where, std

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
    #with open("force-loss2.txt", "r") as file:
    #    for linea in file:
    #        aux = linea.split(', ')
    #        particulas.append(float(aux[1]))
    #        fuerza.append(float(aux[0]))

    #L = []
    #for i in range(len(fuerza)):
    #    L.append(0.5437662932839045/4*(1 + i % 12))

    #m, b = polyfit(fuerza, particulas, 1)
    #p2 = linspace(min(fuerza), max(fuerza), 2)
    #f2 = b + m * p2

    #n, v = polyfit(fuerza, L, 1)
    #p3 = linspace(min(fuerza), max(fuerza), 2)
    #f3 = v + n * p3

    #fig = figure()
    #ax = fig.add_subplot(111)
    #system, = ax.plot(fuerza, particulas, ".")
    #system, = ax.plot(p2, f2)

    #fig = figure()
    #ax = fig.add_subplot(111)
    #system, = ax.plot(fuerza, L, ".")
    #system, = ax.plot(p3, f3)

    L = zeros([12], float)
    for i in range(12):
        L[i] = 0.5437662932839045 / 4 * (1 + i % 12)
    fuerzasL = zeros([12], float)

    with open("force-loss2.txt", "r") as file:
        cont = 0
        for linea in file:
            aux = linea.split(', ')
            fuerzasL[cont % 12] += float(aux[0])/5
            cont += 1

    # Leer los datos del archivo
    data = genfromtxt('force-loss.txt', delimiter=',')
    numbers = data[:, 0]
    categories = data[:, 1]

    indices = where(numbers <= 30)
    numbers = numbers[indices]
    categories = categories[indices]

    # Obtener los números de la derecha únicos
    unique_categories = unique(categories)

    prom = []
    std_dev = []
    # Calcular el promedio y la desviación estándar para cada categoría
    for category in unique_categories:
        indices = where(categories == category)
        numbers_for_category = numbers[indices]
        prom.append(mean(numbers_for_category))
        std_dev.append(std(numbers_for_category))


    m, b = polyfit(unique_categories, prom, 1)
    p2 = linspace(min(unique_categories), max(unique_categories), 2)
    f2 = b + m * p2

    #n, v = polyfit(fuerzasL, L, 1)
    #p3 = linspace(min(fuerzasL), max(fuerzasL), 2)
    #f3 = v + n * p3

    fig, ax = plt.subplots()
    ax.errorbar(unique_categories, prom, yerr=std_dev, fmt=".")
    system, = ax.plot(p2, f2)

    #fig, ax = plt.subplots()
    #ax.plot(fuerzasL, L, ".")
    #ax.plot(p3, f3)

    plt.show()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
