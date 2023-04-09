from numpy import zeros, abs, mean, sqrt
from random import randrange
from matplotlib.pyplot import xlim, ylim, figure, show
import auxiliarFunctions as auxF
from sys import maxsize

# ---------------------------------------- CONSTANTS DEFINITION ----------------------------------------#

iterations = 1000000  # Number of Monte Carlo iterations
nParticles = 99  # Number of particles on the material

with open("r01.txt", "r") as file:
    req = float(file.readline())

N = 17
h = sqrt(3) / 2 * req
a = nParticles // (2 * N - 1)
H = h * (2 * a - 1)

x1 = 0  # wall at x = 0
x2 = (N - 1) * req  # wall at x = 2
y1 = 0  # wall at y = 0
y2 = H  # wall at y = H

ly = y2 - y1
lx = x2 - x1

yi = H
yf = yi - 80 / 100 * ly
dy = ly / 100
dx = dy * lx / (2 * (ly - dy))

expData = int(abs(yi - yf) / dy)  # Number of total iterations of the simulation

eqEnergies = zeros([expData + 1])
eqPressures = zeros([expData + 1])
eqLength = zeros([expData + 1])

U01 = 1  # Lennard - Jones potential constant from material
r01 = 0.5  # Equilibrium radius from material

beta = 100  # Related to temperature constant
accepted = [0]  # Will store the amount of times a change was accepted when dE > 0


position = zeros([2, nParticles])  # Material particles positions
increments = zeros([2, nParticles])  # Material particles increments
energies = zeros([iterations + 1])  # Energy values will be stored here when computed
eq = iterations // 10

# ---------------------------------------- SIMULATION ----------------------------------------#

for i in range(N):  # Starting positions of the material particles
    position[0][i] = req * i
    position[1][i] = H

for i in range(N - 1):
    position[0][i + N] = req * i + req / 2
    position[1][i + N] = H - h

for i in range(nParticles - 2 * N + 1):
    position[0][i + 2 * N - 1] = position[0][i]
    position[1][i + 2 * N - 1] = position[1][i] - 2 * h

E = auxF.interactionEnergy(position, nParticles, U01, r01)
energies[0] = E

pressure = 0
minY = maxsize
for i in range(N):
    minY = min(position[1][i], minY)

aux = []
for i, y in enumerate(position[1]):
    if y > minY:
        aux.append(i)

for e in aux:
    vFef = 0
    for f in range(nParticles):
        if position[1][f] < minY:
            vFef += auxF.verticalInteractionForce(U01, r01, position[0][e], position[0][f], position[1][e], position[1][f])
    pressure += vFef/lx

reqC = zeros([nParticles])
for i in range(nParticles):
    minimum = maxsize
    for j in range(nParticles):
        if i != j:
            r = sqrt((position[0][i] - position[0][j]) ** 2 + (position[1][i] - position[1][j]) ** 2)
            if r < minimum:
                minimum = r
    reqC[i] = minimum

eqEnergies[0] = E
eqPressures[0] = pressure
eqLength[0] = mean(reqC)

for h in range(expData):

    fig1 = figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(position[0], position[1], ".")
    xlim([-1, 10])
    ylim([-0.1, 10.9])
    show()

    auxE = []       # Energy values after eq. will be stored here when computed
    pressures = []  # Pressure values after eq. will be stored here when computed

    print(h + 1)
    for i in range(iterations):     # Computes the changes in the system and shows the simulation
        index = randrange(0, nParticles)
        auxF.sortModifyIncrement(increments, index, r01)
        auxF.boundryControl(position, increments, index, x1, x2, y1, y2)

        dE = auxF.dEi(position, increments, index, nParticles, U01, r01)

        if auxF.probability(beta, dE):
            auxF.modifyPos(position, increments, index)
        else:
            dE = 0
        energies[i+1] = energies[i] + dE

        if i >= iterations - eq:
            pressure = 0
            minY = sys.maxsize
            for e in range(N):
                minY = min(position[1][e], minY)

            aux = []
            for e, y in enumerate(position[1]):
                if y > minY:
                    aux.append(e)

            for e in aux:
                vFef = 0
                for f in range(nParticles):
                    if position[1][f] < minY:
                        vFef += auxF.verticalInteractionForce(U01, r01, position[0][e], position[0][f], position[1][e],
                                                              position[1][f])
                pressure += vFef / lx
            pressures.append(pressure)
            auxE.append(energies[i+1])

        if i % 10000 == 0:
            print(i)

    y2 -= dy
    #x1 -= dx
    #x2 += dx
    auxF.totalBoundryControl(position, x1, x2, y1, y2)

    fig3 = figure()
    ax3 = fig3.add_subplot(111)
    ax3.plot(auxE)

    fig4 = figure()
    ax4 = fig4.add_subplot(111)
    ax4.plot(pressures)

    print([mean(auxE), mean(pressures)])

    eqEnergies[h + 1] = mean(auxE)
    eqPressures[h + 1] = mean(pressures)

    if eqPressures[h] < 0:
        eqPressures[h] = 0

    reqC = zeros([nParticles])
    for i in range(nParticles):
        minimum = maxsize
        for j in range(nParticles):
            if i != j:
                r = sqrt((position[0][i] - position[0][j]) ** 2 + (position[1][i] - position[1][j]) ** 2)
                if r < minimum:
                    minimum = r
        reqC[i] = minimum
    eqLength[h + 1] = mean(reqC)

with open("resultadosEq2.txt", "w") as file:
    str1 = str(list(eqEnergies))
    str2 = str(list(eqPressures))
    str3 = str(list(eqLength))
    file.write(str1 + "\n")
    file.write(str2 + "\n")
    file.write(str3 + "\n")

fig3 = figure()
ax3 = fig3.add_subplot(111)
ax3.plot(eqEnergies)

fig4 = figure()
ax4 = fig4.add_subplot(111)
ax4.plot(eqPressures)

fig5 = figure()
ax5 = fig5.add_subplot(111)
ax5.plot(eqLength)
show()