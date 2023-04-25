from numpy import zeros, abs, mean, sqrt, inf
from random import randrange
from matplotlib.pyplot import xlim, ylim, figure, show
import auxiliarFunctions as auxF
from sys import maxsize

# ---------------------------------------- CONSTANTS DEFINITION ----------------------------------------#
# ---------------------------------------------- GENERAL -----------------------------------------------#

iterations = 1000000  # Number of Monte Carlo iterations
mParticles = 99  # Number of particles on the material
wBase = 8
wParticles = wBase*(wBase + 1)//2
nParticles = mParticles + wParticles
beta = 100  # Related to temperature constant

with open("r01.txt", "r") as file:
    aux = float(file.readline())
    req1 = aux
    req2 = aux
    reqInt = aux

r0Int = 0.5
U0Int = 1

position = zeros([2, nParticles])    # Particles positions
increments = zeros([2, nParticles])  # Particles increments
energies = zeros([iterations + 1])   # Energy values will be stored here when computed

x1 = 0  # wall at x = 0
x2 = 16 * req1  # wall at x = 2
y1 = 0  # wall at y = 0
y2 = inf  # wall at y = H

# ---------------------------------------------- MATERIAL ----------------------------------------------#


a = 17
h1 = sqrt(3) / 2 * req1
b = mParticles // (2 * a - 1)
H = h1 * (2 * b - 1)

r01 = 0.5  # Equilibrium radius from material
U01 = 1  # Lennard - Jones potential constant from material


# ----------------------------------------------- WEDGE -----------------------------------------------#


N = 8        # Wedge base, fixed nodes
h2 = sqrt(3) / 2 * req2

r02 = 0.5    # Equilibrium radius for wedge particles
U02 = 1


#--------------------------------------------- SIMULATION ---------------------------------------------#

# Starting positions of the material particles
for i in range(a):
    position[0][i + wParticles] = req1*i
    position[1][i + wParticles] = H

for i in range(a - 1):
    position[0][i + a + wParticles] = req1*i + req1/2
    position[1][i + a + wParticles] = H - h1

for i in range(mParticles - 2*a + 1):
    position[0][i + 2*a - 1 + wParticles] = position[0][i + wParticles]
    position[1][i + 2*a - 1 + wParticles] = position[1][i + wParticles] - 2*h1

# Starting positions of the wedge particles (triangle form)
dx = (a - 1)*req1/2 - (N - 1)*req2/2
dy = H + 2*reqInt
aux = N
e = 0
f = 1
mod = False
for i in range(wParticles):
    position[0][i] = (f % aux)*req2 + e*req2/2 + dx
    position[1][i] = (N-1) * h2 - e*h2 + dy
    if f % aux == 0 and f // aux == 1:
        aux -= 1
        e += 1
        mod = True
    if mod:
        f = 1
        mod = False
    else:
        f += 1

X1 = position[0][wParticles:]  # Every x material position
Y1 = position[1][wParticles:]  # Every y material position
X2 = position[0][:wParticles]  # Every x wedge position
Y2 = position[1][:wParticles]  # Every y wedge position
fig2 = figure()
ax2 = fig2.add_subplot(111)
system, = ax2.plot(X1, Y1, ".")
system, = ax2.plot(X2, Y2, ".")
xlim([-req1, req1*a])
ylim([-0.3, H + N*h2 + 1])
show()

E = auxF.interactionEnergy2(position, mParticles, wParticles, U01, r01, U02, r02, U0Int, r0Int)
energies[0] = E

parcial = 8
for e in range(parcial):

    for i in range(N):
        position[1][i] -= reqInt / 4

    print(e)

    for i in range(iterations):     # Computes the changes in the system and shows the simulation

        if i % 10000 == 0:
            print(i)

        index = randrange(N, nParticles)
        if index >= wParticles:
            auxF.sortModifyIncrement(increments, index, r01)
        else:
            auxF.sortModifyIncrement(increments, index, r02)
        auxF.boundryControl(position, increments, index, x1, x2, y1, y2)

        dE = auxF.dEi2(position, increments, index, mParticles, wParticles, U01, r01, U02, r02, U0Int, r0Int)

        if auxF.probability(beta, dE):
            auxF.modifyPos(position, increments, index)
        else:
            dE = 0
        energies[i + 1] = energies[i] + dE

    X1 = position[0][wParticles:]  # Every x material position
    Y1 = position[1][wParticles:]  # Every y material position
    X2 = position[0][:wParticles]  # Every x wedge position
    Y2 = position[1][:wParticles]  # Every y wedge position
    fig2 = figure()
    ax2 = fig2.add_subplot(111)
    system, = ax2.plot(X1, Y1, ".")
    system, = ax2.plot(X2, Y2, ".")
    xlim([-req1, req1*a])
    ylim([-0.3, H + N*h2 + 1])
    show()

for e in range(parcial + 4):

    for i in range(N):
        position[1][i] += reqInt / 4

    print(e)

    for i in range(iterations):  # Computes the changes in the system and shows the simulation

        if i % 10000 == 0:
            print(i)

        index = randrange(N, nParticles)
        if index >= wParticles:
            auxF.sortModifyIncrement(increments, index, r01)
        else:
            auxF.sortModifyIncrement(increments, index, r02)
        auxF.boundryControl(position, increments, index, x1, x2, y1, y2)

        dE = auxF.dEi2(position, increments, index, mParticles, wParticles, U01, r01, U02, r02, U0Int, r0Int)

        if auxF.probability(beta, dE):
            auxF.modifyPos(position, increments, index)
        else:
            dE = 0
        energies[i + 1] = energies[i] + dE

    X1 = position[0][wParticles:]  # Every x material position
    Y1 = position[1][wParticles:]  # Every y material position
    X2 = position[0][:wParticles]  # Every x wedge position
    Y2 = position[1][:wParticles]  # Every y wedge position
    fig2 = figure()
    ax2 = fig2.add_subplot(111)
    system, = ax2.plot(X1, Y1, ".")
    system, = ax2.plot(X2, Y2, ".")
    xlim([-req1, req1 * a])
    ylim([-0.3, H + N * h2 + 1])
    show()