from numpy import zeros, sqrt
import auxiliarFunctions as auxF
from random import randrange
from matplotlib.pyplot import xlim, ylim, figure, show

#---------------------------------------- CONSTANTS DEFINITION ----------------------------------------#

iterations = 100000
N = 8        # Wedge base, fixed nodes
nParticles = N * (N + 1) // 2
r02 = 0.1    # Equilibrium radius for wedge particles
U02 = 1
beta = 100
with open("r02.txt", "r") as file:
    req = float(file.readline())
h = sqrt(3) / 2 * req

position = zeros([2, nParticles])            # Material particles positions
increments = zeros([2, nParticles])     # Material particles increments
energies = zeros([iterations + 1])   # Energy values will be stored here when computed


#--------------------------------------------- SIMULATION ---------------------------------------------#


aux = N
e = 0
f = 1
mod = False
for i in range(nParticles):
    position[0][i] = (f % aux)*req + e*req/2
    position[1][i] = (N-1) * h - e*h
    if f % aux == 0 and f // aux == 1:
        aux -= 1
        e += 1
        mod = True
    if mod:
        f = 1
        mod = False
    else:
        f += 1

E = auxF.interactionEnergy(position, nParticles, U02, r02)
energies[0] = E

for e in range(4):
    for i in range(iterations):     # Computes the changes in the system and shows the simulation

        if i % 10000 == 0:
            print(i)

        index = randrange(N, nParticles)
        auxF.sortModifyIncrement(increments, index, r02)

        dE = auxF.dEi(position, increments, index, nParticles, U02, r02)

        if auxF.probability(beta, dE):
            auxF.modifyPos(position, increments, index)
        else:
            dE = 0
        energies[i + 1] = energies[i] + dE

    X = position[0]  # Every x position
    Y = position[1]  # Every y position
    fig2 = figure()
    ax2 = fig2.add_subplot(111)
    system, = ax2.plot(X, Y, ".")
    xlim([-req, N * req])
    ylim([-3 * h, N * h])
    show()

    for i in range(N):
        position[1][i] -= h/2

    print(e)

fig3 = figure()
ax3 = fig3.add_subplot(111)
ax3.plot(energies)
show()
