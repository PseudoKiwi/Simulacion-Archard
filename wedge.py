from numpy import zeros
import auxiliarFunctions as auxF
from random import random, randrange
from matplotlib.pyplot import xlim, ylim, figure, show

#---------------------------------------- CONSTANTS DEFINITION ----------------------------------------#

iterations = 300000
N = 20       # Wedge base, fixed nodes
nParticles = int(N * (N + 1) / 4)
r02 = 0.1   # Equilibrium radius for wedge particles
U02 = 1
beta = 100
accepted = [0]      # Will store the amount of times a change was accepted when dE > 0

acceptance2 = auxF.defineProbability(accepted, beta)        # acceptance is the probability function for change acceptance
interactionType2 = auxF.potential(U02, r02)                       # Interaction potential between wedge particles

position = zeros([2, nParticles])            # Material particles positions
increments = zeros([2, nParticles])     # Material particles increments
energies = zeros([iterations + 1])   # Energy values will be stored here when computed


#--------------------------------------------- SIMULATION ---------------------------------------------#


for i in range(N):
    position[0][i] = i*r02
    position[1][i] = N * r02

for i in range(N, nParticles):     # Starting positions of the material particles
    position[0][i] = random() * (N - 14) * r02 + position[0][0] + 7 * r02
    position[1][i] = position[1][0] - r02

E = auxF.interactionEnergy(position, nParticles, interactionType2)
energies[0] = E

X = position[0]  # Every x position
Y = position[1]  # Every y position

fig = figure()
ax = fig.add_subplot(111)
system, = ax.plot(X, Y, ".")
xlim([-r02, N * r02])
ylim([-0.1, N * r02])

for i in range(iterations):     # Computes the changes in the system and shows the simulation
    index = randrange(N + 1, nParticles)
    auxF.sortModifyIncrement(increments, index, r02)

    dE = auxF.dEi(position, increments, index, nParticles, interactionType2)

    if (acceptance2(dE)):
        auxF.modifyPos(position, increments, index)
    else:
        dE = 0
    energies[i + 1] = energies[i] + dE


X = position[0]  # Every x position
Y = position[1]  # Every y position
fig2 = figure()
ax2 = fig2.add_subplot(111)
system, = ax2.plot(X, Y, ".")
xlim([-r02, N * r02])
ylim([-0.1, N * r02])

fig3 = figure()
ax3 = fig3.add_subplot(111)
ax3.plot(energies)
show()
