from numpy import zeros, sqrt, mean
from matplotlib.pyplot import xlim, ylim, figure, show
from random import randrange
import auxiliarFunctions as auxF

#---------------------------------------- CONSTANTS DEFINITION ----------------------------------------#


iterations = 1000000  # Number of Monte Carlo iterations
nParticles = 99       # Number of particles on the material
U01 = 1               # Lennard - Jones potential constant from material
r01 = 0.5
with open("r01.txt", "r") as file:
    req = float(file.readline())
print(req)
N = 17     # N should be a divisor of nParticles so that the dimensions of the box are correctly specified
a = nParticles//(2*N - 1)
h = sqrt(3)/2*req
H = h*(2*a - 1)
beta = 100            # Related to temperature constant
beta = 100            # Related to temperature constant

x1 = 0              # wall at x = 0
x2 = (N-1)*req      # wall at x = 2
y1 = 0              # wall at y = 0
y2 = 9              # wall at y = 9

position = zeros([2, nParticles])       # Material particles positions
increments = zeros([2, nParticles])     # Material particles increments
energies = zeros([iterations + 1])      # Energy values will be stored here when computed
auxiliar = []

for i in range(N):     # Starting positions of the material particles
    position[0][i] = req*i
    position[1][i] = H

for i in range(N - 1):
    position[0][i + N] = req*i + req/2
    position[1][i + N] = H - h

for i in range(nParticles - 2*N + 1):
    position[0][i + 2*N - 1] = position[0][i]
    position[1][i + 2*N - 1] = position[1][i] - 2*h

X = position[0]  # Every x position
Y = position[1]  # Every y position

fig = figure()
ax = fig.add_subplot(111)
system, = ax.plot(X, Y, ".")
xlim([-0.5, 9.2])
ylim([-0.1, 9])
show()

E = auxF.interactionEnergy(position, nParticles, U01, r01)
energies[0] = E

for i in range(iterations):     # Computes the changes in the system and shows the simulation
    index = randrange(0, nParticles)
    auxF.sortModifyIncrement(increments, index, r01)
    auxF.boundryControl(position, increments, index, x1, x2, y1, y2)

    dE = auxF.dEi(position, increments, index, nParticles, U01, r01)

    if (auxF.probability(beta, dE)):
        auxF.modifyPos(position, increments, index)
    else:
        dE = 0
    energies[i+1] = energies[i] + dE

    if (i > 800000):
        auxiliar.append(energies[i+1])

    if (i % 10000 == 0):
        print(i)

X = position[0]  # Every x position
Y = position[1]  # Every y position
fig = figure()
ax = fig.add_subplot(111)
system, = ax.plot(X, Y, ".")
xlim([-0.5, 9.2])
ylim([-0.1, 9])

fig2 = figure()
ax2 = fig2.add_subplot(111)
ax2.plot(energies)

print(mean(auxiliar))
fig3 = figure()
ax3 = fig3.add_subplot(111)
ax3.plot(auxiliar)
show()