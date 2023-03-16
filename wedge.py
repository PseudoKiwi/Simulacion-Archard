from numpy import pi, sin, cos, sqrt, inf, zeros
from probability import defineProbability
from potential import potential
from random import random, randrange
from matplotlib.pyplot import xlim, ylim, figure, ion, show

#---------------------------------------- CONSTANTS DEFINITION ----------------------------------------#

iterations = 1000
N = 5       # Wedge base, fixed nodes
nParticles = int(N * (N + 1) / 2)
r02 = 0.5   # Equilibrium radius for wedge particles
U02 = 1000
beta = 100
dENeg = [0]         # Will store the amount of times dE < 0
accepted = [0]      # Will store the amount of times a change was accepted when dE > 0

acceptance2 = defineProbability(dENeg, accepted, beta)        # acceptance is the probability function for change acceptance
interactionType2 = potential(U02, r02)                       # Interaction potential between wedge particles

position = zeros([2, nParticles])            # Material particles positions
increments = zeros([2, nParticles])     # Material particles increments
energies = zeros([1, iterations + 1])   # Energy values will be stored here when computed


#----------------------------------------- AUXILIAR FUNCTIONS -----------------------------------------#


def particleInteractionEnergy(pos, inc, i, withIncrements = False):  # Computes the total interaction energy of particle i
    # When increments = True, it considers de future positions
    Eij = []
    if (not withIncrements):
        ri = [pos[0][i], pos[1][i]]
        for j in range(nParticles):
            if (i != j):
                rj = [pos[0][j], pos[1][j]]
                r = sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2)
                Eij.append(interactionType2(r))
    else:
        ri = [pos[0][i] + inc[0][i], pos[1][i] + inc[1][i]]
        for j in range(nParticles):
            if (i != j):
                rj = [pos[0][j] + inc[0][j], pos[1][j] + inc[1][j]]
                r = sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2)
                Eij.append(interactionType2(r))
    return sum(Eij)


def sortModifyIncrement2(inc, i):   # Sorts and modify de increment of particle i
    angle = 2 * pi * random()
    dr = r02/5 * random()
    inc[0][i] = dr * cos(angle)
    inc[1][i] = dr * sin(angle)


def modifyPos(pos, inc, i): # Modifies the position of particle i
    pos[0][i] += inc[0][i]
    pos[1][i] += inc[1][i]


#--------------------------------------------- SIMULATION ---------------------------------------------#


for i in range(N):
    position[0][i] = i*r02
    position[1][i] = N * r02

for i in range(N, nParticles):     # Starting positions of the material particles
    position[0][i] = random() * (N - 5) * r02 + position[0][0] + 2 * r02
    position[1][i] = random() * (4 - N) * r02 + position[1][0] - r02

E = 0
for i in range(nParticles):
    E += particleInteractionEnergy(position, increments, i)

energies[0] = E

X = position[0]  # Every x position
Y = position[1]  # Every y position

ion()
fig = figure()
ax = fig.add_subplot(111)
system, = ax.plot(X, Y, ".")
xlim([-r02, N * r02])
ylim([-0.1, N * r02])

for i in range(iterations):     # Computes the changes in the system and shows the simulation
    index = randrange(N + 1, nParticles)
    sortModifyIncrement2(increments, index)
    Eactual = particleInteractionEnergy(position, increments, index)      # Actual interaction energy of chosen particle
    Emod = particleInteractionEnergy(position, increments, index, True)   # Future interaction energy of chosen particle
    dE = Emod - Eactual
    if (acceptance2(dE)):
        E += dE
        modifyPos(position, increments, index)
        X = position[0]  # Every x position
        Y = position[1]  # Every y position
    energies[0][i+1] = E
    system.set_xdata(X)
    system.set_ydata(Y)
    fig.canvas.draw()
    fig.canvas.flush_events()

fig2 = figure()
ax2 = fig2.add_subplot(111)
ax2.plot(energies[0])
show()

print(dENeg[0])
print(accepted[0])
AR = (accepted[0] / (iterations - dENeg[0]))*100
print(AR)
