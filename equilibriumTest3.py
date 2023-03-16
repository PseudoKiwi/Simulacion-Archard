from numpy import pi, sin, cos, sqrt, inf, zeros
from probability import defineProbability
from potential import potential
from random import random, randrange
from matplotlib.pyplot import xlim, ylim, figure, ion, show


#---------------------------------------- CONSTANTS DEFINITION ----------------------------------------#


iterations = 1000   # Number of Monte Carlo iterations
nParticles = 100    # Number of particles on the material
U01 = 1             # Lennard - Jones potential constant from material
r01 = 0.5           # Equilibrium radius from material
beta = 100          # Related to temperature constant   
dENeg = [0]         # Will store the amount of times dE < 0
accepted = [0]      # Will store the amount of times a change was accepted when dE > 0

acceptance = defineProbability(dENeg, accepted, beta)        # acceptance is the probability function for change acceptance
interactionType1 = potential(U01, r01)                       # Interaction potential between material particles

x1 = -2   # wall at x = -2
x2 = 2    # wall at x = 2
y1 = 0    # wall at y = 0
y2 = 2    # wall at y = 2

position = zeros([2, nParticles])            # Material particles positions
increments = zeros([2, nParticles])     # Material particles increments
energies = zeros([1, iterations + 1])   # Energy values will be stored here when computed


#---------------------------------------- AUXILIAR FUNCTIONS ----------------------------------------#


def particleInteractionEnergy(pos, inc, i, withIncrements = False):  # Computes the total interaction energy of particle i
    # When increments = True, it considers de future positions
    Eij = []
    if (not withIncrements):
        ri = [pos[0][i], pos[1][i]]
        for j in range(nParticles):
            if (i != j):
                rj = [pos[0][j], pos[1][j]]
                r = sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2)
                Eij.append(interactionType1(r))
    else:
        ri = [pos[0][i] + inc[0][i], pos[1][i] + inc[1][i]]
        for j in range(nParticles):
            if (i != j):
                rj = [pos[0][j] + inc[0][j], pos[1][j] + inc[1][j]]
                r = sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2)
                Eij.append(interactionType1(r))
    return sum(Eij)


def sortModifyIncrement1(inc, i):   # Sorts and modify de increment of particle i
    angle = 2 * pi * random()
    dr = r01/5 * random()
    inc[0][i] = dr * cos(angle)
    inc[1][i] = dr * sin(angle)


def modifyPos(pos, inc, i): # Modifies the position of particle i
    pos[0][i] += inc[0][i]
    pos[1][i] += inc[1][i]


def boundryControl(pos, inc, i):    # Controls the particles do not get out of the boundries
    mod = False     # If mod = True it means some value has been modified
    x = pos[0][i]
    dx = inc[0][i]
    y = pos[1][i]
    dy = inc[1][i]
    if (x + dx < x1):
        dx = abs(x1 - x)
        mod = True
    elif (x + dx > x2):
        dx = abs(x2 - x)
        mod = True
    if (y + dy < y1):
        dy = abs(y1 - y)
        mod = True
    elif (y + dy > y2):
        dy = abs(y2 - y)
        mod = True
    if mod:
        inc[0][i] = dx
        inc[1][i] = dy


#---------------------------------------- SIMULATION ----------------------------------------#


for i in range(nParticles):     # Starting positions of the material particles
    position[0][i] = (2 * random() - 1) / 20
    position[1][i] = random() / 10

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
xlim([-2.1, 2.1])
ylim([-0.1, 2.1])

for i in range(iterations):     # Computes the changes in the system and shows the simulation
    index = randrange(0, nParticles)
    sortModifyIncrement1(increments, index)
    boundryControl(position, increments, index)
    Eactual = particleInteractionEnergy(position, increments, index)      # Actual interaction energy of chosen particle
    Emod = particleInteractionEnergy(position, increments, index, True)   # Future interaction energy of chosen particle
    dE = Emod - Eactual
    if (acceptance(dE)):
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
