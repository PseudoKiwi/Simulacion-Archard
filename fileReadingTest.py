from numpy import pi, sin, cos, sqrt, zeros
from probability import defineProbability
from potential import potential
from random import random, randrange
from matplotlib.pyplot import xlim, ylim, figure, show


#---------------------------------------- CONSTANTS DEFINITION ----------------------------------------#


iterations = 1000000   # Number of Monte Carlo iterations
nParticles = 100    # Number of particles on the material
U01 = 1             # Lennard - Jones potential constant from material
r01 = 0.3           # Equilibrium radius from material
beta = 100          # Related to temperature constant
accepted = [0]      # Will store the amount of times a change was accepted when dE > 0

acceptance1 = defineProbability(accepted, beta)        # acceptance is the probability function for change acceptance
interactionType1 = potential(U01, r01)                       # Interaction potential between material particles

x1 = 0      # wall at x = 0
x2 = 2      # wall at x = 2
y1 = 0      # wall at y = 0
y2 = 2      # wall at y = 2

position = zeros([2, nParticles])            # Material particles positions
increments = zeros([2, nParticles])     # Material particles increments
energies = zeros([1, iterations + 1])   # Energy values will be stored here when computed
dEnergies = zeros([1, iterations])


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


def modifyPos(pos, inc, i):  # Modifies the position of particle i
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


with open("ultimoDato.txt", "r") as file:
    a = file.readline()
    b = file.readline()

    a = a.split(',')
    a[0] = a[0].split('[')[1]
    a[-1] = a[-1].split(']')[0]

    b = b.split(',')
    b[0] = b[0].split('[')[1]
    b[-1] = b[-1].split(']')[0]

    position[0] = a
    position[1] = b

E = 0
for i in range(nParticles):
    E += particleInteractionEnergy(position, increments, i)

energies[0] = E

for i in range(iterations):     # Computes the changes in the system and shows the simulation
    index = randrange(0, nParticles)
    sortModifyIncrement1(increments, index)
    boundryControl(position, increments, index)
    Eactual = particleInteractionEnergy(position, increments, index)      # Actual interaction energy of chosen particle
    Emod = particleInteractionEnergy(position, increments, index, True)   # Future interaction energy of chosen particle
    dE = Emod - Eactual
    if (acceptance1(dE)):
        E += dE
        modifyPos(position, increments, index)
        dEnergies[0][i] = dE
    energies[0][i+1] = E


X = position[0]  # Every x position
Y = position[1]  # Every y position
fig = figure()
ax = fig.add_subplot(111)
system, = ax.plot(X, Y, ".")
xlim([-0.1, 2.1])
ylim([-0.1, 2.1])

fig2 = figure()
ax2 = fig2.add_subplot(111)
ax2.plot(energies[0])
show()

fig3 = figure()
ax3 = fig3.add_subplot(111)
ax3.plot(dEnergies[0])
show()

with open("ultimoDato.txt", "w") as file:
    str1 = str(list(X))
    str2 = str(list(Y))
    file.write(str1 + "\n")
    file.write(str2 + "\n")

print(accepted[0])
AR = accepted[0] / iterations * 100
print(AR)
