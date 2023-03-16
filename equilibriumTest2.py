from numpy import array, pi, sin, cos, sqrt, inf, abs
from probability import defineProbability
from potential import potential
from materialParticle import MaterialParticle
from random import random, randrange
from matplotlib.pyplot import xlim, ylim, figure, ion, show

iterations = 1000   # Number of Monte Carlo iterations
nParticles = 100    # Number of particles on the material
U01 = 1             # Lennard - Jones potential constant from material
r01 = 0.5           # Equilibrium radius from material
beta = 100          # Related to temperature constant
material = []       # Material particles will be stored here when created
Energies = []       # Energy values will be stored here when computed
dENeg = [0]         # Will store the amount of times dE < 0
accepted = [0]      # Will store the amount of times a change was accepted when dE > 0

acceptance = defineProbability(dENeg, accepted, beta)        # acceptance is the probability function for change acceptance
interactionType1 = potential(U01, r01)      # Interaction potential between material particles

x1 = -2  # wall at x = -2
x2 = 2   # wall at x = 2
y1 = 0   # wall at y = 0
y2 = 2   # wall at y = inf


def getX(p):  # Returs de x position of a particle p
    return p.getPos()[0]


def getY(p):  # Returs de y position of a particle p
    return p.getPos()[1]


def sortModifyIncrement1(p):  # Sorts and modify de increment of particle p
    angle = 2 * pi * random()
    dr = r01 / 5 * random()
    p.setIncrement(array([dr * cos(angle), dr * sin(angle)], float))


def modifyPos(p):  # Modifies the position of particle p
    p.modifyPos()


def systemEnergy(material, increments=False):  # Computes the total interaction energy of the system
                                               # When increments = True, it considers de future positions
    Eij = []
    if (not increments):
        for i, p in enumerate(material):
            for j in range(i + 1, len(material)):
                if (i != j):
                    ri = p.getPos()
                    rj = material[j].getPos()
                    r = sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2)
                    Eij.append(interactionType1(r))
    else:
        for i, p in enumerate(material):
            for j in range(nParticles):
                if (i != j):
                    ri = p.getPos() + p.getIncrement()
                    rj = material[j].getPos() + material[j].getIncrement()
                    r = sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2)
                    Eij.append(interactionType1(r))
    return sum(Eij)


def boundryControl(p):  # Controls the particles do not get out of the boundries
    mod = False         # If mod=True it means some value has been modified
    pos = p.getPos()
    inc = p.getIncrement()
    x = pos[0]
    dx = inc[0]
    y = pos[1]
    dy = inc[1]
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
        p.setIncrement(array([dx, dy], float))  # Sets the new increments of the particle p


for i in range(nParticles):     # Creation of the material particles
    newP = MaterialParticle(array([(2 * random() - 1) / 20, random() / 10]), 1)
    material.append(newP)

E = systemEnergy(material)
Energies.append(E)

X = list(map(getX, material))
Y = list(map(getY, material))

ion()
fig = figure()
ax = fig.add_subplot(111)
system, = ax.plot(X, Y, ".")
xlim([-2.1, 2.1])
ylim([-0.1, 2])

for i in range(iterations):     # Computes the changes in the system and shows the simulation
    index = randrange(0, nParticles-1)
    p = material[index]
    sortModifyIncrement1(p)
    boundryControl(p)
    Emod = systemEnergy(material, True)
    dE = Emod - E
    if (acceptance(dE)):
        E += dE
        modifyPos(p)
        X = list(map(getX, material))
        Y = list(map(getY, material))
    Energies.append(E)
    system.set_xdata(X)
    system.set_ydata(Y)
    fig.canvas.draw()
    fig.canvas.flush_events()

fig2 = figure()
ax2 = fig2.add_subplot(111)
ax2.plot(Energies)
show()

print(dENeg[0])
print(accepted[0])
AR = (accepted[0] / (iterations - dENeg[0]))*100
print(AR)