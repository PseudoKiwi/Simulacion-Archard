from numpy import zeros, abs, mean
from random import randrange
from matplotlib.pyplot import xlim, ylim, figure, ion, show
import auxiliarFunctions as auxF


#---------------------------------------- CONSTANTS DEFINITION ----------------------------------------#


iterations = 10000   # Number of Monte Carlo iterations
nParticles = 100    # Number of particles on the material
U01 = 1             # Lennard - Jones potential constant from material
r01 = 0.5           # Equilibrium radius from material
beta = 100          # Related to temperature constant
accepted = [0]      # Will store the amount of times a change was accepted when dE > 0

acceptance1 = auxF.defineProbability(accepted, beta)        # acceptance is the probability function for change acceptance
interactionType1 = auxF.potential(U01, r01)                       # Interaction potential between material particles
vForce = auxF.verticalInteractionForce(U01, r01)

x1 = -2     # wall at x = 0
x2 = 4      # wall at x = 2
y1 = 0      # wall at y = 0
y2 = 5      # wall at y = 5

yi = y2
yf = y2 - 3*r01

ly = y2 - y1
lx = x2 - x1
dy = ly/100
dx = dy*lx/(2*(ly - dy))
expData = int(abs(yi - yf)/dy)       # Number of total iterations of the simulation

boundryP = int((x2 - x1)//r01) + 1
position = zeros([2, nParticles + boundryP])        # Material particles positions
increments = zeros([2, nParticles + boundryP])      # Material particles increments
energies = zeros([1, iterations + 1])               # Energy values will be stored here when computed
pressures = zeros([1, iterations])                  # Pressure values will be stored here when computed
eqEnergies = zeros([1, expData])
eqPressures = zeros([1, expData])


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

    position[0][boundryP:] = a
    position[1][boundryP:] = b

for i in range(boundryP):
    position[0][i] = r01 * i + ((x2 - x1) % r01) + x1
    position[1][i] = y2

E = auxF.interactionEnergy(position, nParticles + boundryP, interactionType1)
energies[0] = E

X = position[0]  # Every x position
Y = position[1]  # Every y position
fig1 = figure()
ax1 = fig1.add_subplot(111)
system, = ax1.plot(X, Y, ".")
xlim([-3, 5])
ylim([-0.1, 5])

for h in range(expData):
    print(h + 1)
    for i in range(iterations):     # Computes the changes in the system and shows the simulation
        index = randrange(0, nParticles + boundryP)
        auxF.sortModifyIncrement(increments, index, r01, boundryP)
        auxF.boundryControl(position, increments, index, x1, x2, y1, y2)

        dE = auxF.dEi(position, increments, index, nParticles + boundryP, interactionType1)

        if (acceptance1(dE)):
            auxF.modifyPos(position, increments, index)
        else:
            dE = 0
        energies[0][i+1] = energies[0][i] + dE

        vF = zeros(boundryP)
        if (i >= iterations - int(iterations/100)):
            for f in range(boundryP):
                vFie = []
                for e in range(boundryP, nParticles + boundryP):
                    vFie.append(vForce(position[0][f], position[0][e], position[1][f], position[1][e]))
                vF[f] = sum(vFie)
            totalForce = sum(vF)
            pressure = totalForce/lx

            pressures[0][i] = pressure

        if (i % 9999 == 0):
            print(i+1)

    y2 -= dy
    x1 -= dx
    x2 += dx
    auxF.roofParticlesMovement(position, boundryP, -dy)
    auxF.totalBoundryControl(position, x1, x2, y1, y2)

    eqEnergies[0][h] = mean(energies[0][-int(iterations/100):])
    eqPressures[0][h] = mean(pressures[0][-int(iterations/100):])


X = position[0]  # Every x position
Y = position[1]  # Every y position
fig2 = figure()
ax2 = fig2.add_subplot(111)
system, = ax2.plot(X, Y, ".")
xlim([-3, 5])
ylim([-0.1, 5])

fig3 = figure()
ax3 = fig3.add_subplot(111)
ax3.plot(eqEnergies[0])

fig4 = figure()
ax4 = fig4.add_subplot(111)
ax4.plot(eqPressures[0])
show()