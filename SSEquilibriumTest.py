from numpy import zeros, mean, inf
from random import randrange
from matplotlib.pyplot import xlim, ylim, figure, ion, show
import auxiliarFunctions as auxF


#---------------------------------------- CONSTANTS DEFINITION ----------------------------------------#


iterations = 2000000  # Number of Monte Carlo iterations
nParticles = 100      # Number of particles on the material
U01 = 1               # Lennard - Jones potential constant from material
r01 = 0.5             # Equilibrium radius from material
beta = 100            # Related to temperature constant

x1 = -2     # wall at x = 0
x2 = 4      # wall at x = 2
y1 = 0      # wall at y = 0
y2 = inf      # wall at y = inf

position = zeros([2, nParticles])            # Material particles positions
increments = zeros([2, nParticles])     # Material particles increments
energies = zeros([iterations + 1])   # Energy values will be stored here when computed
auxiliar = []

#---------------------------------------- SIMULATION ----------------------------------------#


aux = 2 + 2/9
for i in range(nParticles):     # Starting positions of the material particles
    position[0][i] = 2/9 * (i % 10)
    if (i % 10 == 0):
        aux -= 2/9
    position[1][i] = aux

E = auxF.interactionEnergy(position, nParticles, U01, r01)
energies[0] = E

X = position[0]  # Every x position
Y = position[1]  # Every y position

#ion()
fig = figure()
ax = fig.add_subplot(111)
system, = ax.plot(X, Y, ".")
xlim([-3, 5])
ylim([-0.1, 5])
show()

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

    if (i > 700000):
        auxiliar.append(energies[i+1])

    if (i % 10000 == 0):
        print(i)

    #    X = position[0]  # Every x position
    #    Y = position[1]  # Every y position
    #system.set_xdata(X)
    #system.set_ydata(Y)
    #fig.canvas.draw()
    #fig.canvas.flush_events()


X = position[0]  # Every x position
Y = position[1]  # Every y position
fig = figure()
ax = fig.add_subplot(111)
system, = ax.plot(X, Y, ".")
xlim([-3, 5])
ylim([-0.1, 5])

fig2 = figure()
ax2 = fig2.add_subplot(111)
ax2.plot(energies)
show()

with open("materialEq.txt", "w") as file:
    str1 = str(list(X))
    str2 = str(list(Y))
    str3 = str(list(energies))
    file.write(str1 + "\n")
    file.write(str2 + "\n")
    file.write(str3 + "\n")

print(mean(auxiliar))
fig3 = figure()
ax3 = fig3.add_subplot(111)
ax3.plot(auxiliar)
show()