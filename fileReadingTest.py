from numpy import zeros
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

x1 = -2     # wall at x = 0
x2 = 4      # wall at x = 2
y1 = 0      # wall at y = 0
y2 = 6      # wall at y = 5

position = zeros([2, nParticles])            # Material particles positions
increments = zeros([2, nParticles])     # Material particles increments
energies = zeros([1, iterations + 1])   # Energy values will be stored here when computed


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

E = auxF.interactionEnergy(position, nParticles, interactionType1)
energies[0] = E

for i in range(iterations):     # Computes the changes in the system and shows the simulation
    index = randrange(0, nParticles)
    auxF.sortModifyIncrement(increments, index, r01)
    auxF.boundryControl(position, increments, index, x1, x2, y1, y2)

    dE = auxF.dEi(position, increments, index, nParticles, interactionType1)

    if (acceptance1(dE)):
        auxF.modifyPos(position, increments, index)
    else:
        dE = 0
    energies[0][i+1] = energies[0][i] + dE


X = position[0]  # Every x position
Y = position[1]  # Every y position
fig = figure()
ax = fig.add_subplot(111)
system, = ax.plot(X, Y, ".")
xlim([-3, 5])
ylim([-0.1, 5])

fig2 = figure()
ax2 = fig2.add_subplot(111)
ax2.plot(energies[0])
show()

#with open("ultimoDato.txt", "w") as file:
#    str1 = str(list(X))
#    str2 = str(list(Y))
#    file.write(str1 + "\n")
#    file.write(str2 + "\n")

print(accepted[0])
AR = accepted[0] / iterations * 100
print(AR)