from numpy import zeros, sqrt, inf
from random import randrange
from matplotlib.pyplot import xlim, ylim, figure, show
import auxiliarFunctions as auxF

# ---------------------------------------- CONSTANTS DEFINITION ----------------------------------------#
# ---------------------------------------------- GENERAL -----------------------------------------------#


iterations = 1000000  # Number of Monte Carlo iterations
mParticles = 99  # Number of particles on the material
wBase = 8   # Number of particles forming the wedge base
wParticles = wBase*(wBase + 1)//2   # Number of wedge particles
a = 17  # a is chosen such that a divides mParticles, so there are an integer number of "levels"
b = mParticles // (2 * a - 1)  # Number of material "levels"
fParticles = 2*b + a   # Number of fixed material particles
nParticles = mParticles + wParticles + fParticles
beta = 100  # Related to temperature constant

with open("r01.txt", "r") as file:
    aux = float(file.readline())
    req1 = aux
    req2 = aux
    reqInt = aux

r0Int = 0.5     # Interaction radius material - wedge
U0Int = 3       # Lennard - Jones potential constant from interaction

position = zeros([2, nParticles])    # Particles positions
increments = zeros([2, nParticles])  # Particles increments
energies = zeros([iterations + 1])   # Energy values will be stored here when computed

# x1 = 0  # wall at x = 0
# x2 = 16 * req1  # wall at x = 2
# y1 = 0  # wall at y = 0
# y2 = inf  # wall at y = H


# ---------------------------------------------- MATERIAL ----------------------------------------------#


h1 = sqrt(3) / 2 * req1     # Equilibrium material triangle height
H = h1 * (2 * b - 1)        # Total height of material structure simulated

r01 = 0.5  # Interaction radius from material
U01 = 1    # Lennard - Jones potential constant from material


# ----------------------------------------------- WEDGE -----------------------------------------------#


h2 = sqrt(3) / 2 * req2

r02 = 0.5    # Interaction radius for wedge particles
U02 = 5      # Lennard - Jones potential constant from wedge


#--------------------------------------------- SIMULATION ---------------------------------------------#


# Starting positions of the wedge particles (triangle form)
dx = (a - 1)*req1/2 - (wBase - 1)*req2/2
dy = H + 2*reqInt
aux = wBase
e = 0
f = 1
mod = False
for i in range(wParticles):
    position[0][i] = (f % aux)*req2 + e*req2/2 + dx
    position[1][i] = (wBase-1) * h2 - e*h2 + dy
    if f % aux == 0 and f // aux == 1:
        aux -= 1
        e += 1
        mod = True
    if mod:
        f = 1
        mod = False
    else:
        f += 1


# Starting positions of the material particles
for i in range(a):
    position[0][i + wParticles] = req1*i
    position[1][i + wParticles] = H

for i in range(a - 1):
    position[0][i + a + wParticles] = req1*i + req1/2
    position[1][i + a + wParticles] = H - h1

for i in range(mParticles - 2*a + 1):
    position[0][i + 2*a - 1 + wParticles] = position[0][i + wParticles]
    position[1][i + 2*a - 1 + wParticles] = position[1][i + wParticles] - 2*h1


# Starting positions of fixed material particles
for i in range(b):
    position[0][i + wParticles + mParticles] = position[0][wParticles + a] - req1
    position[1][i + wParticles + mParticles] = position[1][a * (i + 1) + i * (a - 1) + wParticles]
    position[0][i + wParticles + mParticles + b] = position[0][wParticles + 2*a - 2] + req1
    position[1][i + wParticles + mParticles + b] = position[1][a * (i + 1) + i * (a - 1) + wParticles]

for i in range(a):
    position[0][i + wParticles + mParticles + 2*b] = position[0][wParticles + i]
    position[1][i + wParticles + mParticles + 2*b] = position[1][wParticles + i] - 2*b*h1


# Initial positions graph
X1 = position[0][wParticles:]  # Every x material position
Y1 = position[1][wParticles:]  # Every y material position
X2 = position[0][-fParticles:]  # Every x fixed position
Y2 = position[1][-fParticles:]  # Every y fixed position
X3 = position[0][:wParticles]  # Every x wedge position
Y3 = position[1][:wParticles]  # Every x wedge position
fig1 = figure()
ax2 = fig1.add_subplot(111)
system, = ax2.plot(X1, Y1, ".")
system, = ax2.plot(X2, Y2, ".")
system, = ax2.plot(X3, Y3, ".")
xlim([-req1, req1*a])
ylim([-2*req1, H + wBase*h2 + 1.5])
show()

# Initial system energy
E = auxF.interactionEnergy2(position, nParticles, wParticles, U01, r01, U02, r02, U0Int, r0Int)
energies[0] = E

partial = 12
force = zeros([3*partial + 1])
F = 0
for k in range(wParticles):
    for w in range(mParticles):
        F += auxF.totalForce(U0Int, r0Int, position[0][k], position[0][w + wParticles], position[1][k], position[1][w + wParticles])
force[0] = F

for e in range(3*partial):

    # Conditional to move down or up de wedge
    if e < partial:
        for i in range(wBase):
            position[1][i] -= req2 / 4
    else:
        for i in range(wBase):
            position[1][i] += req2 / 4


    print(e)    # Control print

    for i in range(iterations):     # Computes the changes in the system and shows the simulation

        if i % 50000 == 0:
            print(i)    # Control print

        index = randrange(wBase, nParticles - fParticles)   # Sort particle index
        if index >= wParticles:
            auxF.sortModifyIncrement(increments, index, r01)    # Material particle index action
        else:
            auxF.sortModifyIncrement(increments, index, r02)
        # auxF.boundryControl(position, increments, index, x1, x2, y1, y2)    # Wedge particle index action

        # Energy change caused by possible movement
        dE = auxF.dEi2(position, increments, index, nParticles, wParticles, U01, r01, U02, r02, U0Int, r0Int)

        # The probability function result gives the acceptance of the movement
        if auxF.probability(beta, dE):
            auxF.modifyPos(position, increments, index)
        else:
            dE = 0
        energies[i + 1] = energies[i] + dE

    F = 0
    for k in range(wParticles):
        for w in range(mParticles):
            F += auxF.totalForce(U0Int, r0Int, position[0][k], position[0][w + wParticles], position[1][k],
                                 position[1][w + wParticles])
    force[e + 1] = F

    #   Control graph
    X1 = position[0][wParticles:]  # Every x material position
    Y1 = position[1][wParticles:]  # Every y material position
    X2 = position[0][-fParticles:]  # Every x fixed position
    Y2 = position[1][-fParticles:]  # Every y fixed position
    X3 = position[0][:wParticles]  # Every x wedge position
    Y3 = position[1][:wParticles]  # Every x wedge position
    fig1 = figure()
    ax1 = fig1.add_subplot(111)
    system, = ax1.plot(X1, Y1, ".")
    system, = ax1.plot(X2, Y2, ".")
    system, = ax1.plot(X3, Y3, ".")
    xlim([-req1, req1 * a])
    ylim([-2 * req1, H + wBase * h2 + 1.5])
    show()

loss = 0
for i in range(mParticles):
    if position[1][i + wParticles] >= H + reqInt:
        loss += 1

fig2 = figure()
ax2 = fig2.add_subplot(111)
system, = ax2.plot(force)
show()
print(max(force))
print(loss)

with open("force-loss.txt", "a") as file:
    file.write(f"{max(force)}, {loss}\n")