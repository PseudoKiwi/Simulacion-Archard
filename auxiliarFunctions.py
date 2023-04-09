from numpy import pi, sin, cos, sqrt, e
from random import random


#----------------------------------------- AUXILIAR FUNCTIONS -----------------------------------------#


def interactionEnergy(pos, nParticles, U0, r0):
    Eij = []
    for i in range(nParticles):
        ri = [pos[0][i], pos[1][i]]
        for j in range(nParticles):
            if (i != j):
                rj = [pos[0][j], pos[1][j]]
                r = sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2)
                Eij.append(potential(U0, r0, r))
    return sum(Eij)/2

def interactionEnergy2(pos, mParticles, wParticles, U01, r01, U02, r02, U0Int, r0Int):
    Eij = []
    for i in range(mParticles + wParticles):
        ri = [pos[0][i], pos[1][i]]
        for j in range(mParticles + wParticles):
            if (i != j):
                rj = [pos[0][j], pos[1][j]]
                r = sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2)

                if i >= wParticles and j >= wParticles:
                    Eij.append(potential(U01, r01, r))
                elif (i >= wParticles > j) or (i < wParticles <= j):
                    Eij.append(potential(U0Int, r0Int, r))
                else:
                    Eij.append(potential(U02, r02, r))

    return sum(Eij)/2

def dEi(pos, inc, i, nParticles, U0, r0):
    dEi = []
    for j in range(nParticles):
        if (i != j):
            r1 = sqrt((pos[0][i] - pos[0][j]) ** 2 + (pos[1][i] - pos[1][j]) ** 2)
            r2 = sqrt((pos[0][i] + inc[0][i] - pos[0][j]) ** 2 + (pos[1][i] + inc[1][i] - pos[1][j]) ** 2)
            dEi.append(potential(U0, r0, r2) - potential(U0, r0, r1))
    return sum(dEi)

def dEi2(pos, inc, i, mParticles, wParticles, U01, r01, U02, r02, U0Int, r0Int):
    dEi = []
    for j in range(mParticles + wParticles):
        if (i != j):
            r1 = sqrt((pos[0][i] - pos[0][j]) ** 2 + (pos[1][i] - pos[1][j]) ** 2)
            r2 = sqrt((pos[0][i] + inc[0][i] - pos[0][j]) ** 2 + (pos[1][i] + inc[1][i] - pos[1][j]) ** 2)

            if i >= wParticles and j >= wParticles:
                dEi.append(potential(U01, r01, r2) - potential(U01, r01, r1))
            elif (i >= wParticles > j) or (i < wParticles <= j):
                dEi.append(potential(U0Int, r0Int, r2) - potential(U0Int, r0Int, r1))
            else:
                dEi.append(potential(U02, r02, r2) - potential(U02, r02, r1))

    return sum(dEi)

def sortModifyIncrement(inc, i, r0):   # Sorts and modify de increment of particle i
    angle = 2 * pi * random()
    dr = r0/5 * random()
    inc[0][i] = dr * cos(angle)
    inc[1][i] = dr * sin(angle)


def modifyPos(pos, inc, i): # Modifies the position of particle i
    pos[0][i] += inc[0][i]
    pos[1][i] += inc[1][i]


def boundryControl(pos, inc, i, x1, x2, y1, y2):    # Controls the particles do not get out of the boundries
    mod = False     # If mod = True it means some value has been modified
    x = pos[0][i]
    dx = inc[0][i]
    y = pos[1][i]
    dy = inc[1][i]
    if (x + dx < x1):
        dx = x1 - x
        mod = True
    elif (x + dx > x2):
        dx = x - x2
        mod = True
    if (y + dy < y1):
        dy = y1 - y
        mod = True
    elif (y + dy > y2):
        dy = y - y2
        mod = True
    if mod:
        inc[0][i] = dx
        inc[1][i] = dy


def totalBoundryControl(pos, x1, x2, y1, y2):
    for i, p in enumerate(pos[0]):
        if (pos[0][i] < x1):
            pos[0][i] = x1
        elif (pos[0][i] > x2):
            pos[0][i] = x2
        if (pos[1][i] < y1):
            pos[1][i] = y1
        elif (pos[1][i] > y2):
            pos[1][i] = y2


# The defineProbability function receives the parameter beta, such that 1/beta=kB*T. It then returns
# the probability function to be used during the simulation asuming negligible changes in temperature.
def probability(beta, dE):
    if( dE <= 0 ):
        return True
    else:
        r = random()
        aux = r < e**(-beta*dE)
        return aux


# The potential function receives the parameters of the Lennard-Jones
# potential energy and returns the function to aply with only the
# relative distances between particles.

# U0: Measures how strong particles atract each other
# r0: Equilibrium radius

def potential(U0, r0, r):
    return U0*((r0/r)**12 - (r0/r)**6)


# Vertical component of the interaction force between 2 particles
# Used to compute pressure over the top or botton boundries
def verticalInteractionForce(U0, r0, x, x0, y, y0):
    r = sqrt((x-x0)**2 + (y-y0)**2)
    return U0*(12*(r0**12)/(r**14) - 6*(r0**6)/(r**8))*(y-y0)


# Horizontal component of the interaction force between 2 particles
# Used to compute pressure over the left or right boundries
def horizontalInteractionForce(U0, r0, x, x0, y, y0):
    r = sqrt((x - x0) ** 2 + (y - y0) ** 2)
    return U0*(12*(r0**12)/(r**14) - 6*(r0**6)/(r**8))*(x-x0)
