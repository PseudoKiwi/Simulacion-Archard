from numpy import pi, sin, cos, sqrt, e
from random import random


#----------------------------------------- AUXILIAR FUNCTIONS -----------------------------------------#


def interactionEnergy(pos, nParticles, interaction):  # Computes the total interaction energy of particle i
    # When increments = True, it considers de future positions
    Eij = []
    for i in range(nParticles):
        ri = [pos[0][i], pos[1][i]]
        for j in range(nParticles):
            if (i != j):
                rj = [pos[0][j], pos[1][j]]
                r = sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2)
                Eij.append(interaction(r))
    return sum(Eij)/2




def dEi(pos, inc, i, nParticles, interaction):
    dEi = []
    for j in range(nParticles):
        if (i != j):
            r1 = sqrt((pos[0][i] - pos[0][j]) ** 2 + (pos[1][i] - pos[1][j]) ** 2)
            r2 = sqrt((pos[0][i] + inc[0][i] - pos[0][j]) ** 2 + (pos[1][i] + inc[1][i] - pos[1][j]) ** 2)
            dEi.append(interaction(r2) - interaction(r1))
    return sum(dEi)

def sortModifyIncrement1(inc, i, r0):   # Sorts and modify de increment of particle i
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


def roofMovement(pos, n, distance):
    for i in range(n):
        pos[1][i] += distance


# The defineProbability function receives the parameter beta, such that 1/beta=kB*T. It then returns
# the probability function to be used during the simulation asuming negligible changes in temperature.
def defineProbability(accepted, beta):
    def acceptance(dE):
        if( dE <= 0 ):
            accepted[0] += 1
            return True
        else:
            r = random()
            aux = r < e**(-beta*dE)
            if aux:
                accepted[0] += 1
            return aux
    return acceptance


# The potential function receives the parameters of the Lennard-Jones
# potential energy and returns the function to aply with only the
# relative distances between particles.

# U0: Measures how strong particles atract each other
# r0: Equilibrium radius

def potential(U0, r0):
    def energy(r):
        return U0*((r0/r)**12 - (r0/r)**6)
    return energy


# Vertical component of the interaction force between 2 particles
# Used to compute pressure over the top or botton boundries
def verticalInteractionForce(U0, r0):
    def vForce(x, x0, y, y0):
        r = sqrt((x-x0)**2 + (y-y0)**2)
        return U0*(12*(r0**12)/(r**14) - 6*(r0**6)/(r**8))*(y-y0)
    return vForce


# Horizontal component of the interaction force between 2 particles
# Used to compute pressure over the left or right boundries
def horizontalInteractionForce(U0, r0):
    def hForce(x, x0, y, y0):
        r = sqrt((x - x0) ** 2 + (y - y0) ** 2)
        return U0*(12*(r0**12)/(r**14) - 6*(r0**6)/(r**8))*(x-x0)
    return hForce