import numpy as np
from libcpp cimport bool
from random import randrange, random
from matplotlib.pyplot import xlim, ylim, figure, show


#----------------------------------------- USED FUNCTIONS -----------------------------------------#


cpdef float interactionEnergy(X,  Y, int nParticles, float U0, float r0):  # Computes the total interaction energy of particle i
    # When increments = True, it considers de future positions
    cdef float Eij = 0
    cdef float ri[2]
    cdef float rj[2]
    cdef float r
    cdef int i = 0
    cdef int j = 0
    cdef int cont = 0

    while i < nParticles:
        ri[0] = X[i]
        ri[1] = Y[i]
        while i < nParticles:
            if (i != j):
                rj[0] = X[j]
                rj[1] = Y[j]
                r = np.sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2)
                Eij += potential(U0, r0, r)
            j += 1
            cont += 1
        i += 1
    return Eij/2


cpdef float dEi(X,  Y,  incX,  incY, int i, int nParticles, float U0, float r0):
    cdef float dEi = 0
    cdef int j = 0
    cdef float r1
    cdef float r2

    while j < nParticles:
        if (i != j):
            r1 = np.sqrt((X[i] - X[j]) ** 2 + (Y[i] - Y[j]) ** 2)
            r2 = np.sqrt((X[i] + incX[i] - X[j]) ** 2 + (Y[i] + incY[i] - Y[j]) ** 2)
            dEi += potential(U0, r0, r2) - potential(U0, r0, r1)
        j += 1
    return dEi


cpdef void sortModifyIncrement( incX,  incY, int i, float r0, int boundryP):   # Sorts and modify de increment of particle i
    cdef float angle = 2 * np.pi * random()
    cdef float dr = r0/5 * random()
    incX[i] = dr * np.cos(angle)
    if (boundryP == 0):
        incY[i] = dr * np.sin(angle)
    else:
        if (i < boundryP):
            incY[i] = 0


cpdef void modifyPos( X,  Y,  incX,  incY, int i): # Modifies the position of particle i
    X[i] += incX[i]
    Y[i] += incY[i]


cpdef void boundryControl( X,  incX,  Y,  incY, int i, float x1, float x2, float y1, float y2):    # Controls the particles do not get out of the boundries
    cdef bool mod = False     # If mod = True it means some value has been modified
    cdef float x = X[i]
    cdef float dx = incX[i]
    cdef float y = Y[i]
    cdef float dy = incY[i]

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
        incX[i] = dx
        incY[i] = dy


cpdef void totalBoundryControl( X,  Y, float x1, float x2, float y1, float y2, int nParticles):
    cdef int i = 0

    while i < nParticles:
        if (X[i] < x1):
            X[i] = x1
        elif (X[i] > x2):
            X[i] = x2
        if (Y[i] < y1):
            Y[i] = y1
        elif (Y[i] > y2):
            Y[i] = y2
        i += 1


cpdef void roofParticlesMovement( Y, int boundryP, float distance):
    cdef int i = 0

    while i < boundryP:
        Y[i] += distance
        i += 1


# The defineProbability function receives the parameter beta, such that 1/beta=kB*T. It then returns
# the probability function to be used during the simulation asuming negligible changes in temperature.
cpdef bool defineProbability(accepted, float beta, float dE):
    cdef float r
    cdef bool aux

    if( dE <= 0 ):
        accepted[0] += 1
        return True
    else:
        r = random()
        aux = r < np.e**(-beta*dE)
        if aux:
            accepted[0] += 1
        return aux


# The potential function receives the parameters of the Lennard-Jones
# potential energy and returns the function to aply with only the
# relative distances between particles.

# U0: Measures how strong particles atract each other
# r0: Equilibrium radius

cpdef float potential(float U0, float r0, float r):
    return U0*((r0/r)**12 - (r0/r)**6)


# Vertical component of the interaction force between 2 particles
# Used to compute pressure over the top or botton boundries
cpdef float verticalInteractionForce(float U0, float r0, float x, float x0, float y, float y0):
    cdef float r = np.sqrt((x-x0)**2 + (y-y0)**2)
    return U0*(12*(r0**12)/(r**14) - 6*(r0**6)/(r**8))*(y-y0)


# Horizontal component of the interaction force between 2 particles
# Used to compute pressure over the left or right boundries
cpdef float horizontalInteractionForce(float U0, float r0, float x, float x0, float y, float y0):
    cdef float r = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
    return U0*(12*(r0**12)/(r**14) - 6*(r0**6)/(r**8))*(x-x0)


cpdef void simulation():
    #---------------------------------------- CONSTANTS DEFINITION ----------------------------------------#

    DEF iterations = 1000000  # Number of Monte Carlo iterations
    DEF nParticles = 100  # Number of particles on the material
    cdef float U01 = 1  # Lennard - Jones potential constant from material
    cdef float r01 = 0.5  # Equilibrium radius from material
    cdef float beta = 100  # Related to temperature constant
    cdef int accepted[1] # Will store the amount of times a change was accepted when dE > 0
    accepted[0] = 0

    cdef float x1 = -2  # wall at x = 0
    cdef float x2 = 4  # wall at x = 2
    cdef float y1 = 0  # wall at y = 0
    cdef float y2 = 5  # wall at y = 5

    cdef float X[nParticles]      # Material particles positions
    cdef float Y[nParticles]
    cdef float incrementsX[nParticles]    # Material particles increments
    cdef float incrementsY[nParticles]
    cdef float energies[iterations + 1]     # Energy values will be stored here when computed
    cdef int auxiliar[200000]

    #---------------------------------------- SIMULATION ----------------------------------------#

    cdef float aux = 2 + 2 / 9
    cdef int i = 0
    while i < nParticles:  # Starting positions of the material particles
        X[i] = 2 / 9 * (i % 10)
        if (i % 10 == 0):
            aux -= 2 / 9
        Y[i] = aux
        i += 1

    cdef float E = interactionEnergy(X, Y, nParticles, U01, r01)
    energies[0] = E

    #ion()
    fig = figure()
    ax = fig.add_subplot(111)
    system, = ax.plot(X, Y, ".")
    xlim([-3, 5])
    ylim([-0.1, 5])
    show()

    cdef int index

    i = 0
    while i < iterations:  # Computes the changes in the system and shows the simulation
        index = randrange(0, nParticles)
        sortModifyIncrement(incrementsX, incrementsY, index, r01, 0)
        boundryControl(X, Y, incrementsX, incrementsY, index, x1, x2, y1, y2)

        dE = dEi(X, Y, incrementsX, incrementsY, index, nParticles, U01, r01)

        if (defineProbability(accepted, beta, dE)):
            modifyPos(X, Y, incrementsX, incrementsY, index)
        else:
            dE = 0
        energies[i + 1] = energies[i] + dE

        if (i > 800000):
            auxiliar.append(energies[i + 1])

        if (i % 10000 == 0):
            print(i)

        i += 1

        #    X = position[0]  # Every x position
        #    Y = position[1]  # Every y position
        #system.set_xdata(X)
        #system.set_ydata(Y)
        #fig.canvas.draw()
        #fig.canvas.flush_events()

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

    print(accepted[0])
    AR = accepted[0] / iterations * 100
    print(AR)

    print(np.mean(auxiliar))
    fig3 = figure()
    ax3 = fig3.add_subplot(111)
    ax3.plot(auxiliar)
    show()