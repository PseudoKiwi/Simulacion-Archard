from numpy import e, abs
from random import random

# The defineProbability function receives the parameter beta, such that 1/beta=kB*T. It then returns
# the probability function to be used during the simulation asuming negligible changes in temperature.


def defineProbability(dENeg, accepted, beta):
    def acceptance(dE):
        if( dE < 0 ):   # ¿ < o <= ?
            dENeg[0] += 1
            return True
        else:
            r = random()
            print([e**(-beta*dE), r, dE])
            aux = r < e**(-beta*dE)
            if aux:
                accepted[0] += 1
            return aux
    return acceptance