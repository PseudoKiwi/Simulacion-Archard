
# The potential function receives the parameters of the Lennard-Jones
# potential energy and returns the function to aply with only the
# relative distances between particles.

# U0: Measures how strong particles atract each other
# r0: Equilibrium radius

# 2*r0 is defined here as the critical radius where the interaction becomes
# negligible
def potential(U0, r0):
    def energy(r):
        return U0*((r0/r)**12 - (r0/r)**6)
    return energy