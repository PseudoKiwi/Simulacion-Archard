from numpy import array

# The particle class permits to distinguish one another by the type atribute.
# Different materials will be defined by the different types.
# The position is needed to compute the interaction between particles.

class MaterialParticle:

    def __init__(self, pos, particleType):
        self.pos = pos              # numpy array
        self.type = particleType    # Integer value
        self.increment = array([0,0], float)

    def set_pos(self, newPos):
        self.pos = newPos

    def getPos(self):
        return self.pos

    def getType(self):
        return self.type

    def modifyPos(self):
        self.pos += self.increment

    def setIncrement(self, newIncrement):
        self.increment = newIncrement

    def getIncrement(self):
        return self.increment