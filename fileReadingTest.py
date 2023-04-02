from numpy import zeros, sqrt, mean
from sys import maxsize


nParticles = 100
position = zeros([2, nParticles])            # Material particles positions


with open("materialEq.txt", "r") as file:
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


req = zeros([nParticles])
for i in range(nParticles):
    minimum = maxsize
    for j in range(nParticles):
        if  i != j:
            r = sqrt((position[0][i] - position[0][j]) ** 2 + (position[1][i] - position[1][j]) ** 2)
            if r < minimum:
                minimum = r
    req[i] = minimum
req = mean(req)
print(req)

with open("r01.txt", "w") as file:
    file.write(str(req))