import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

import numpy as np

lineCount = 0

startPts = []
endPts = []

f = open("displayTest.txt", 'r')

vertex = ''
pathEnd = ''

for line in f:
    line = line.strip()
    lineCount += 1
    if (lineCount == 1):
        vertex = line
        continue
    if (lineCount == 2):
        pathEnd = line
        continue
    if lineCount % 2 == 1:
        startPts.append(line)
    else:
        endPts.append(line)

f.close()

for i, elem1 in enumerate(startPts):
    elem2 = endPts[i]
    x = [float(elem1.split()[0]), float(elem2.split()[0])]
    y = [float(elem1.split()[1]), float(elem2.split()[1])]
    plt.plot(x,y, color='blue')

pathX = [float(vertex.split()[0]), float(pathEnd.split()[0])]
pathY = [float(vertex.split()[1]), float(pathEnd.split()[1])]

plt.plot(pathX, pathY, color='red')








f = open("hitLines.txt", "r")

lineCount = 0

startPts = []
endPts = []

for line in f:
    line = line.strip()
    lineCount += 1
    if lineCount % 2 == 1:
        startPts.append(line)
    else:
        endPts.append(line)

for i, elem1 in enumerate(startPts):
    elem2 = endPts[i]
    x = [float(elem1.split()[0]), float(elem2.split()[0])]
    y = [float(elem1.split()[1]), float(elem2.split()[1])]
    plt.plot(x,y, color='orange')


f.close()


f = open("displayPoints.txt", 'r')

x = []
y = []
z = []

for line in f:
    line = line.strip()
    x.append(float(line.split()[0]))
    y.append(float(line.split()[1]))
    z.append(float(line.split()[2]))
f.close()
for i, elem  in enumerate(z):
    z[i] = elem / max(z)

plt.scatter(x,y, marker='D', c=z, cmap='viridis')

plt.show()
plt.close("all")

fig, ax = plt.subplots()
patches = []

f = open("cellVertices.txt", 'r')

lineCount = 0

x = []
y = []
for line in f:
    line = line.strip()
    strList = line.split()
    for elem in strList:
        if (lineCount % 2 == 0):
            x.append(float(elem))
        else:
            y.append(float(elem))
    if (lineCount % 2 == 1):
        # myVerts = np.array([x, y])
        # myVerts.reshape((len(x), 2))
        # print(myVerts)
        # print(myVerts.shape)
        ax.add_patch(Polygon(zip(x,y), facecolor='r', alpha=0.5))
        x = []
        y = []
    lineCount += 1
    #if lineCount == 9:
    #    break

lineCount = 0

startPts = []
endPts = []

f = open("displayTest.txt", 'r')

vertex = ''
pathEnd = ''

for line in f:
    line = line.strip()
    lineCount += 1
    if (lineCount == 1):
        vertex = line
        continue
    if (lineCount == 2):
        pathEnd = line
        continue
    if lineCount % 2 == 1:
        startPts.append(line)
    else:
        endPts.append(line)

f.close()

for i, elem1 in enumerate(startPts):
    elem2 = endPts[i]
    x = [float(elem1.split()[0]), float(elem2.split()[0])]
    y = [float(elem1.split()[1]), float(elem2.split()[1])]
    plt.plot(x,y, color='blue')


plt.show()


plt.clf()

plt.close("all")
