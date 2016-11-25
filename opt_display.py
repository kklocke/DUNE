import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

import numpy as np

lineCount = 0

startPts = []
endPts = []

f = open("opt_vertex_wires.txt", 'r')

vertex = ''
pathEnd = ''

for line in f:
    line = line.strip()
    lineCount += 1
    if (lineCount <= 2):
        if (lineCount % 2 == 1):
            vertex = line
        else:
            pathEnd = line
        continue
    if lineCount % 2 == 1:
        startPts.append(line)
    else:
        endPts.append(line)
f.close()

gridX = []
gridY = []
for i, elem1 in enumerate(startPts):
    elem2 = endPts[i]
    x = [float(elem1.split()[0]), float(elem2.split()[0])]
    gridX.append(x)
    y = [float(elem1.split()[1]), float(elem2.split()[1])]
    gridY.append(y)
    plt.plot(x,y, color='blue')


plt.plot([vertex.split()[0], pathEnd.split()[0]], [vertex.split()[1], pathEnd.split()[1]], color='r')

f = open("opt_hitlines.txt", "r")

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


f = open("opt_displayPoints.txt", 'r')

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
