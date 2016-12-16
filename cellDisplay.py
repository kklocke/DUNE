import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

import numpy as np

f = open("cellTest_cells.txt")

fig, ax = plt.subplots()
patches = []

lineCount = 0
x = []
y = []
for line in f:
    line = line.strip().split()
    z = 0
    elemCount = 0
    if (lineCount % 2 == 0):
        for elem in line:
            if (elemCount % 2 == 0):
                x.append(float(elem))
            else:
                y.append(float(elem))
            elemCount += 1
    else:
        z1 = float(line[0])
        z2 = float(line[1])
        ax.add_patch(Polygon(zip(x,y), facecolor='r', alpha=(z1 / 10.)))
        ax.add_patch(Polygon(zip(x,y), facecolor='b', alpha=(z2 / 20.)))
        x = []
        y = []
    lineCount += 1

f.close()

f = open('cellTest_path.txt')

lineCount = 0
x = []
y = []
for line in f:
    if (lineCount == 0):
        for elem in line.strip().split():
            x.append(float(elem))
    else:
        for elem in line.strip().split():
            y.append(float(elem))
    lineCount += 1
f.close()
plt.plot(x, y)

f = open('cellTest_wires.txt')

lineCount = 0
x = []
y = []
for line in f:
    line = line.strip().split()
    if (lineCount % 2 == 0):
        for elem in line:
            x.append(float(elem))
    else:
        for elem in line:
            y.append(float(elem))
        plt.plot(x, y, color='blue')
        x = []
        y = []
    lineCount += 1

f.close()
plt.show()







plt.show()
