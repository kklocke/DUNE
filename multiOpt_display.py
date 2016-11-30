import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

import numpy as np

f = open("multiOpt_layer.txt")

for line in f:
    line = line.strip()
    points = line.split()
    plt.plot([float(points[0]), float(points[2])],[float(points[1]), float(points[3])], color='blue')
f.close()

x = []
y = []

f = open("multiOpt_path.txt")

for line in f:
    line = line.strip()
    line = line.split()
    x.append(float(line[0]))
    x.append(float(line[2]))
    y.append(float(line[1]))
    y.append(float(line[3]))
f.close()

plt.plot(x, y, color='red', marker='D')

f = open("multiOpt_charges.txt")

x = []
y = []
val = []
avgVal = 0
lineCount = 0

for line in f:
    line = line.strip()
    line = line.split()
    x.append(float(line[0]))
    y.append(float(line[1]))
    val.append(float(line[2]))
    lineCount += 1
    avgVal += float(line[2])

f.close()

avgVal /= float(lineCount)

for i, v in enumerate(val):
    val[i] = v/avgVal

plt.scatter(x, y, cmap='viridis', c=val, marker='D')

plt.show()
