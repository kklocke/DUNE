import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

import numpy as np

import seaborn as sns
sns.set()


# f = open("cellTest_cells.txt")
f = open("toy_cells.txt")

fig, ax = plt.subplots(2, 2)

patches = []
z1 = []
z2 = []
z3 = []
z4 = []

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

        # z1 = float(line[0])
        # z2 = float(line[1])
        # z3 = float(line[2])
        # z4 = float(line[3])
        # ax[0, 0].add_patch(Polygon(zip(x,y), facecolor='r', alpha=(z1 / 100.)))
        # ax[0, 1].add_patch(Polygon(zip(x,y), facecolor='r', alpha=(z2 / 100.)))
        # ax[1, 0].add_patch(Polygon(zip(x,y), facecolor='r', alpha=(z3 / 100.)))
        # ax[1, 1].add_patch(Polygon(zip(x,y), facecolor='r', alpha=(z4 / 100.)))
        patches.append(zip(x,y))
        z1.append(float(line[0]))
        z2.append(float(line[1]))
        z3.append(float(line[2]))
        z4.append(float(line[3]))
        x = []
        y = []
    lineCount += 1

f.close()

min1 = min(z1)
min2 = min(z2)
min3 = min(z3)
min4 = min(z4)
max1 = max(z1)
max2 = max(z2)
max3 = max(z3)
max4 = max(z4)

for i in range(len(z1)):
    a1 = 1
    if (max1 != min1):
        a1 = (z1[i] - min1) / (max1 - min1)
    if (min1 < 0):
        if (z1[i] < 0):
            a1 = z1[i] / min1
        else:
            a1 = z1[i] / max1
    a2 = 1
    if (max2 != min2):
        a2 = (z2[i] - min2) / (max2 - min2)
    if (min2 < 0):
        if (z2[i] < 0):
            a2 = z2[i] / min2
        else:
            a2 = z2[i] / max2
    a3 = 1
    if (max3 != min3):
        a3 = (z3[i] - min3) / (max3 - min3)
    if (min3 < 0):
        if (z3[i] < 0):
            a3 = z3[i] / min3
        else:
            a3 = z3[i] / max3
    a4 = 1
    if (max4 != min4):
        a4 = (z4[i] - min4) / (max4 - min4)
    if (min4 < 0):
        if (z4[i] < 0):
            a4 = z4[i] / min4
        else:
            a4 = z4[i] / max4

    if z1[i] > 0:
        ax[0, 0].add_patch(Polygon(patches[i], facecolor='r', alpha=a1))
    elif z1[i] == 0:
        ax[0, 0].add_patch(Polygon(patches[i], facecolor='white', alpha=a1))
    else:
        ax[0, 0].add_patch(Polygon(patches[i], facecolor='orange', alpha=a1))
    if z2[i] > 0:
        ax[0, 1].add_patch(Polygon(patches[i], facecolor='r', alpha=a2))
    elif z2[i] == 0:
        ax[0, 1].add_patch(Polygon(patches[i], facecolor='white', alpha=a2))
    else:
        ax[0, 1].add_patch(Polygon(patches[i], facecolor='orange', alpha=a2))
    if z3[i] > 0:
        ax[1, 0].add_patch(Polygon(patches[i], facecolor='r', alpha=a3))
    elif z3[i] == 0:
        ax[1, 0].add_patch(Polygon(patches[i], facecolor='white', alpha=a3))
    else:
        ax[1, 0].add_patch(Polygon(patches[i], facecolor='orange', alpha=a3))
    if z4[i] > 0:
        ax[1, 1].add_patch(Polygon(patches[i], facecolor='r', alpha=a4))
    elif z4[i] == 0:
        ax[1, 1].add_patch(Polygon(patches[i], facecolor='white', alpha=a4))
    else:
        ax[1, 1].add_patch(Polygon(patches[i], facecolor='orange', alpha=a4))




# f = open('cellTest_path.txt')
f = open('toy_path.txt')

lineCount = 0
x = []
y = []
for line in f:
    line = line.strip()
    if (len(line) == 0):
        ax[0, 0].plot(x, y, marker='D', color='green', lw=3, markerfacecolor='cyan', ms=10)
        ax[0, 1].plot(x, y, marker='D', color='green', lw=3, markerfacecolor='cyan', ms=10)
        ax[1, 0].plot(x, y, marker='D', color='green', lw=3, markerfacecolor='cyan', ms=10)
        ax[1, 1].plot(x, y, marker='D', color='green', lw=3, markerfacecolor='cyan', ms=10)
        x = []
        y = []
    else:
        line = line.split()
        x.append(float(line[0]))
        y.append(float(line[1]))

    # if (lineCount % 2 == 0):
    #     for elem in line.strip().split():
    #         x.append(float(elem))
    # else:
    #     for elem in line.strip().split():
    #         y.append(float(elem))
    #     #plt.plot(x, y, marker='D', color='green')
    #     ax[0, 0].plot(x, y, marker='D', color='green', lw=3)
    #     ax[0, 1].plot(x, y, marker='D', color='green', lw=3)
    #     ax[1, 0].plot(x, y, marker='D', color='green', lw=3)
    #     ax[1, 1].plot(x, y, marker='D', color='green', lw=3)
    #     x = []
    #     y = []
    # lineCount += 1
f.close()

ax[0, 0].plot(x, y, marker='D', color='green', lw=3, markerfacecolor='cyan', ms=10)
ax[0, 1].plot(x, y, marker='D', color='green', lw=3, markerfacecolor='cyan', ms=10)
ax[1, 0].plot(x, y, marker='D', color='green', lw=3, markerfacecolor='cyan', ms=10)
ax[1, 1].plot(x, y, marker='D', color='green', lw=3, markerfacecolor='cyan', ms=10)
x = []
y = []

# f = open('cellTest_wires.txt')
f = open('toy_wires.txt')

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
        ax[0, 0].plot(x, y, color='blue')
        ax[0, 1].plot(x, y, color='blue')
        ax[1, 0].plot(x, y, color='blue')
        ax[1, 1].plot(x, y, color='blue')
        # plt.plot(x, y, color='blue')
        x = []
        y = []
    lineCount += 1

f.close()



f = open('cellTest_scores.txt')

line = f.readline().strip().split()
#line2 = f.readline().strip().split()

f.close()

# t1 = "True Charge on Cells\nScore: " + line[0] + "   Cost: " + line2[0] + "      Color Scale: " + str(min1) + " to " +  str(max1)
# t2 = "Matrix Inversion\nScore: " + line[1] + "   Cost: " + line2[1] + "     Color Scale: " + str(min2) + " to " + str(max2)
# t3 = "Genetic Algorithm\nScore: " + line[2] + "   Cost: " + line2[2] + "     Color Scale: " + str(min3) + " to " +  str(max3)
# t4 = "Naive Randomization\nScore: " + line[3] + "   Cost: " + line2[3] + "     Color Scale: " + str(min4) + " to " +  str(max4)

t1 = "True Charge on Cells\nCost: " + line[0] + "      Color Scale: " + str(min1) + " to " +  str(max1)
t2 = "Matrix Inversion\nCost: " + line[1] + "     Color Scale: " + str(min2) + " to " + str(max2)
t3 = "Genetic Algorithm\nCost: " + line[2] + "     Color Scale: " + str(min3) + " to " +  str(max3)
t4 = "Naive Randomization\nCost: " + line[3] + "     Color Scale: " + str(min4) + " to " +  str(max4)

ax[0, 0].set_title(t1)
ax[0, 1].set_title(t2)
ax[1, 0].set_title(t3)
ax[1, 1].set_title(t4)



plt.show()
