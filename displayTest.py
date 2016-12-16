import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

import numpy as np

lineCount = 0

startPts = []
endPts = []

f = open("displayTest.txt", 'r')

vertex = []
pathEnd = []
numPaths = 1

for line in f:
    line = line.strip()
    lineCount += 1
    if (lineCount == 1):
        numPaths = int(line)
        continue
    if (lineCount <= 2*numPaths + 1):
        if (lineCount % 2 == 0):
            vertex.append(line)
        else:
            pathEnd.append(line)
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
    #plt.plot(x,y, color='blue')

pathX = []
pathY = []
for i, elem1 in enumerate(vertex):
    elem2 = pathEnd[i]
    x = [float(elem1.split()[0]), float(elem2.split()[0])]
    y = [float(elem1.split()[1]), float(elem2.split()[1])]
    pathX.append(x)
    pathY.append(y)


#plt.plot(pathX, pathY, color='red')








f = open("hitLines.txt", "r")

# lineCount = 0
#
# startPts = []
# endPts = []
#
# for line in f:
#     line = line.strip()
#     lineCount += 1
#     if lineCount % 2 == 1:
#         startPts.append(line)
#     else:
#         endPts.append(line)
#
# for i, elem1 in enumerate(startPts):
#     elem2 = endPts[i]
#     x = [float(elem1.split()[0]), float(elem2.split()[0])]
#     y = [float(elem1.split()[1]), float(elem2.split()[1])]
#     plt.plot(x,y, color='orange')
#
#
# f.close()


lineCount = 0
startPts = []
endPts = []
secCount = 0
startPtsAll = []
endPtsAll = []

for line in f:
    line = line.strip()
    if line == "NEXT PARTITION":
        secCount += 1
        startPtsAll.append(startPts)
        endPtsAll.append(endPts)
        startPts = []
        endPts = []
        lineCount = 0
        continue
    lineCount += 1
    if lineCount % 2 == 1:
        startPts.append(line)
    else:
        endPts.append(line)
f.close()


f = open("displayPoints.txt", 'r')

xAll = []
yAll = []
zAll = []
x = []
y = []
z = []
for line in f:
    line = line.strip()
    if (line == "NEXT PARTITION"):
        xAll.append(x)
        yAll.append(y)
        zAll.append(z)
        x = []
        y = []
        z = []
        continue
    x.append(float(line.split()[0]))
    y.append(float(line.split()[1]))
    z.append(float(line.split()[2]))
f.close()


# now display
for i in range(numPaths):
    # plot the wire array
    # plot the point
    # plot the path
    plt.plot(pathX[i], pathY[i], color='red')
    for j in range(len(gridX)):
        plt.plot(gridX[j], gridY[j], color='blue')
    for j in range(len(startPtsAll[i])):
        elem1 = startPtsAll[i][j]
        elem2 = endPtsAll[i][j]
        x = [float(elem1.split()[0]), float(elem2.split()[0])]
        y = [float(elem1.split()[1]), float(elem2.split()[1])]
        plt.plot(x, y, color='orange')
    myElems = zip(xAll[i], yAll[i])
    for j, elem in enumerate(myElems):
        if elem == max(xAll[i]):
            plt.scatter(elem, yAll[i][j], marker='D', color='green')
    # for j in range(len(xAll[i])):
    #     elem1 = xAll[i][j]
    #     elem2 = yAll[i][j]
    #     plt.scatter(elem1, elem2, marker='D', color='green')
    plt.show()


#
# f = open("displayPoints.txt", 'r')
#
# x = []
# y = []
# z = []
#
# for line in f:
#     line = line.strip()
#     x.append(float(line.split()[0]))
#     y.append(float(line.split()[1]))
#     z.append(float(line.split()[2]))
# f.close()
# for i, elem  in enumerate(z):
#     z[i] = elem / max(z)
# d
# plt.scatter(x,y, marker='D', c=z, cmap='viridis')
#
# plt.show()
# plt.close("all")

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
f.close()
for j in range(len(gridX)):
    plt.plot(gridX[j], gridY[j], color='blue')

for i in range(numPaths):
    plt.plot(pathX[i], pathY[i], color='green', marker='D')

plt.show()
# lineCount = 0
#
# startPts = []
# endPts = []
#
# f = open("displayTest.txt", 'r')
#
# vertex = ''
# pathEnd = ''
#
# for line in f:
#     line = line.strip()
#     lineCount += 1
#     if (lineCount == 1):
#         vertex = line
#         continue
#     if (lineCount == 2):
#         pathEnd = line
#         continue
#     if lineCount % 2 == 1:
#         startPts.append(line)
#     else:
#         endPts.append(line)
#
# f.close()
#
# for i, elem1 in enumerate(startPts):
#     elem2 = endPts[i]
#     x = [float(elem1.split()[0]), float(elem2.split()[0])]
#     y = [float(elem1.split()[1]), float(elem2.split()[1])]
#     plt.plot(x,y, color='blue')
#
#
# plt.show()
#
#
# plt.clf()
#
# plt.close("all")
