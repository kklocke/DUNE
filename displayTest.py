import matplotlib.pyplot as plt

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

for line in f:
    line = line.strip()
    x = float(line.split()[0])
    y = float(line.split()[1])
    plt.plot(x,y, marker='o', color='green')
f.close()


plt.show()
