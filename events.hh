#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <assert.h>
#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <armadillo>
#include "cell.hh"

using namespace std;
using namespace arma;


class Path
{
public:
	float length;
	float theta; // radial about the beam direction (0 to pi)
	float phi; // in plane with the beam / perp to the wire planes (0 to pi/2)
	Point vertex;
	Path(float r, float theta, float phi, Point v) : length(r), theta(theta), phi(phi), vertex(v) {}
	Path() : length(1.), theta(0.), phi(0.), vertex(Point(0., 0., 0.)) {}
	~Path() {};
	Point path2vec() {
		Point myVec;
		myVec.x = length * cos(phi);
		myVec.y = length * sin(phi) * sin(theta);
		myVec.z = length * sin(phi) * cos(theta);
		return myVec;
	}
	Path scale(float f) {
		return Path(length * f, theta, phi, vertex);
	}
};



Path randPath(float maxLeng) {
	// generate a random length between 0 and maxlength
	float myLength = maxLeng * (float(rand()) / float(RAND_MAX));
	float myTheta = M_PI * (float(rand()) / float(RAND_MAX));
	float myPhi = 0.5 * M_PI * (float(rand()) / float(RAND_MAX));
	Point orig = Point(0., 0., 0.);
	return Path(myLength, myTheta, myPhi, orig);
}










/* This class will describe naive events, wherein there is a
 * single well defined event vertex and two variable length
 * paths leading off from this and their angular orientation.
 */
class NaiveEvent
{
	public:
		Point vertex;
		Path path1;
		Path path2;
		NaiveEvent(Point v, Path p1, Path p2) : vertex(v), path1(p1), path2(p2) {}
		NaiveEvent() : vertex(Point()), path1(Path()), path2(Path()) {}
		~NaiveEvent() {};
		vector <vector <Path>> partitionEvent(int zPlanes, int zDim);
		// some attributes and functions here here
	private:
		// more attributes and functions here
};

NaiveEvent randEvent(int dims[]) {
	float x = float(0.1 * dims[0]) * float(rand()) / float(RAND_MAX);
	float y = float(0.1 * dims[1]) * float(rand()) / float(RAND_MAX);
	float z = float(0.1 * dims[2]) * float(rand()) / float(RAND_MAX);
	Point v = Point(x, y, z);
	float maxX = float(dims[0]) - x;
	float maxY = float(dims[1]) - y;
	float maxZ = float(dims[2]) - z;
	float maxLeng = min(min(maxX, maxY), maxZ); // Figure out which dim is appropriate later
	NaiveEvent myEvent = NaiveEvent(v, randPath(maxLeng), randPath(maxLeng));
	myEvent.path2.vertex = v;
	myEvent.path1.vertex = v;
	return myEvent;
}

vector <vector <Path>> NaiveEvent::partitionEvent(int zPlanes, int zDim) {
	vector <float> zPosns;
	float zIncr = float(zDim) / float(zPlanes - 1);
	// cout << "zIncr: " << zIncr << endl;
	float currPos = 0;
	while (currPos <= zDim)
	{
		// cout << "zPosn: " << currPos << endl;
		zPosns.push_back(currPos);
		currPos += zIncr;
	}
	// split path 1
	Point endPt1 = vertex + path1.path2vec();
	Point startPt1;
	int dirSign = 1;
	if (vertex.z < endPt1.z) {
		startPt1 = vertex;
	}
	else {
		startPt1 = endPt1;
		endPt1 = vertex;
		dirSign = -1;
	}
	int firstQuad = int(startPt1.z / zIncr);
	int lastQuad = int(endPt1.z / zIncr);
	vector <Path> p1_split;
	Point segStart = startPt1;
	float zLength = fabs(path1.path2vec().z);
	while (firstQuad <= lastQuad)
	{
		float endZ = min((firstQuad + 1)*zIncr, endPt1.z);
		// cout << "END Z: " << endZ << endl;
		// cout << "SegStart z: " << segStart.z << endl;
		Path myPath = path1.scale(dirSign *((endZ - segStart.z) / zLength));
		myPath.vertex = segStart;
		p1_split.push_back(myPath);
		segStart = segStart + myPath.path2vec();
		firstQuad++;
	}
	// split path 2
	// ignoring path 2 for now
	vector <vector <Path>> splitPaths;
	splitPaths.push_back(p1_split);
	return splitPaths;
}
