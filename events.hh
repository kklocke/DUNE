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
#include "blob.hh"

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
		// Point myVec;
		float myX = length * cos(phi);
		float myY = length * sin(phi) * sin(theta);
		float myZ = length * sin(phi) * cos(theta);
		return Point(myX, myY, myZ);
	}
	Path scale(float f) {
		return Path(length * f, theta, phi, vertex);
	}
	void operator=(const Path &p) {
		this->length = p.length;
		this->theta = p.theta;
		this->phi = p.phi;
		this->vertex = p.vertex;
	}
	bool operator!=(const Path &p) {
		if (fabs(this->length - p.length) > 0.0001) {
			return true;
		}
		if (fabs(this->theta - p.theta) > 0.0001) {
			return true;
		}
		if (fabs(this->phi - p.phi) > 0.0001) {
			return true;
		}
		if (this->vertex != p.vertex) {
			return true;
		}
		return false;
	}
	bool operator==(const Path &p) {
		return !(*this != p);
	}
};



Path randPath(float maxLeng) {
	// generate a random length between 0 and maxlength
	float myLength = maxLeng * (float(rand()) / float(RAND_MAX));
	float myTheta = 0.5 * M_PI * (float(rand()) / float(RAND_MAX));
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
	vector <vector <Path>> allBins(2);
	float zStart1 = max(float(0.), path1.vertex.z);
	float zStart2 = max(float(0.), path2.vertex.z);
	float binSize = float(zDim) / float(zPlanes);
	float zLength1 = path1.path2vec().z;
	float zLength2 = path2.path2vec().z;

	vector <Path> p1Bin(zPlanes);
	vector <Path> p2Bin(zPlanes);

	for (int i = 0; i < zPlanes; i++) {
		// vector <Path> thisBin(0);

		float binStart = float(i) * binSize;
		float binEnd = float(i + 1) * binSize;

		cout << "Bin: " << binStart << " --> " << binEnd << endl;
		cout << "zStart1: " << zStart1 << "\tzStart2: " << zStart2 << endl;

		if ((zStart1 >= binStart) && (zStart1 <= binEnd)) {
			cout << "true 1\n";
			float zEnd1 = min(path1.vertex.z + zLength1, binEnd);
			Point startPt1 = path1.vertex + path1.path2vec() * ((zStart1 - path1.vertex.z) / zLength1);
			// Point endPt1 = startPt1 + path1.path2vec() * ((zEnd1 - zStart1) / zLength1);
			Path newP = path1;
			newP.vertex = startPt1;
			newP.length *= (zEnd1 - zStart1) / zLength1;
			p1Bin[i] = newP;
			// thisBin.push_back(newP);
			zStart1 = zEnd1;
		}

		if ((zStart2 >= binStart) && (zStart2 <= binEnd)) {
			cout << "true 2\n";
			float zEnd2 = min(path2.vertex.z + zLength2, binEnd);
			Point startPt2 = path2.vertex + path2.path2vec() * ((zStart2 - path2.vertex.z) / zLength2);
			// Point endPt1 = startPt1 + path1.path2vec() * ((zEnd1 - zStart1) / zLength1);
			Path newP = path2;
			newP.vertex = startPt2;
			newP.length *= (zEnd2 - zStart2) / zLength2;
			p2Bin[i] = newP;
			// thisBin.push_back(newP);
			zStart2  = zEnd2;
		}
		// allBins[i] = thisBin;
	}
	allBins[0] = p1Bin;
	allBins[1] = p2Bin;
	return allBins;
}










//
// vector <vector <Path>> NaiveEvent::partitionEvent(int zPlanes, int zDim) {
// 	vector <float> zPosns;
// 	float zIncr = float(zDim) / float(zPlanes - 1);
// 	// cout << "zIncr: " << zIncr << endl;
// 	float currPos = 0;
// 	while (currPos <= zDim)
// 	{
// 		// cout << "zPosn: " << currPos << endl;
// 		zPosns.push_back(currPos);
// 		currPos += zIncr;
// 	}
// 	// split path 1
// 	Point endPt1 = vertex + path1.path2vec();
// 	Point startPt1;
// 	int dirSign = 1;
// 	if (vertex.z < endPt1.z) {
// 		startPt1 = vertex;
// 	}
// 	else {
// 		startPt1 = endPt1;
// 		endPt1 = vertex;
// 		dirSign = -1;
// 	}
// 	int firstQuad = int(startPt1.z / zIncr);
// 	int lastQuad = int(endPt1.z / zIncr);
// 	vector <Path> p1_split;
// 	Point segStart = startPt1;
// 	float zLength = fabs(path1.path2vec().z);
// 	while (firstQuad <= lastQuad)
// 	{
// 		float endZ = min((firstQuad + 1)*zIncr, endPt1.z);
// 		// cout << "END Z: " << endZ << endl;
// 		// cout << "SegStart z: " << segStart.z << endl;
// 		Path myPath = path1.scale(dirSign *((endZ - segStart.z) / zLength));
// 		myPath.vertex = segStart;
// 		p1_split.push_back(myPath);
// 		segStart = segStart + myPath.path2vec();
// 		firstQuad++;
// 	}
// 	// split path 2
// 	// ignoring path 2 for now
// 	vector <vector <Path>> splitPaths;
// 	splitPaths.push_back(p1_split);
// 	return splitPaths;
// }
