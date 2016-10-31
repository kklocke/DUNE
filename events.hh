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

using namespace std;
using namespace arma;

class Point
{
public:
	float x;
	float y;
	float z;
	Point(float x, float y, float z) : x(x), y(y), z(z) {}
	Point() : x(0.), y(0.), z(0.) {}
	~Point() {}
	void printPt() {cout << "(" << x << ", " << y << ", " << z << ")";};
	Point operator+(const Point& p) {
		Point sumPt;
		sumPt.x = this->x + p.x;
		sumPt.y = this->y + p.y;
		sumPt.z = this->z + p.z;
		return sumPt;
	}
	Point operator-(const Point& p) {
		Point diffPt;
		diffPt.x = this->x - p.x;
		diffPt.y = this->y - p.y;
		diffPt.z = this->z - p.z;
		return diffPt;
	}
	void operator=(const Point& p) {
		this->x = p.x;
		this->y = p.y;
		this->z = p.z;
	}
	bool operator<(const Point & p) const {
		if (x < p.x) {
			return true;
		}
		else if (x == p.x) {
			if (y < p.y) {
				return true;
			}
			else if (y == p.y) {
				if (z < p.z) {
					return true;
				}
			}
		}
		return false;
	}
	bool operator==(const Point& p) {
		return ((fabs(this->x - p.x) < 0.001) && (fabs(this->y - p.y) < 0.001) && (fabs(this->z - p.z) < 0.001));
	}
	bool operator!=(const Point& p) {
		return !(*this == p);
	}
	void round(int decPlace) {
		float delta = .5 * pow(10, -decPlace);
		float newX = ceil((x - delta)*pow(10, decPlace))/ pow(10, decPlace);
		float newY = ceil((y - delta)*pow(10, decPlace))/ pow(10, decPlace);
		float newZ = ceil((z - delta)*pow(10, decPlace))/ pow(10, decPlace);
		x = newX;
		y = newY;
		z = newZ;
	}
	Point scalarMult(float f) {
		return Point(x * f, y * f, z*f);
	}
};

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


class wireArray
{
public:
	float pitch; // mm
	float angle; // angle with respect to x axis
	float height; // y dim
	float length; // x dim
	float z; // z posn
    vector<Point> startPoints;
    vector<Point> endPoints;
	wireArray(float pitch, float angle, float l, float h, float z) : pitch(pitch), angle(angle), height(h), length(l), z(z) {
		float xPos = 0., yPos = h;
		float xIncr = pitch / fabs(sin(M_PI * angle / 180));
		float yIncr = pitch / fabs(cos(M_PI * angle / 180));
		double slope = tan(M_PI * angle / 180.);
		if (int(angle) == 90)
		{
			while (xPos < l)
			{
				startPoints.push_back(Point(xPos, h, z));
				endPoints.push_back(Point(xPos, 0, z));
				xPos += xIncr;
			}
		}
		else if (angle == 0.)
		{
			while (yPos > 0)
			{
				startPoints.push_back(Point(l, yPos, z));
				endPoints.push_back(Point(0, yPos, z));
				yPos -= yIncr;
			}
		}
		else if (slope > 0)
		{
			while (xPos < l)
			{
				startPoints.push_back(Point(xPos, yPos, z));
				float yInt = yPos - xPos * slope;
				if ((yInt >= 0) && (yInt <= h))
				{
					endPoints.push_back(Point(0., yInt, z));
				}
				else
				{
					float xInt = xPos - (yPos / slope);
					endPoints.push_back(Point(xInt, 0., z));
				}
				xPos += xIncr;
			}
			xPos = l;
			yPos = h;
			while (yPos > 0)
			{
				startPoints.push_back(Point(xPos, yPos, z));
				float xInt = l - yPos / slope;
				if ((xInt > 0) && (xInt < l))
				{
					endPoints.push_back(Point(xInt, 0., z));
				}
				else
				{
					float yInt = yPos - l*slope;
					endPoints.push_back(Point(0., yInt, z));
				}
				yPos -= yIncr;
			}
		}
		else if (slope < 0)
		{
			while (xPos < h)
			{
				startPoints.push_back(Point(xPos, yPos, z));
				float xInt = xPos - h / slope;
				if ((xInt > 0) && (xInt < h))
				{
					endPoints.push_back(Point(xInt, 0., z));
				}
				else
				{
					float yInt = yPos  + slope * (h - xPos);
					endPoints.push_back(Point(l, yInt, z));
				}
				xPos += xIncr;
			}
			xPos = 0.;
			yPos = 0.;
			while (yPos < l)
			{
				startPoints.push_back(Point(0, yPos, z));
				float xInt = xPos - yPos / slope;
				if ((xInt > 0) && (xInt < l))
				{
					endPoints.push_back(Point(xInt, 0., z));
				}
				else
				{
					float yInt = yPos - slope * h;
					endPoints.push_back(Point(l, yInt, z));
				}
				yPos += yIncr;
			}
		}
	}
	wireArray() : pitch(0.), angle(0.), height(0.), length(0.), z(0.) {};
};



vector <int> crossings(wireArray myArr, Point p1, Point p2)
{
    // figure out which lines are crossed in the wire array
	if (p1.x > p2.x)
	{
		Point p3 = p1;
		p1 = p2;
		p2 = p3;
	}
    vector <int> myCross;
	int firstWire = -1;
	int lastWire = -1;
	// cout << "my points to cross: \n\t";
	// p1.printPt();
	// cout << "\n\t";
	// p2.printPt();
	// cout << endl;
	for (int i = 0; i < (int)myArr.endPoints.size(); i++)
	{
		Point refPt1 = myArr.startPoints[i];
		Point refPt2 = myArr.endPoints[i];
		float slope = (refPt2.y - refPt1.y) / (refPt2.x - refPt1.x);
		// y at given x
		float refY_1 = refPt2.y + (slope * (p1.x - refPt2.x));
		float refY_2 = refPt2.y + (slope * (p2.x - refPt2.x));
		if ((slope > 0) && (p1.y >= refY_1) && (firstWire == -1))
		{
			firstWire = i;
		}
		else if ((slope < 0) && (p1.y <= refY_1) && (firstWire == -1))
		{
			firstWire = i;
		}
		if ((slope > 0) && (p2.y > refY_2) && (lastWire == -1))
		{
			lastWire = i;
		}
		else if ((slope < 0) && (p2.y < refY_2) && (lastWire == -1))
		{
			lastWire = i;
		}
		if ((firstWire != -1) && (lastWire != -1))
		{
			break;
		}
	}
	if (lastWire == -1) {
		lastWire = (int)myArr.endPoints.size();
	}
	// cout << "FIRST WIRE: " << firstWire << "\tLAST WIRE: " << lastWire << endl;
	for (int i = min(firstWire, lastWire); i < max(lastWire, firstWire); i++)
	{
		myCross.push_back(i);
		// cout << i << "\t";
	}
	// cout << endl;
    return myCross;
}


Point intersection(Point p1_s, Point p1_e, Point p2_s, Point p2_e) {
	Point vec1 = p1_e - p1_s;
	Point vec2 = p2_e - p2_s;
	Point diff_s = p1_s - p2_s;
	if (fabs(vec1.x) < 0.0001) {
		Point newPt = Point(p1_s.x, 0, p1_s.z);
		float scale = fabs(diff_s.x / vec2.x);
		newPt.y = p2_s.y + scale * vec2.y;
		return newPt;
	}
	if (fabs(vec2.x) < 0.0001) {
		Point newPt = Point(p2_s.x, 0, p2_s.z);
		float scale = fabs(diff_s.x / vec1.x);
		newPt.y = p1_s.y + scale * vec1.y;
		return newPt;
	}
	if (vec1.y == 0) {
		// there may be a sign error here
		float yDiff = diff_s.y;
		float yFrac = yDiff / vec2.y;
		return p2_s + vec2.scalarMult(yFrac);
	}
	if (vec2.y == 0) {
		// there may be a sign error here
		float yDiff = diff_s.y;
		float yFrac = yDiff / vec1.y;
		return p1_s + vec1.scalarMult(yFrac);
	}
	float t = ((diff_s.y/vec2.y) - (diff_s.x/vec2.x))/((vec1.x / vec2.x) - (vec1.y/vec2.y));
	float s = (diff_s.x + t * vec1.x) / vec2.x;
	if (p1_s + vec1.scalarMult(t) == p2_s + vec2.scalarMult(s)) {
		return p1_s + vec1.scalarMult(t);
	}
	return Point(-1., -1., -1.);
}

class wireLayer
{
public:
    wireArray vertical; // theta = Pi/2
    wireArray lSlant; // theta = Pi/4
    wireArray rSlant; // theta = -Pi/4
	float pitch;
	float length;
	float height;
	float z;
	vector <Point> grids;
	wireLayer(float h, float l, float z, float pitch) : pitch(pitch), length(l), height(h), z(z) {
		vertical = wireArray(pitch, 90., l, h, z);
		lSlant = wireArray(pitch, 45., l, h, z);
		rSlant = wireArray(pitch, -45., l, h, z);
	}
	vector <vector <int>> allCrossings(Point p1, Point p2) {
		vector <vector <int>> all3;
		all3.push_back(crossings(vertical, p1, p2));
		all3.push_back(crossings(lSlant, p1, p2));
		all3.push_back(crossings(rSlant, p1, p2));
		return all3;
	}
	vector <vector <int>> geoMatrix()
	{
		int vSize = vertical.startPoints.size();
		int lSize = lSlant.startPoints.size();
		int rSize = rSlant.startPoints.size();
		int numRows = vSize + lSize + rSize;
		vector <vector <int>> myMatrix;
		map<Point, vector<int>> gridPts;

		for (int i = 0; i < (int)lSlant.startPoints.size(); i++)
		{
			Point p1 = lSlant.startPoints[i];
			Point p2 = lSlant.endPoints[i];
			vector <int> wires = crossings(vertical, p1, p2);
			for (int j = 0; j < (int)wires.size(); j++) {
				Point refp1 = vertical.startPoints[wires[j]];
				Point refp2 = vertical.endPoints[wires[j]];
				Point intPoint = intersection(p1, p2, refp1, refp2);
				if (intPoint != Point(-1., -1., -1.)) {
					intPoint.round(2);
					// cout << "Point found: ";
					// intPoint.printPt();
					// cout << endl;
					if (gridPts.find(intPoint) == gridPts.end())
					{
						int wireRefs [3] = {wires[j], i, -1};
						vector<int> newEntry;
						newEntry.assign(wireRefs, wireRefs+3);
						gridPts[intPoint] = newEntry;
					}
					else {
						if (gridPts[intPoint][0] == -1) {
							gridPts[intPoint][0] = wires[j];
						}
						if (gridPts[intPoint][1] == -1) {
							gridPts[intPoint][1] = i;
						}
					}
				}
			}
		}
		for (int i = 0; i < (int)rSlant.startPoints.size(); i++)
		{
			Point p1 = rSlant.startPoints[i];
			Point p2 = rSlant.endPoints[i];
			vector <int> wires = crossings(vertical, p1, p2);
			for (int j = 0; j < (int)wires.size(); j++) {
				Point refp1 = vertical.startPoints[wires[j]];
				Point refp2 = vertical.endPoints[wires[j]];
				Point intPoint = intersection(p1, p2, refp1, refp2);
				if (intPoint != Point(-1., -1., -1.)) {
					intPoint.round(2);
					// cout << "Point found: ";
					// intPoint.printPt();
					// cout << endl;
					if (gridPts.find(intPoint) == gridPts.end())
					{
						int wireRefs [3] = {wires[j], -1, i};
						vector<int> newEntry;
						newEntry.assign(wireRefs, wireRefs+3);
						gridPts[intPoint] = newEntry;
					}
					else {
						if (gridPts[intPoint][0] == -1) {
							gridPts[intPoint][0] = wires[j];
						}
						if (gridPts[intPoint][2] == -1) {
							gridPts[intPoint][2] = i;
						}
					}
				}
			}
		}
		for (int i = 0; i < (int)rSlant.startPoints.size(); i++)
		{
			Point p1 = rSlant.startPoints[i];
			Point p2 = rSlant.endPoints[i];
			vector <int> wires = crossings(lSlant, p1, p2);
			for (int j = 0; j < (int)wires.size(); j++) {
				Point refp1 = lSlant.startPoints[wires[j]];
				Point refp2 = lSlant.endPoints[wires[j]];
				Point intPoint = intersection(p1, p2, refp1, refp2);
				if (intPoint != Point(-1., -1., -1.)) {
					intPoint.round(2);
					// cout << "Point found: ";
					// intPoint.printPt();
					// cout << endl;
					if ((intPoint.x <= 0) || (intPoint.y <= 0)) {
						continue;
					}
					if (gridPts.find(intPoint) == gridPts.end())
					{
						int wireRefs [3] = {-1, wires[j], i};
						vector<int> newEntry;
						newEntry.assign(wireRefs, wireRefs+3);
						gridPts[intPoint] = newEntry;
					}
					else {
						if (gridPts[intPoint][1] == -1) {
							gridPts[intPoint][1] = wires[j];
						}
						if (gridPts[intPoint][2] == -1) {
							gridPts[intPoint][2] = i;
						}
					}
				}
			}
		}
		int numCols = gridPts.size();
		for (int i = 0; i < numRows; i++) {
			vector <int> temp(numCols, 0);
			myMatrix.push_back(temp);
		}
		int count = 0;
		this->grids.clear();
		for (map<Point, vector<int>>::iterator j = gridPts.begin(); j != gridPts.end(); j++) {
			//j->first = key
			//j->second = vertical
			Point myPt = j->first;
			// cout << "Intersection Point: ";
			// myPt.printPt();
			// cout <<"\n";
			this->grids.push_back(myPt);
			vector<int> my3 = j->second;
			if (my3[0] != -1) {
				myMatrix[my3[0]][count] = 1;
			}
			if (my3[1] != -1) {
				myMatrix[my3[1] + vSize][count] = 1;
			}
			if (my3[0] != -1) {
				myMatrix[my3[2] + vSize + lSize][count] = 1;
			}
			count++;
		}
		return myMatrix;
	}
	vector <float> signalVec(Path p) {
		int vDim = this->vertical.startPoints.size();
		int lDim = this->lSlant.startPoints.size();
		int rDim = this->rSlant.startPoints.size();
		vector <int> crossV = crossings(this->vertical, p.vertex, (p.vertex + p.path2vec()));
		vector <int> crossL = crossings(this->lSlant, p.vertex, (p.vertex + p.path2vec()));
		vector <int> crossR = crossings(this->rSlant, p.vertex, (p.vertex + p.path2vec()));
		vector <float> sigVec(vDim + lDim + rDim, 0.);
		for (int i = 0; i < (int)crossV.size(); i++) {
			sigVec[crossV[i]] = 1.;
		}
		for (int i = 0; i < (int)crossL.size(); i++) {
			sigVec[crossL[i] + vDim] = 1.;
		}
		for (int i = 0; i < (int)crossR.size(); i++) {
			sigVec[crossR[i] + vDim + lDim] = 1.;
		}
		return sigVec;
	}
};


vector <float> solveTrue(vector <float> mySig, vector <vector <int>> geoMat) {
	mat geo(geoMat.size(), geoMat[0].size());
	vec sig(mySig.size());
	for (int i = 0; i < (int)geoMat.size(); i++) {
		for (int j = 0; j < (int)geoMat[i].size(); j++) {
			geo(i, j) = geoMat[i][j];
		}
	}
	for (int i = 0; i < (int)mySig.size(); i++) {
		sig[i] = mySig[i];
	}
	vec trueSig = solve(geo, sig);
	vector <float> trueSignal;
	for (int i = 0; i < (int)geoMat[0].size(); i++) {
		trueSignal.push_back(trueSig[i]);
	}
	return trueSignal;
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
	float x = float(dims[0]) * float(rand()) / float(RAND_MAX);
	float y = float(dims[1]) * float(rand()) / float(RAND_MAX);
	float z = float(dims[2]) * float(rand()) / float(RAND_MAX);
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
