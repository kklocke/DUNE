#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <assert.h>
#include <iostream>
#include <vector>

using namespace std;


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
	bool operator==(const Point& p) {
		return ((this->x == p.x) && (this->y == p.y) && (this->z == p.z));
	}
};

class Path
{
public:
	float length;
	float theta; // radial about the beam direction (0 to pi)
	float phi; // in plane with the beam / perp to the wire planes (0 to pi/2)
	Path(float r, float theta, float phi) : length(r), theta(theta), phi(phi) {}
	Path() : length(1.), theta(0.), phi(0.) {}
	~Path() {};
	Point path2vec() {
		Point myVec;
		myVec.x = length * cos(phi);
		myVec.y = length * sin(phi) * sin(theta);
		myVec.z = length * sin(phi) * cos(theta);
		return myVec;
	}
};

Path randPath(float maxLeng) {
	// generate a random length between 0 and maxlength
	float myLength = maxLeng * (float(rand()) / float(RAND_MAX));
	float myTheta = M_PI * (float(rand()) / float(RAND_MAX));
	float myPhi = 0.5 * M_PI * (float(rand()) / float(RAND_MAX));
	return Path(myLength, myTheta, myPhi);
}


class wireArray
{
public:
	int pitch; // mm
	float angle; // angle with respect to x axis
    vector<Point> startPoints;
    vector<Point> endPoints;
	wireArray(float pitch, float angle, float l, float h, float z) : pitch(pitch), angle(angle)
	{
		float xPos = 0., yPos = h;
		float xIncr = pitch / sin(M_PI * angle / 180);
		float yIncr = pitch / cos(M_PI * angle / 180);
		double slope = tan(M_PI * angle / 180.);
		if (slope > 0)
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
	};
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
		if ((slope < 0) && (p1.y <= refY_1) && (firstWire == -1))
		{
			firstWire = i;
		}
		if ((slope > 0) && (p2.y > refY_2) && (lastWire == -1))
		{
			lastWire = i;
		}
		if ((slope < 0) && (p2.y < refY_2) && (lastWire == -1))
		{
			lastWire = i;
		}
		if ((firstWire != -1) && (lastWire != -1))
		{
			break;
		}
	}
	cout << "First Wire " << firstWire << endl;
	cout << "Last Wire " << lastWire << endl;

	for (int i = min(firstWire, lastWire); i < max(lastWire, firstWire); i++) {
		myCross.push_back(i);
	}
    return myCross;
}

class wireLayer
{
    wireArray vertical;
    wireArray lSlant;
    wireArray rSlant;
};

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
		// some attributes and functions here here
	private:
		// more attributes and functions here
};

NaiveEvent randEvent(int dims[]) {
	float x = float(dims[0]) * float(rand()) / float(RAND_MAX);
	float y = float(dims[1]) * float(rand()) / float(RAND_MAX);
	float z = float(dims[2]) * float(rand()) / float(RAND_MAX);
	Point v = Point(x, y, z);
	float maxLeng = float(dims[0])  - x; // Figure out which dim is appropriate later
	return NaiveEvent(v, randPath(maxLeng), randPath(maxLeng));
}
