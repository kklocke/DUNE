#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <assert.h>
#include <iostream>

using namespace std;


class Point
{
public:
	float x;
	float y;
	float z;
	Point(float x, float y, float z) : x(x), y(y), z(z) {}
	Point() : x(0.), y(0.), z(0.) {}
	~Point() {};
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
};

Path randPath(float maxLeng) {
	// generate a random length between 0 and maxlength
	srand(time(NULL));
	float myLength = maxLeng * (float(rand()) / float(RAND_MAX));
	float myTheta = M_PI * (float(rand()) / float(RAND_MAX));
	float myPhi = 0.5 * M_PI * (float(rand()) / float(RAND_MAX));
	return Path(myLength, myTheta, myPhi);
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
		// some attributes and functions here here
	private:
		// more attributes and functions here
};

NaiveEvent randEvent(int dims[]) {
	// generate random vertex
	srand(time(NULL));
	float x = float(dims[0]) * float(rand()) / float(RAND_MAX);
	float y = float(dims[1]) * float(rand()) / float(RAND_MAX);
	float z = float(dims[2]) * float(rand()) / float(RAND_MAX);
	Point v = Point(x, y, z);
	float maxLeng = float(dims[0])  - x; // Figure out which dim is appropriate later
	return NaiveEvent(v, randPath(maxLeng), randPath(maxLeng));
}