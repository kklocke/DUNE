#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <assert.h>
#include <iostream>
#include <cstdio>

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
	void operator+=(const Point& p) {
		this->x += p.x;
		this->y += p.y;
		this->z += p.z;
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
	bool pointMatch(const Point& p) {
		return ((fabs(this->x - p.x) < 0.01) && (fabs(this->y - p.y) < 0.01));
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

Point intersection(Point p1_s, Point p1_e, Point p2_s, Point p2_e) {
	Point vec1 = p1_e - p1_s;
	Point vec2 = p2_e - p2_s;
	Point diff_s = p1_s - p2_s;
	if (fabs(vec1.x) < 0.0001) {
		Point newPt = Point(p1_s.x, 0, p1_s.z);
		float slope = vec2.y / vec2.x;
		float myY = p2_s.y + slope*(p1_s.x - p2_s.x);
		// float scale = fabs(diff_s.x / vec2.x);
		// newPt.y = p2_s.y + scale * vec2.y;
		newPt.y = myY;
		return newPt;
	}
	if (fabs(vec2.x) < 0.0001) {
		Point newPt = Point(p2_s.x, 0, p2_s.z);
		float slope = vec1.y / vec1.x;
		float myY = p1_s.y + slope*(p2_s.x - p1_s.x);
		// float scale = fabs(diff_s.x / vec1.x);
		// newPt.y = p1_s.y + scale * vec1.y;
		newPt.y = myY;
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
