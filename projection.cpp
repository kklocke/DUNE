#include "events.hh"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

int main () {
	srand(time(NULL));

	// create a random event
	// project it --> time slices
	// figure out how to dispay it
	// figure out how to apply the current method to it.
	int dims [3] = {50, 50, 50};
	NaiveEvent myEvent = randEvent(dims);
	// print vertex
	cout << "Vertex: " << myEvent.vertex.x << "\t" << myEvent.vertex.y << "\t" << myEvent.vertex.z << endl;
	cout << "Path 1:" << endl;
	cout << "\tLength: " << myEvent.path1.length << endl;
	cout << "\tTheta: " << myEvent.path1.theta << endl;
	cout << "\tPhi: " << myEvent.path1.phi << endl;
	cout << "Path 2:" << endl;
	cout << "\tLength: " << myEvent.path2.length << endl;
	cout << "\tTheta: " << myEvent.path2.theta << endl;
	cout << "\tPhi: " << myEvent.path2.phi << endl;
	// print path1
	// print path2
	wireArray test = wireArray(3., 45., 50., 50., 0.);
	for (int i = 0; i < int(test.startPoints.size()); i++) {
		cout << "Wire " << i << endl;
		cout << "\t";
		test.startPoints[i].printPt();
		cout << "\n\t";
		test.endPoints[i].printPt();
		cout << "\n";
	}
	Point endPath = myEvent.vertex + myEvent.path1.path2vec();
	cout << "My end path: ";
	endPath.printPt();

	cout << "\nLines crossed by the first path from my event" << endl;

	vector <int> myCrosses = crossings(test, myEvent.vertex, endPath);
	for (int i = 0; i < (int)myCrosses.size(); i++) {
		cout << myCrosses[i] << endl;
	}

	ofstream myfile;
	myfile.open("displayTest.txt");
	myfile << myEvent.vertex.x << "\t" << myEvent.vertex.y << "\n";
	myfile << endPath.x << "\t" << endPath.y << "\n";
	for (int i = 0; i < (int)test.startPoints.size(); i++)
	{
		myfile << test.startPoints[i].x << "\t" << test.startPoints[i].y << "\n";
		myfile << test.endPoints[i].x << "\t" << test.endPoints[i].y << "\n";
	}
	myfile.close();
	return 0;
}
