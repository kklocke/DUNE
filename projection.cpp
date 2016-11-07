#include "wires.hh"
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
	int dims [3] = {100, 100, 50};
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
	wireArray test = wireArray(3., 45., 100., 100., 0.);
	wireLayer testLayer = wireLayer(100., 100., 0., 3.);
	Point endPath = myEvent.vertex + myEvent.path1.path2vec();
	cout << "My end path: ";
	endPath.printPt();
	cout << endl;


	vector <vector <int>> allCross = testLayer.allCrossings(myEvent.vertex, endPath);

	// vector <vector <Path>> myPaths = myEvent.partitionEvent(5, 20);
	// for (int i = 0; i < (int)myPaths[0].size(); i++) {
	// 	cout << "Path: ";
	// 	myPaths[0][i].vertex.printPt();
	// 	cout << "\t";
	// 	(myPaths[0][i].vertex + myPaths[0][i].path2vec()).printPt();
	// 	cout << "\n";
	// }

	ofstream myfile;
	myfile.open("displayTest.txt");
	myfile << myEvent.vertex.x << "\t" << myEvent.vertex.y << "\n";
	myfile << endPath.x << "\t" << endPath.y << "\n";
	for (int i = 0; i < (int)testLayer.vertical.startPoints.size(); i++)
	{
		myfile << testLayer.vertical.startPoints[i].x << "\t" << testLayer.vertical.startPoints[i].y << "\n";
		myfile << testLayer.vertical.endPoints[i].x << "\t" << testLayer.vertical.endPoints[i].y << "\n";
	}
	for (int i = 0; i < (int)testLayer.lSlant.startPoints.size(); i++)
	{
		myfile << testLayer.lSlant.startPoints[i].x << "\t" << testLayer.lSlant.startPoints[i].y << "\n";
		myfile << testLayer.lSlant.endPoints[i].x << "\t" << testLayer.lSlant.endPoints[i].y << "\n";
	}
	for (int i = 0; i < (int)testLayer.rSlant.startPoints.size(); i++)
	{
		myfile << testLayer.rSlant.startPoints[i].x << "\t" << testLayer.rSlant.startPoints[i].y << "\n";
		myfile << testLayer.rSlant.endPoints[i].x << "\t" << testLayer.rSlant.endPoints[i].y << "\n";
	}
	myfile.close();

	vector <vector <int>> myMatrix = testLayer.geoMatrix();
	// for (int i = 0; i < (int)myMatrix.size(); i++) {
	// 	for (int j = 0; j < (int)myMatrix[i].size(); j++) {
	// 		if (myMatrix[i][j] == 1) {
	// 			cout << "0" << " ";
	// 		}
	// 		else {
	// 			cout << "- ";
	// 		}
	// 	}
	// 	cout << endl;
	// }
	ofstream sigFile;
	sigFile.open("hitLines.txt");
	vector <float> signals = testLayer.signalVec(myEvent.path1);
	for (int i = 0; i < (int)signals.size(); i++) {
		if (signals[i] > 0) {
			if (i < (int)testLayer.vertical.startPoints.size()) {
				sigFile << testLayer.vertical.startPoints[i].x << "\t" << testLayer.vertical.startPoints[i].y << "\n";
				sigFile << testLayer.vertical.endPoints[i].x << "\t" << testLayer.vertical.endPoints[i].y << "\n";
			}
			else if (i < (int)testLayer.vertical.startPoints.size() + (int)testLayer.lSlant.startPoints.size()) {
				int j = i - (int)testLayer.vertical.startPoints.size();
				sigFile << testLayer.lSlant.startPoints[j].x << "\t" << testLayer.lSlant.startPoints[j].y << "\n";
				sigFile << testLayer.lSlant.endPoints[j].x << "\t" << testLayer.lSlant.endPoints[j].y << "\n";
			}
			else {
				int j = i - (int)testLayer.vertical.startPoints.size() - (int)testLayer.lSlant.startPoints.size();
				sigFile << testLayer.rSlant.startPoints[j].x << "\t" << testLayer.rSlant.startPoints[j].y << "\n";
				sigFile << testLayer.rSlant.endPoints[j].x << "\t" << testLayer.rSlant.endPoints[j].y << "\n";
			}
		}
		// cout << signals[i] << " ";
	}
	sigFile.close();
	cout << endl;
	vector <float> trueSig = solveTrue(signals, myMatrix);
	// cout << "True Signal\n";
	ofstream pointsFile;
	pointsFile.open("displayPoints.txt");
	for (int i = 0; i < (int)trueSig.size(); i++) {
		if (trueSig[i] > 0.1) {
			pointsFile << testLayer.grids[i].x << "\t" << testLayer.grids[i].y << endl;
			cout << trueSig[i] << " ";
		}
	}
	pointsFile.close();
	cout << endl;

	ofstream cellFile;
	cellFile.open("cellVertices.txt");

	vector <Cell> allCells = tiling(testLayer, 100., 100.);
	for (int i = 0; i < (int)allCells.size(); i++) {
		for (int j = 0; j < (int)allCells[i].vertices.size(); j++) {
			cellFile << allCells[i].vertices[j].x << "\t";
		}
		cellFile << "\n";
		for (int j = 0; j < (int)allCells[i].vertices.size(); j++) {
			cellFile << allCells[i].vertices[j].y << "\t";
		}
		cellFile << "\n";
	}
	cellFile.close();
	return 0;
}
