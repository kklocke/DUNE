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
    int dims [3] = {200, 200, 200};
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
    wireArray test = wireArray(5., 60., 200., 200., 0.);
    wireLayer testLayer = wireLayer(200., 200., 0., 5.);
    Point endPath = myEvent.vertex + myEvent.path1.path2vec();
    cout << "My end path: ";
    endPath.printPt();
    cout << endl;

	ofstream myfile;
	myfile.open("opt_vertex_wires.txt");
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

	ofstream sigFile;
	sigFile.open("opt_hitlines.txt");
	vector <float> signals = testLayer.signalVec(myEvent.path1);
	vector <float> mySig;
	for (int i = 0; i < (int)signals.size(); i++) {
		if (signals[i] > 0) {
			mySig.push_back(signals[i]);
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
	vector <vector <int>> myMatrix = testLayer.geoMatrix_test(signals);
	cout << myMatrix.size();
	vector <Cell> allCells = tiling(testLayer, 200., 200., testLayer.signalVec(myEvent.path1));
	vector <float> trueSig = solveCharge(mySig, myMatrix, testLayer, allCells);
	cout << "True Signal\n";
	ofstream pointsFile;
	pointsFile.open("opt_displayPoints.txt");

	cout << "TRUE SIG SIZE: " << trueSig.size() << endl;
	cout << "MY SIG SIZE: " << mySig.size() << endl;
	cout << "SIGNALS SIZE: " << signals.size() << endl;
	cout << "GRIDS SIZE: " << testLayer.grids.size() << endl;
	cout << "GEOMAT DIMS: " << myMatrix.size() << "\t" <<  myMatrix[0].size() << endl;

	for (int i = 0; i < (int)trueSig.size(); i++) {
		if (trueSig[i] > 3) {
			pointsFile << testLayer.grids[i].x << "\t" << testLayer.grids[i].y << "\t" << trueSig[i] << endl;
			cout << trueSig[i] << " ";
		}
	}
	pointsFile.close();
	cout << endl;




	// find all the crossings
	// construct the geoMatrix
	// solve true
	// display the optimum in a heatmap style

	return 0;
}
