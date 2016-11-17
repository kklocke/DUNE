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
	wireArray test = wireArray(3., 45., 200., 200., 0.);
	wireLayer testLayer = wireLayer(200., 200., 0., 3.);
	Point endPath = myEvent.vertex + myEvent.path1.path2vec();
	cout << "My end path: ";
	endPath.printPt();
	cout << endl;

	// cout << "PARTITIONING PATH\n";
	// vector <vector <Path>> myPartEvent = myEvent.partitionEvent(20, 200);
	// for (int i = 0; i < (int)myPartEvent.size(); i++) {
	// 	for (int j = 0; j < (int)myPartEvent[i].size(); j++) {
	// 		cout << "Vertex: " << myPartEvent[i][j].vertex.x << "\t" << myPartEvent[i][j].vertex.y << "\t" << myPartEvent[i][j].vertex.z << endl;
	// 		cout << "\tLength: " << myPartEvent[i][j].length << endl;
	// 		cout << "\tTheta: " << myPartEvent[i][j].theta << endl;
	// 		cout << "\tPhi: " << myPartEvent[i][j].phi << endl;
	// 		Point myEndPt = myPartEvent[i][j].vertex + myPartEvent[i][j].path2vec();
	// 		cout << "\tEnd Point: ";
	// 		myEndPt.printPt();
	// 		cout << endl;
	// 	}
	// }
	//
	// vector <vector <int>> allCross = testLayer.allCrossings(myEvent.vertex, endPath);
	//
	// // vector <vector <Path>> myPaths = myEvent.partitionEvent(5, 20);
	// // for (int i = 0; i < (int)myPaths[0].size(); i++) {
	// // 	cout << "Path: ";
	// // 	myPaths[0][i].vertex.printPt();
	// // 	cout << "\t";
	// // 	(myPaths[0][i].vertex + myPaths[0][i].path2vec()).printPt();
	// // 	cout << "\n";
	// // }
	//
	// ofstream myfile;
	// myfile.open("displayTest.txt");
	// myfile << myEvent.vertex.x << "\t" << myEvent.vertex.y << "\n";
	// myfile << endPath.x << "\t" << endPath.y << "\n";
	// for (int i = 0; i < (int)testLayer.vertical.startPoints.size(); i++)
	// {
	// 	myfile << testLayer.vertical.startPoints[i].x << "\t" << testLayer.vertical.startPoints[i].y << "\n";
	// 	myfile << testLayer.vertical.endPoints[i].x << "\t" << testLayer.vertical.endPoints[i].y << "\n";
	// }
	// for (int i = 0; i < (int)testLayer.lSlant.startPoints.size(); i++)
	// {
	// 	myfile << testLayer.lSlant.startPoints[i].x << "\t" << testLayer.lSlant.startPoints[i].y << "\n";
	// 	myfile << testLayer.lSlant.endPoints[i].x << "\t" << testLayer.lSlant.endPoints[i].y << "\n";
	// }
	// for (int i = 0; i < (int)testLayer.rSlant.startPoints.size(); i++)
	// {
	// 	myfile << testLayer.rSlant.startPoints[i].x << "\t" << testLayer.rSlant.startPoints[i].y << "\n";
	// 	myfile << testLayer.rSlant.endPoints[i].x << "\t" << testLayer.rSlant.endPoints[i].y << "\n";
	// }
	// myfile.close();
	//
	// // for (int i = 0; i < (int)myMatrix.size(); i++) {
	// // 	for (int j = 0; j < (int)myMatrix[i].size(); j++) {
	// // 		if (myMatrix[i][j] == 1) {
	// // 			cout << "0" << " ";
	// // 		}
	// // 		else {
	// // 			cout << "- ";
	// // 		}
	// // 	}
	// // 	cout << endl;
	// // }
	// ofstream sigFile;
	// sigFile.open("hitLines.txt");
	// vector <float> signals = testLayer.signalVec(myEvent.path1);
	// vector <float> mySig;
	// for (int i = 0; i < (int)signals.size(); i++) {
	// 	if (signals[i] > 0) {
	// 		mySig.push_back(signals[i]);
	// 		if (i < (int)testLayer.vertical.startPoints.size()) {
	// 			sigFile << testLayer.vertical.startPoints[i].x << "\t" << testLayer.vertical.startPoints[i].y << "\n";
	// 			sigFile << testLayer.vertical.endPoints[i].x << "\t" << testLayer.vertical.endPoints[i].y << "\n";
	// 		}
	// 		else if (i < (int)testLayer.vertical.startPoints.size() + (int)testLayer.lSlant.startPoints.size()) {
	// 			int j = i - (int)testLayer.vertical.startPoints.size();
	// 			sigFile << testLayer.lSlant.startPoints[j].x << "\t" << testLayer.lSlant.startPoints[j].y << "\n";
	// 			sigFile << testLayer.lSlant.endPoints[j].x << "\t" << testLayer.lSlant.endPoints[j].y << "\n";
	// 		}
	// 		else {
	// 			int j = i - (int)testLayer.vertical.startPoints.size() - (int)testLayer.lSlant.startPoints.size();
	// 			sigFile << testLayer.rSlant.startPoints[j].x << "\t" << testLayer.rSlant.startPoints[j].y << "\n";
	// 			sigFile << testLayer.rSlant.endPoints[j].x << "\t" << testLayer.rSlant.endPoints[j].y << "\n";
	// 		}
	// 	}
	// 	// cout << signals[i] << " ";
	// }
	// sigFile.close();
	// cout << endl;
	// vector <vector <int>> myMatrix = testLayer.geoMatrix_test(signals);
	// cout << myMatrix.size();
	// vector <float> trueSig = solveTrue(mySig, myMatrix);
	// // cout << "True Signal\n";
	// ofstream pointsFile;
	// pointsFile.open("displayPoints.txt");
	// for (int i = 0; i < (int)trueSig.size(); i++) {
	// 	if (trueSig[i] > 1) {
	// 		pointsFile << testLayer.grids[i].x << "\t" << testLayer.grids[i].y << "\t" << trueSig[i] << endl;
	// 		cout << trueSig[i] << " ";
	// 	}
	// }
	// pointsFile.close();
	// cout << endl;

	vector <vector <Path>> myPartEvent = myEvent.partitionEvent(20, 200);
	for (int i = 0; i < (int)myPartEvent.size(); i++) {
	    for (int j = 0; j < (int)myPartEvent[i].size(); j++) {
	        Point myEndPt = myPartEvent[i][j].vertex + myPartEvent[i][j].path2vec();
	        myEndPt.printPt();
	    }
	}

	ofstream myfile;
	myfile.open("displayTest.txt");
	myfile << myPartEvent[0].size() << "\n";
	for (int i = 0; i < (int)myPartEvent[0].size(); i++) {
	    Point endPath = myPartEvent[0][i].vertex + myPartEvent[0][i].path2vec();
	    vector <vector <int>> layerCross = testLayer.allCrossings(myPartEvent[0][i].vertex, endPath);
	    myfile << myPartEvent[0][i].vertex.x << "\t" << myPartEvent[0][i].vertex.y << "\n";
	    myfile << endPath.x << "\t" << endPath.y << "\n";
	}
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
	sigFile.open("hitLines.txt");
	ofstream pointsFile;
	pointsFile.open("displayPoints.txt");
	cout << "\n\nNUMBER OF PATHS: " << myPartEvent[0].size() << endl;
	for (int i = 0; i < (int)myPartEvent[0].size(); i++) {
	    vector <float> signals = testLayer.signalVec(myPartEvent[0][i]);
	    vector <float> mySig;
	    for (int j = 0; j < (int)signals.size(); j++) {
	        if (signals[j] > 0) {
	            mySig.push_back(signals[j]);
	            if (j < (int)testLayer.vertical.startPoints.size()) {
	                sigFile << testLayer.vertical.startPoints[j].x << "\t" << testLayer.vertical.startPoints[j].y << "\n";
	                sigFile << testLayer.vertical.endPoints[j].x << "\t" << testLayer.vertical.endPoints[j].y << "\n";
	            }
	            else if (j < (int)testLayer.vertical.startPoints.size() + (int)testLayer.lSlant.startPoints.size()) {
	                int k = j - (int)testLayer.vertical.startPoints.size();
	                sigFile << testLayer.lSlant.startPoints[k].x << "\t" << testLayer.lSlant.startPoints[k].y << "\n";
	                sigFile << testLayer.lSlant.endPoints[k].x << "\t" << testLayer.lSlant.endPoints[k].y << "\n";
	            }
	            else {
	                int k  = j - (int)testLayer.vertical.startPoints.size() - (int)testLayer.lSlant.startPoints.size();
	                sigFile << testLayer.rSlant.startPoints[k].x << "\t" << testLayer.rSlant.startPoints[k].y << "\n";
	                sigFile << testLayer.rSlant.endPoints[k].x << "\t" << testLayer.rSlant.endPoints[k].y << "\n";
	            }
	        }
	    }
	    sigFile << "NEXT PARTITION\n";
	    vector <vector <int>> myMatrix = testLayer.geoMatrix_test(signals);
	    vector <float> trueSig = solveTrue(mySig, myMatrix);
	    for (int j = 0; j < (int)trueSig.size(); j++) {
	        if (trueSig[i] > 1) {
	            pointsFile << testLayer.grids[j].x << "\t" << testLayer.grids[j].y << "\t" << trueSig[j] << endl;
	        }
	    }
		cout << "\n";
	    pointsFile << "NEXT PARTITION\n";
	}
	pointsFile.close();
	sigFile.close();

	cout << endl;




	ofstream cellFile;
	cellFile.open("cellVertices.txt");

	vector <Cell> allCells = tiling(testLayer, 200., 200.);
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
