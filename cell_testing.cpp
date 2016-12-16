#include "wires_rewrite.hh"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

int main () {
    srand(time(NULL));

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
	wireArray test = wireArray(3., 54.3, 200., 200., 0.);
	wireLayer testLayer = wireLayer(200., 200., 0., 3.);
	Point endPath = myEvent.vertex + myEvent.path1.path2vec();
	cout << "My end path: ";
	endPath.printPt();
	cout << endl;

    ofstream pathFile;
    pathFile.open("cellTest_path.txt");
    pathFile << myEvent.vertex.x << " " << endPath.x << "\n";
    pathFile << myEvent.vertex.y << " " << endPath.y << "\n";
    pathFile.close();

    vector <float> signals = testLayer.signalVec(myEvent.path1);


    cout << "MySignal: " << endl;
    for (int i = 0; i < (int)signals.size(); i++) {
        cout << signals[i] << " ";
    }
    cout << endl;


    int vDim = testLayer.vertical.startPoints.size();
    int lDim = testLayer.lSlant.startPoints.size();
    int rDim = testLayer.rSlant.startPoints.size();

    ofstream wiresFile;
    wiresFile.open("cellTest_wires.txt");
    vector <float> mySig;
    for (int i = 0; i < vDim - 1; i++) {
        if (signals[i] > 0) {
            mySig.push_back(signals[i]);
            wiresFile << testLayer.vertical.startPoints[i].x << " " << testLayer.vertical.endPoints[i].x << "\n";
            wiresFile << testLayer.vertical.startPoints[i].y << " " << testLayer.vertical.endPoints[i].y << "\n";
        }
    }
    for (int i = 0; i < lDim - 1; i++) {
        if (signals[i + vDim - 1] > 0) {
            mySig.push_back(signals[i + vDim - 1]);
            wiresFile << testLayer.lSlant.startPoints[i].x << " " << testLayer.lSlant.endPoints[i].x << "\n";
            wiresFile << testLayer.lSlant.startPoints[i].y << " " << testLayer.lSlant.endPoints[i].y << "\n";
        }
    }
    for (int i = 0; i < rDim - 1; i++) {
        if (signals[i + vDim + lDim - 2] > 0) {
            mySig.push_back(signals[i + vDim + lDim - 2]);
            wiresFile << testLayer.rSlant.startPoints[i].x << " " << testLayer.rSlant.endPoints[i].x << "\n";
            wiresFile << testLayer.rSlant.startPoints[i].y << " " << testLayer.rSlant.endPoints[i].y << "\n";
        }
    }
    wiresFile.close();

    vector <vector <int>> myMatrix = testLayer.geoMatrix(signals);

    cout << "mySig size: " << mySig.size() << endl;
    cout << "Num cells: " << testLayer.cells.size() << endl;
    cout << "Geomatrix length: " << myMatrix[0].size() << endl;
    cout << "Geomatrix height: " << myMatrix.size() << endl;
    cout << "Sig Size: " << signals.size() << endl;

    vector <float> actualCharge = path2true(myEvent.path1, testLayer);

    vector <float> trueSig = solveTrue(signals, myMatrix);

    cout << "\nTrue sig size: " << trueSig.size() << endl;

    // cout << "GEOMATRIX: " << endl;
    // for (int i = 0; i < (int)myMatrix.size(); i++) {
    //     for (int j = 0; j < (int)myMatrix[i].size(); j++) {
    //         cout << myMatrix[i][j] << " ";
    //     }
    //     cout << "\n";
    // }

    vector <float> optSoln = solveCharge(signals, myMatrix, testLayer);

    ofstream cellFile;
    cellFile.open("cellTest_cells.txt");
    for (int i = 0; i < (int)testLayer.cells.size(); i++) {
        for (int j = 0; j < (int)testLayer.cells[i].vertices.size(); j++) {
            cellFile << testLayer.cells[i].vertices[j].x << " " << testLayer.cells[i].vertices[j].y << " ";
        }
        cellFile << "\n";
        cellFile << optSoln[i] << " " << actualCharge[i] << "\n";
        // cellFile << trueSig[i] << "\n";
    }
    cellFile.close();





    return 0;
}
