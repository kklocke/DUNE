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

    ofstream layerFile;
    layerFile.open("multiOpt_layer.txt");
    for (int i = 0; i < (int)testLayer.vertical.startPoints.size(); i++) {
        layerFile << testLayer.vertical.startPoints[i].x << "\t" << testLayer.vertical.startPoints[i].y;
        layerFile << "\t" << testLayer.vertical.endPoints[i].x << "\t" << testLayer.vertical.endPoints[i].y << "\n";
    }
    for (int i = 0; i < (int)testLayer.lSlant.startPoints.size(); i++) {
        layerFile << testLayer.lSlant.startPoints[i].x << "\t" << testLayer.lSlant.startPoints[i].y;
        layerFile << "\t" << testLayer.lSlant.endPoints[i].x << "\t" << testLayer.lSlant.endPoints[i].y << "\n";
    }
    for (int i = 0; i < (int)testLayer.rSlant.startPoints.size(); i++) {
        layerFile << testLayer.rSlant.startPoints[i].x << "\t" << testLayer.rSlant.startPoints[i].y;
        layerFile << "\t" << testLayer.rSlant.endPoints[i].x << "\t" << testLayer.rSlant.endPoints[i].y << "\n";
    }
    layerFile.close();
    // project the event
    vector <vector <Path>> myPartEvent = myEvent.partitionEvent(20, 200);
    ofstream myFile;
    myFile.open("multiOpt_path.txt");
    cout << "\nPartitioned End Points\n";
    for (int i = 0;  i < (int)myPartEvent[0].size(); i++) {
        Point endPt = myPartEvent[0][i].vertex + myPartEvent[0][i].path2vec();
        endPt.printPt();
        cout << "\n";
        myFile << myPartEvent[0][i].vertex.x << "\t" << myPartEvent[0][i].vertex.y << "\t";
        myFile << endPath.x << "\t" << endPath.y << "\n";
    }
    myFile.close();

    vector <wireLayer> allLayers;

    // For each partition:
    for (int i = 0; i < (int)myPartEvent[0].size(); i++) {
        Point endPath = myPartEvent[0][i].vertex + myPartEvent[0][i].path2vec();
        vector <vector <int>> layerCross = testLayer.allCrossings(myPartEvent[0][i].vertex, endPath);
        wireLayer tempLayer = testLayer;
        allLayers.push_back(tempLayer);
        // Save these to a file
    }

    vector <vector <float>> allSignals;
    vector <vector <float>> allMySignals;
    vector <vector <vector <int>>> allGeoMats;
    vector <vector <Cell>> allTiles;
    for (int i = 0; i < (int)myPartEvent[0].size(); i++) {
        vector <float> signals = testLayer.signalVec(myPartEvent[0][i]);
        vector <float> mySig;
        for (int j = 0; j < (int)signals.size(); j++) {
            if (signals[j] > 0) {
                mySig.push_back(signals[j]);
                // The stuff about saving to a file
            }
        }
        vector <vector <int>> myMatrix = allLayers[i].geoMatrix_test(signals);
        allSignals.push_back(signals);
        allMySignals.push_back(mySig);
        allGeoMats.push_back(myMatrix);
        vector <Cell> myTiles = tiling(allLayers[i], 200., 200., allLayers[i].signalVec(myPartEvent[0][i]));
        vector <Blob> myBlobs = tile2blob(myTiles);
        allTiles.push_back(myBlobs[0].cells);
    }
    vector <vector <float>> solvedCharges = solveChargeMulti(allMySignals, allGeoMats, allLayers, allTiles);

    ofstream chargeFile;
    chargeFile.open("multiOpt_charges.txt");
    for (int i = 0; i < (int)solvedCharges.size(); i++) {
        for (int j = 0; j < (int)solvedCharges[i].size(); j++) {
            if (solvedCharges[i][j] > 3) {
                chargeFile << allLayers[i].grids[j].x << "\t" << allLayers[i].grids[j].y << "\t";
                chargeFile << solvedCharges[i][j] << "\n";
            }
        }
    }
    chargeFile.close();
    return 0;
}
