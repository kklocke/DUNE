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


    ////////////////////////////////////////////////////////////////////////////
    //                           Single Time Slice                            //
    ////////////////////////////////////////////////////////////////////////////

    // ofstream pathFile;
    // pathFile.open("cellTest_path.txt");
    // pathFile << myEvent.vertex.x << " " << endPath.x << "\n";
    // pathFile << myEvent.vertex.y << " " << endPath.y << "\n";
    // pathFile.close();
    //
    // vector <float> signals = testLayer.signalVec(myEvent.path1);
    //
    //
    // cout << "MySignal: " << endl;
    // for (int i = 0; i < (int)signals.size(); i++) {
    //     cout << signals[i] << " ";
    // }
    // cout << endl;
    //
    //
    // int vDim = testLayer.vertical.startPoints.size();
    // int lDim = testLayer.lSlant.startPoints.size();
    // int rDim = testLayer.rSlant.startPoints.size();
    //
    // ofstream wiresFile;
    // wiresFile.open("cellTest_wires.txt");
    // vector <float> mySig;
    // for (int i = 0; i < vDim - 1; i++) {
    //     if (signals[i] > 0) {
    //         mySig.push_back(signals[i]);
    //         wiresFile << testLayer.vertical.startPoints[i].x << " " << testLayer.vertical.endPoints[i].x << "\n";
    //         wiresFile << testLayer.vertical.startPoints[i].y << " " << testLayer.vertical.endPoints[i].y << "\n";
    //         wiresFile << testLayer.vertical.startPoints[i + 1].x << " " << testLayer.vertical.endPoints[i + 1].x << "\n";
    //         wiresFile << testLayer.vertical.startPoints[i + 1].y << " " << testLayer.vertical.endPoints[i + 1].y << "\n";
    //
    //     }
    // }
    // for (int i = 0; i < lDim - 1; i++) {
    //     if (signals[i + vDim - 1] > 0) {
    //         mySig.push_back(signals[i + vDim - 1]);
    //         wiresFile << testLayer.lSlant.startPoints[i].x << " " << testLayer.lSlant.endPoints[i].x << "\n";
    //         wiresFile << testLayer.lSlant.startPoints[i].y << " " << testLayer.lSlant.endPoints[i].y << "\n";
    //         wiresFile << testLayer.lSlant.startPoints[i + 1].x << " " << testLayer.lSlant.endPoints[i + 1].x << "\n";
    //         wiresFile << testLayer.lSlant.startPoints[i + 1].y << " " << testLayer.lSlant.endPoints[i + 1].y << "\n";
    //     }
    // }
    // for (int i = 0; i < rDim - 1; i++) {
    //     if (signals[i + vDim + lDim - 2] > 0) {
    //         mySig.push_back(signals[i + vDim + lDim - 2]);
    //         wiresFile << testLayer.rSlant.startPoints[i].x << " " << testLayer.rSlant.endPoints[i].x << "\n";
    //         wiresFile << testLayer.rSlant.startPoints[i].y << " " << testLayer.rSlant.endPoints[i].y << "\n";
    //         wiresFile << testLayer.rSlant.startPoints[i + 1].x << " " << testLayer.rSlant.endPoints[i + 1].x << "\n";
    //         wiresFile << testLayer.rSlant.startPoints[i + 1].y << " " << testLayer.rSlant.endPoints[i + 1].y << "\n";
    //     }
    // }
    // wiresFile.close();
    //
    // vector <vector <int>> myMatrix = testLayer.geoMatrix(signals);
    //
    // cout << "mySig size: " << mySig.size() << endl;
    // cout << "Num cells: " << testLayer.cells.size() << endl;
    // cout << "Geomatrix length: " << myMatrix[0].size() << endl;
    // cout << "Geomatrix height: " << myMatrix.size() << endl;
    // cout << "Sig Size: " << signals.size() << endl;
    //
    // vector <float> actualCharge = path2true(myEvent.path1, testLayer);
    //
    // // vector <float> trueSig = solveTrue(signals, myMatrix);
    //
    // vector <float> geneticSig = solveCharge_genetic(signals, myMatrix, testLayer);
    //
    // // cout << "\nTrue sig size: " << trueSig.size() << endl;
    //
    // // cout << "GEOMATRIX: " << endl;
    // // for (int i = 0; i < (int)myMatrix.size(); i++) {
    // //     for (int j = 0; j < (int)myMatrix[i].size(); j++) {
    // //         cout << myMatrix[i][j] << " ";
    // //     }
    // //     cout << "\n";
    // // }
    //
    // // vector <float> optSoln = solveCharge(signals, myMatrix, testLayer);
    //
    // ofstream cellFile;
    // cellFile.open("cellTest_cells.txt");
    // for (int i = 0; i < (int)testLayer.cells.size(); i++) {
    //     for (int j = 0; j < (int)testLayer.cells[i].vertices.size(); j++) {
    //         cellFile << testLayer.cells[i].vertices[j].x << " " << testLayer.cells[i].vertices[j].y << " ";
    //     }
    //     cellFile << "\n";
    //     cellFile << geneticSig[i] << " " << actualCharge[i] << "\n";
    //     // cellFile << trueSig[i] << " " << actualCharge[i] << "\n";
    // }
    // cellFile.close();
    //




    ////////////////////////////////////////////////////////////////////////////
    //                      Path Partitioned Version                          //
    ////////////////////////////////////////////////////////////////////////////
    vector <vector <Path>> myPartEvent = myEvent.partitionEvent(20, 200);
    for (int i = 0; i < (int)myPartEvent.size(); i++) {
        for (int j = 0; j < (int)myPartEvent[i].size(); j++) {
            cout << "Partitioned Path part " << j << endl;
            cout << "\t";
            myPartEvent[i][j].vertex.printPt();
            cout << " --> ";
            Point endPt = myPartEvent[i][j].vertex + myPartEvent[i][j].path2vec();
            endPt.printPt();
            cout << endl;
        }
    }
    cout << "\n";

    vector <wireLayer> allLayers;
    vector <vector <vector <int>>> allGeoMats;
    vector <vector <float>> allSigs;
    vector <vector <float>> allActual;

    // Save the partitioned path
    ofstream pathFile;
    pathFile.open("cellTest_path.txt");
    ofstream wiresFile;
    wiresFile.open("cellTest_wires.txt");

    for (int i = 0; i < (int)myPartEvent[0].size(); i++) {
        Point endPt = myPartEvent[0][i].vertex + myPartEvent[0][i].path2vec();
        pathFile << myPartEvent[0][i].vertex.x << "\t" << endPt.x << "\n";
        pathFile << myPartEvent[0][i].vertex.y << "\t" << endPt.y << "\n";
        vector <float> signals = testLayer.signalVec(myPartEvent[0][i]);
        cout << "Signal vec for partition " << i << endl;
        for (int j = 0; j < (int)signals.size(); j++) {
            cout << signals[j] << " ";
        }
        cout << endl;
        allSigs.push_back(signals);
        int vDim = testLayer.vertical.startPoints.size();
        int lDim = testLayer.lSlant.startPoints.size();
        int rDim = testLayer.rSlant.startPoints.size();
        int chargeCount = 0;
        for (int j = 0; j < vDim - 1; j++) {
            if (signals[j] > 0) {
                chargeCount++;
                wiresFile << testLayer.vertical.startPoints[j].x << " " << testLayer.vertical.endPoints[j].x << "\n";
                wiresFile << testLayer.vertical.startPoints[j].y << " " << testLayer.vertical.endPoints[j].y << "\n";
                wiresFile << testLayer.vertical.startPoints[j + 1].x << " " << testLayer.vertical.endPoints[j + 1].x << "\n";
                wiresFile << testLayer.vertical.startPoints[j + 1].y << " " << testLayer.vertical.endPoints[j + 1].y << "\n";
            }
        }
        for (int j = 0; j < lDim - 1; j++) {
            if (signals[j + vDim - 1] > 0) {
                chargeCount++;
                wiresFile << testLayer.lSlant.startPoints[j].x << " " << testLayer.lSlant.endPoints[j].x << "\n";
                wiresFile << testLayer.lSlant.startPoints[j].y << " " << testLayer.lSlant.endPoints[j].y << "\n";
                wiresFile << testLayer.lSlant.startPoints[j + 1].x << " " << testLayer.lSlant.endPoints[j + 1].x << "\n";
                wiresFile << testLayer.lSlant.startPoints[j + 1].y << " " << testLayer.lSlant.endPoints[j + 1].y << "\n";
            }
        }
        for (int j = 0; j < rDim - 1; j++) {
            if (signals[j + vDim + lDim - 2] > 0) {
                chargeCount++;
                wiresFile << testLayer.rSlant.startPoints[j].x << " " << testLayer.rSlant.endPoints[j].x << "\n";
                wiresFile << testLayer.rSlant.startPoints[j].y << " " << testLayer.rSlant.endPoints[j].y << "\n";
                wiresFile << testLayer.rSlant.startPoints[j + 1].x << " " << testLayer.rSlant.endPoints[j + 1].x << "\n";
                wiresFile << testLayer.rSlant.startPoints[j + 1].y << " " << testLayer.rSlant.endPoints[j + 1].y << "\n";
            }
        }
        cout << "Number of charged wires in slice " << i <<": " << chargeCount << endl;
        // finished writing the wires, generate the geomatrix and push
        wireLayer tempLayer = testLayer;
        vector <vector <int>> myGeoMat = tempLayer.geoMatrix(signals);
        vector <float> actualCharge = path2true(myPartEvent[0][i], tempLayer);
        allActual.push_back(actualCharge);
        allLayers.push_back(tempLayer);
        allGeoMats.push_back(myGeoMat);
    }
    pathFile.close();
    wiresFile.close();

    // vector <vector <float>> solvedCharges = solveChargeMulti(allSigs, allGeoMats, allLayers);
    vector <vector <float>> solvedCharges = solveChargeMulti_genetic(allSigs, allGeoMats, allLayers);

    cout << "SOLVED CHARGES MULTI: " << endl;
    for (int i = 0; i < (int)solvedCharges.size(); i++) {
        cout << "PARTITION: " << i << endl;
        for (int j = 0; j < (int)solvedCharges[i].size(); j++) {
            cout << solvedCharges[i][j] << " ";
        }
        cout << "\n";
    }

    ofstream cellFile;
    cellFile.open("cellTest_cells.txt");
    for (int i = 0; i < (int)allLayers.size(); i++) {
        for (int j = 0; j < (int)allLayers[i].cells.size(); j++) {
            for (int k = 0; k < (int)allLayers[i].cells[j].vertices.size(); k++) {
                cellFile << allLayers[i].cells[j].vertices[k].x << " " << allLayers[i].cells[j].vertices[k].y << " ";
            }
            cellFile << "\n";
            cellFile << solvedCharges[i][j] << " " << allActual[i][j] << "\n";
        }
    }
    cellFile.close();

    return 0;
}
