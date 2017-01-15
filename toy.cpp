#include "wires_rewrite.hh"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

int main (int argc, char *argv[]) {

    // -m : sets multi-slice mode (single slice is default)
    // -d : sets double path mode (single path is default)

    bool singleSlice = true;
    bool singlePath = true;

    srand(time(NULL));

    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            if (strcmp(argv[i], "-m") == 0) {
                singleSlice = false;
            }
            else if (strcmp(argv[i], "-d") == 0) {
                singlePath = false;
            }
            else {
                cout << "Usage: ./toy [-m] [-d]\n";
                exit(1);
            }
        }
    }

    // Set the detector size
    int dims [3] = {200, 200, 200};

    // Generate a random event consisting of two paths emanating from a common vetex
    NaiveEvent myEvent = randEvent(dims);

    // Display information about the event that was generated
    cout << "////////////////////////////////////////////////////////////\n";
    cout << "//                   Event  Information                   //\n";
    cout << "////////////////////////////////////////////////////////////\n";

    cout << "\tVertex: ";
    myEvent.vertex.printPt();
    cout << "\n\tPath1: \n";
    cout << "\t\tLength: " << myEvent.path1.length << endl;
    cout << "\t\tTheta: " << myEvent.path1.theta << endl;
    cout << "\t\tPhi: " << myEvent.path1.phi << endl;
    cout << "\t\tEnd Point: ";
    Point path1End = myEvent.path1.vertex + myEvent.path1.path2vec();
    path1End.printPt();
    cout << "\n";
    cout << "\n\tPath2: \n";
    cout << "\t\tLength: " << myEvent.path2.length << endl;
    cout << "\t\tTheta: " << myEvent.path2.theta << endl;
    cout << "\t\tPhi: " << myEvent.path2.phi << endl;
    cout << "\t\tEnd Point: ";
    Point path2End = myEvent.path2.vertex + myEvent.path2.path2vec();
    path2End.printPt();
    cout << "\n";

    // Generate the wire layer based on the dimensions given earlier
    wireLayer myLayer = wireLayer(200., 200., 0., 3.);

    // partition the events appropriately
    int numSlices = 1;
    if (!singleSlice) {
        numSlices = 10;
    }
    vector <vector <Path>> eventSlices = myEvent.partitionEvent(numSlices, 200);

    // if you only want the one path, then just remove all of the 2nd path's contributions
    // if (singlePath) {
        // for (int i = 0; i < numSlices; i++ ) {
        //     if (eventSlices[i].size() == 2) {
        //         eventSlices[i].erase(eventSlices[i].begin() + 1);
        //     }
        // }
    // }


    // Pull the info from the eventSlices vector to reconstruct the two paths separately
    // vector <Path> path1_parts;
    // vector <Path> path2_parts;
    // for (int i = 0; i < numSlices; i++) {
    //     if (eventSlices[i].size() == 2) {
    //         path1_parts.push_back(eventSlices[i][0]);
    //         path2_parts.push_back(eventSlices[i][1]);
    //     }
    //     else if (eventSlices[i].size() == 1) {
    //         path1_parts.push_back(eventSlices[i][0]);
    //     }
    // }

    // Use this information to save the information about path partitions into a file
    ofstream pathFile;
    pathFile.open("toy_path.txt");
    Path testBound = Path();
    for (int i = 0; i < (int)eventSlices[0].size(); i++) {
        if (eventSlices[0][i] != testBound) {
            pathFile << eventSlices[0][i].vertex.x << "\t" << eventSlices[0][i].vertex.y << "\n";
            eventSlices[0][i].vertex.printPt();
            cout << " --> ";
        }
    }
    Point endPt1 = myEvent.path1.vertex + myEvent.path1.path2vec();
    endPt1.printPt();
    cout << endl;
    pathFile << endPt1.x << "\t" << endPt1.y << "\n\n";

    if (!singlePath) {
        for (int i = 0; i < (int)eventSlices[1].size(); i++) {
            if (eventSlices[1][i] != testBound) {
                pathFile << eventSlices[1][i].vertex.x << "\t" << eventSlices[1][i].vertex.y << "\n";
                eventSlices[1][i].vertex.printPt();
                cout << " --> ";
            }
        }
        Point endPt2 = myEvent.path2.vertex + myEvent.path2.path2vec();
        endPt2.printPt();
        cout << endl;
        pathFile << endPt2.x << "\t" << endPt2.y << "\n";
    }

    // Now let us build the geometry matrices and do path projection
    vector <vector <float>> allSignals(numSlices);
    vector <vector <vector <int>>> allGeoMats(numSlices);
    vector <vector <float>> allActual(numSlices);
    vector <wireLayer> allLayers(numSlices);

    ofstream wiresFile;
    wiresFile.open("toy_wires.txt");

    vector <vector <Path>> zippedSlices = zip1D(eventSlices[0], eventSlices[1]);

    if (singlePath) {
        for (int j = 0; j < (int)zippedSlices.size(); j++) {
            zippedSlices[j].erase(zippedSlices[j].begin() + 1);
        }
    }


    for (int i = 0; i < numSlices; i++) {

        if (zippedSlices[i].size() == 0) {
            continue;
        }


        vector <float> signals = myLayer.signalVec(zippedSlices[i]);

        allSignals[i] = signals;

        int vDim = myLayer.vertical.startPoints.size();
        int lDim = myLayer.lSlant.startPoints.size();
        int rDim = myLayer.rSlant.startPoints.size();

        for (int j = 0; j < vDim - 1; j++) {
            if (signals[j] > 0) {
                wiresFile << myLayer.vertical.startPoints[j].x << " " << myLayer.vertical.endPoints[j].x << "\n";
                wiresFile << myLayer.vertical.startPoints[j].y << " " << myLayer.vertical.endPoints[j].y << "\n";
                wiresFile << myLayer.vertical.startPoints[j + 1].x << " " << myLayer.vertical.endPoints[j + 1].x << "\n";
                wiresFile << myLayer.vertical.startPoints[j + 1].y << " " << myLayer.vertical.endPoints[j + 1].y << "\n";
            }
        }
        for (int j = 0; j < lDim - 1; j++) {
            if (signals[j + vDim - 1] > 0) {
                wiresFile << myLayer.lSlant.startPoints[j].x << " " << myLayer.lSlant.endPoints[j].x << "\n";
                wiresFile << myLayer.lSlant.startPoints[j].y << " " << myLayer.lSlant.endPoints[j].y << "\n";
                wiresFile << myLayer.lSlant.startPoints[j + 1].x << " " << myLayer.lSlant.endPoints[j + 1].x << "\n";
                wiresFile << myLayer.lSlant.startPoints[j + 1].y << " " << myLayer.lSlant.endPoints[j + 1].y << "\n";
            }
        }
        for (int j = 0; j < rDim - 1; j++) {
            if (signals[j + vDim + lDim - 2] > 0) {
                wiresFile << myLayer.rSlant.startPoints[j].x << " " << myLayer.rSlant.endPoints[j].x << "\n";
                wiresFile << myLayer.rSlant.startPoints[j].y << " " << myLayer.rSlant.endPoints[j].y << "\n";
                wiresFile << myLayer.rSlant.startPoints[j + 1].x << " " << myLayer.rSlant.endPoints[j + 1].x << "\n";
                wiresFile << myLayer.rSlant.startPoints[j + 1].y << " " << myLayer.rSlant.endPoints[j + 1].y << "\n";
            }
        }

        wireLayer tempLayer = myLayer;
        vector <vector <int>> geoMat = tempLayer.geoMatrix(signals);
        allGeoMats[i] = geoMat;
        vector <float> actualCharge = path2true(zippedSlices[i], tempLayer);
        allActual[i] = actualCharge;
        allLayers[i] = tempLayer;
    }
    wiresFile.close();

    // remove unambiguous and reduce dimensionality


    pair <vector <vector <vector <int>>>, vector <vector <float>>> allReduced = removeUnambigMulti(allSignals, allGeoMats, &allLayers);

    vector <vector <vector <int>>> allRedGeoMats = get<0>(allReduced);
    vector <vector <float>> allRedSigs = get<1>(allReduced);

    for (int i = 0; i < (int)allSignals.size(); i++) {
        cout << i << "\n";
        cout << "\tSignals: " << allSignals[i].size() << "\tGeoMats: " << allGeoMats[i].size() << endl;
        cout << "\tCells: " << allLayers[i].cells.size() << "\tUnambig: " << allLayers[i].unambig.size() << endl;
    }




    cout << "Running solve charge multi" << endl;
    vector <vector <float>> solvedCharges_monte = solveChargeMulti(allRedSigs, allRedGeoMats, allLayers);
    cout << "Running solve charge genetic" << endl;
    vector <vector <float>> solvedCharges_genetic = solveChargeMulti_genetic(allRedSigs, allRedGeoMats, allLayers);
    // cout << "Running solve true multi" << endl;
    // vector <vector <float>> solveTrue_matrix = solveTrue_multi(allRedSigs, allRedGeoMats);


    ofstream cellFile;
    cellFile.open("toy_cells.txt");

    for (int i = 0; i < (int)allLayers.size(); i++) {
        // the ambiguous cells
        for (int j = 0; j < (int)allLayers[i].cells.size(); j++) {
            for (int k = 0; k < (int)allLayers[i].cells[j].vertices.size(); k++) {
                cellFile << allLayers[i].cells[j].vertices[k].x << " " << allLayers[i].cells[j].vertices[k].y << " ";
            }
            cellFile << "\n" << allActual[i][j] << " " << 0 << " " << solvedCharges_genetic[i][j] << " " << solvedCharges_monte[i][j] << "\n";
            // cellFile << "\n" << allActual[i][j] << " " << solveTrue_matrix[i][j] << " " << solvedCharges_genetic[i][j] << " " << solvedCharges_monte[i][j] << "\n";
        }
        for (auto j = allLayers[i].unambig.begin(); j != allLayers[i].unambig.end(); j++) {
            for (int k = 0; k < (int)j->first.vertices.size(); k++) {
                cellFile << j->first.vertices[k].x << " " << j->first.vertices[k].y << " ";
            }
            cellFile << "\n" << j->second << " " << j->second << " " << j->second << " " << j->second << "\n";
        }

    }
    cellFile.close();


    for (int i = 0; i < (int)allActual.size(); i++) {
        float trueChargeTot = 0.;
        for (int j = 0; j < (int)allActual[i].size(); j++) {
            trueChargeTot += allActual[i][j];
        }
        cout << "True charge in layer " << i << ": " << trueChargeTot << endl;
    }

    for (int i = 0; i < (int)solvedCharges_monte.size(); i++) {
        float monteChargeTot = 0.;
        for (int j = 0; j < (int)solvedCharges_monte[i].size(); j++) {
            monteChargeTot += solvedCharges_monte[i][j];
        }
        cout << "Monte charge in layer " << i << ": " << monteChargeTot << endl;
    }



    vector <vector <float>> allRemGenetic = charge2remain(solvedCharges_genetic, allRedSigs, allRedGeoMats, allLayers);
    vector <vector <vector <float>>> geneticBestTrial = zip2D(solvedCharges_genetic, allRemGenetic);
    float costGenetic = computeCostMulti(geneticBestTrial, allRedSigs, allRedGeoMats, allLayers);

    vector <vector <float>> allRemMonte = charge2remain(solvedCharges_monte, allRedSigs, allRedGeoMats, allLayers);
    vector <vector <vector <float>>> monteBestTrial = zip2D(solvedCharges_monte, allRemMonte);
    float costMonte = computeCostMulti(monteBestTrial, allRedSigs, allRedGeoMats, allLayers);

    // vector <vector <float>> allRemMatrix = charge2remain(solveTrue_matrix, allRedSigs, allRedGeoMats, allLayers);
    // vector <vector <vector <float>>> matrixTrial = zip2D(solveTrue_matrix, allRemMatrix);
    // float costMatrix = computeCostMulti(matrixTrial, allRedSigs, allRedGeoMats, allLayers);

    vector <vector <float>> allRemActual = charge2remain(allActual, allRedSigs, allRedGeoMats, allLayers);
    vector <vector <vector <float>>> actualTrial = zip2D(allActual, allRemActual);
    float costActual = computeCostMulti(actualTrial, allRedSigs, allRedGeoMats, allLayers);


    // add in stuff for scoring later

    ofstream scoreFile;
    scoreFile.open("cellTest_scores.txt");
    // // scoreFile << sumSqrDiff(allActual, allActual) << " ";
    // // scoreFile << sumSqrDiff(solveTrue_matrix, allActual) << " ";
    // // scoreFile << sumSqrDiff(solvedCharges_genetic, allActual) << " ";
    // // scoreFile << sumSqrDiff(solvedCharges_monte, allActual) << "\n";
    // scoreFile << costActual << " " << costMatrix << " " << costGenetic << " " << costMonte << "\n";
    scoreFile << costActual << " " << 0 << " " << costGenetic << " " << costMonte << "\n";
    scoreFile.close();



    cout << "True Solution Score: " << costActual;

    int count = 0;
    for (int i = 0; i < 100; i++) {
        vector <vector <vector <float>>> mutatedTrial = mutateChargeMulti(actualTrial, allRedGeoMats);
        float mutatedCost = computeCostMulti(mutatedTrial, allRedSigs, allRedGeoMats, allLayers);
        if (mutatedCost > costActual) {
            count++;
        }
    }

    cout << "\nPercentage of time mutated is worse than true: " << count << endl;


    return 0;
}
