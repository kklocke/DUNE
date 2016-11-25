#include <limits>

// Arguments:
//     Wire layer
//     True charge on wires
// Returns:
//     Best placement of charges

// Method:
//     Identify allowed points and their associated wires
//     Iteratively:
//         Generate a random charge distribution
//         Compute potential or cost or whatever
//         Update best

// Some other ideas: try using genetic algorithms for mutating good solutions


vector <vector <float>> randChargeDistrib(vector <float> mySig, vector <vector <int>> geoMat) {
    vector <float> remain = mySig;
    vector <float> myCharge = (geoMat[0].size(), 0);
    for (int i = 0; i < (int)myCharge.size(); i++) {
        float minRem = numeric_limits<float>::infinity();
        vector <int> toDec;
        for (int j = 0; j < (int)mySig.size(); j++) {
            if (geoMat[j][i] > 0.) {
                minRem = min(remain[j], minRem);
                toDec.push_back(j);
            }
        }
        float charge = float(rand() * minRem) / float(RAND_MAX);
        myCharge[i] = charge;
        for (int j = 0; j < (int)toDec.size(); j++) {
            remain[toDec[j]] -= charge;
        }
    }
    vector <vector <float>> toReturn;
    toReturn.push_back(myCharge);
    toReturn.push_back(remain);
    return toReturn;
}

vector <float> solveCharge(vector <float> mySig, vector <vector <int>> geoMat) {
    int tryCount = 0;
    float minCost = numeric_limits<float>::infinity();
    vector <float> bestDistrib (geoMat[0].size(), -1);
    while (tryCount < 10000):
        vector <vector <float>> myTrial = randChargeDistrib(mySig, geoMat);
        float trialCost = computeCost(myTrial, mySig, geoMat);
        if (trialCost < minCost):
            minCost = trialCost;
            bestDistrib = myTrial[0];
            tryCount = 0;
        }
        tryCount++;
    }
    return bestDistrib;
}

float computeCost(vector <vector <float>> trial, vector <float> mySig,vector<vector <int>> geoMat) {
    float myCost = 0;
    for (int i = 0; i < (int)mySig.size(); i++) {
        cost += (trial[1][i]^2)/mySig[i];
    }

    

    // add stuff about spacial packing

    return myCost;
}




// MIGHT BE USEFUL TO TRY A DISCRETE VERSION OF THIS WHERE I JUST SELECT
// A RANDOM GROUP OF THE POINTS AND GIVE THEM THE UNIT CHARGE. BUT I WILL
// TRY BOTH WAYS AN SEE WHCH ONE GIVES BETTER RESULTS FOR THE TIME BEING
// THE ONE WITH FLOATS WILL PROBABLY BETTER MINIMIZE THE COST BUT THE
// RESULTS FROM A DISCRETE CHARGE PACKET VERSION WILL BE MORE PHYSICALLY MEANINGFUL
