#include "cell.hh"
#include <queue>

using namespace std;

class Blob{
public:
    vector <Cell> cells;
    int vLims [2];
    int lLims [2];
    int rLims [2];
    void addCell(Cell c) {
        this->cells.push_back(c);
        vLims[0] = min(vLims[0], c.wireNums[0]);
        vLims[1] = max(vLims[1], c.wireNums[0]+1);
        lLims[0] = min(lLims[0], c.wireNums[1]);
        lLims[1] = max(lLims[1], c.wireNums[1]+1);
        rLims[0] = min(rLims[0], c.wireNums[2]);
        rLims[1] = max(rLims[1], c.wireNums[2]+1);
    }
    Blob() {};
    //Cell blobCell;
};

vector <Blob> tile2blob(vector <Cell> c) {
    int numVisited = 0;
    int numCells = (int)c.size();
    int visited[numCells];
    // make a quick map with int vecs to cells
    map <string, int> myMap;
    for (int i = 0; i < numCells; i++) {
        string id = to_string(c[i].wireNums[0]) + "_" + to_string(c[i].wireNums[1]) + "_" + to_string(c[i].wireNums[2]);
        myMap[id] = i;
        visited[i] = 0;
    }

    vector <Blob> allBlobs;
    while (numVisited != numCells) {
        Blob b;
        int seedRef = 0;
        for (int i = 0; i < numCells; i++) {
            if (!visited[i]) {
                seedRef = i;
                break;
            }
        }
        visited[seedRef] = 1;
        numVisited++;
        queue<int> myQ;
        myQ.push(seedRef);
        while(!myQ.empty()) {
            int f = myQ.front();
            b.addCell(c[f]);
            myQ.pop();
            int v = c[f].wireNums[0];
            int l = c[f].wireNums[1];
            int r = c[f].wireNums[2];
            map<string, int>::iterator it;
            // check the adjacent tiles and put them on the stack if they are in
            // the tiling set of cells

            // search for v+1
            string vPlus = to_string(v+1) + "_" + to_string(l) + "_" + to_string(r);
            it = myMap.find(vPlus);
            if (it != myMap.end()) {
                if (!visited[it->second]) {
                    myQ.push(it->second);
                    visited[it->second] = 1;
                    numVisited++;
                }
            }

            // search for v-1
            string vMin = to_string(v-1) + "_" + to_string(l) + "_" + to_string(r);
            it = myMap.find(vMin);
            if (it != myMap.end()) {
                if (!visited[it->second]) {
                    myQ.push(it->second);
                    visited[it->second] = 1;
                    numVisited++;
                }
            }

            // search for l+1
            string lPlus = to_string(v) + "_" + to_string(l+1) + "_" + to_string(r);
            it = myMap.find(lPlus);
            if (it != myMap.end()) {
                if (!visited[it->second]) {
                    myQ.push(it->second);
                    visited[it->second] = 1;
                    numVisited++;
                }
            }

            // search for l-1
            string lMin = to_string(v) + "_" + to_string(l-1) + "_" + to_string(r);
            it = myMap.find(lMin);
            if (it != myMap.end()) {
                if (!visited[it->second]) {
                    myQ.push(it->second);
                    visited[it->second] = 1;
                    numVisited++;
                }
            }

            // search for r+1
            string rPlus = to_string(v) + "_" + to_string(l) + "_" + to_string(r+1);
            it = myMap.find(rPlus);
            if (it != myMap.end()) {
                if (!visited[it->second]) {
                    myQ.push(it->second);
                    visited[it->second] = 1;
                    numVisited++;
                }
            }

            // search for r-1
            string rMin = to_string(v) + "_" + to_string(l) + "_" + to_string(r-1);
            it = myMap.find(rMin);
            if (it != myMap.end()) {
                if (!visited[it->second]) {
                    myQ.push(it->second);
                    visited[it->second] = 1;
                    numVisited++;
                }
            }
        }
        allBlobs.push_back(b);
    }
    return allBlobs;
}
