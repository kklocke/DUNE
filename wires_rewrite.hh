#include "events.hh"
#include <algorithm>

using namespace std;

class wireArray
{
public:
    float pitch; // mm
    float angle; // angle wrt x axis
    float height; // y dim
    float length; // x dim
    float z; // z posn
    vector <Point> startPoints; // points defining start of wire boundaries
    vector <Point> endPoints; // points defining end of wire boundaries
    void operator=(const wireArray &w) {
		this->pitch = w.pitch;
		this->angle = w.angle;
		this->height = w.height;
		this->length = w.length;
		this->z = w.z;
        this->endPoints.clear();
        this->startPoints.clear();
		for (int i = 0; i < (int)w.startPoints.size(); i++) {
			this->startPoints.push_back(w.startPoints[i]);
		}
		for (int i = 0; i < (int)w.endPoints.size(); i++) {
			this->endPoints.push_back(w.endPoints[i]);
		}
	}
    wireArray(float pitch, float angle, float l, float h, float z) : pitch(pitch), angle(angle), height(h), length(l), z(z) {
		float xPos = 0., yPos = h;
		float xIncr = pitch / fabs(sin(M_PI * angle / 180));
		float yIncr = pitch / fabs(cos(M_PI * angle / 180));
		double slope = tan(M_PI * angle / 180.);
		if (int(angle) == 90)
		{
			while (xPos < l)
			{
				startPoints.push_back(Point(xPos, h, z));
				endPoints.push_back(Point(xPos, 0, z));
				xPos += xIncr;
			}
		}
		else if (angle == 0.)
		{
			while (yPos > 0)
			{
				startPoints.push_back(Point(l, yPos, z));
				endPoints.push_back(Point(0, yPos, z));
				yPos -= yIncr;
			}
		}
		else if (slope > 0)
		{
			while (xPos < l)
			{
				startPoints.push_back(Point(xPos, yPos, z));
				float yInt = yPos - xPos * slope;
				if ((yInt >= 0) && (yInt <= h))
				{
					endPoints.push_back(Point(0., yInt, z));
				}
				else
				{
					float xInt = xPos - (yPos / slope);
					endPoints.push_back(Point(xInt, 0., z));
				}
				xPos += xIncr;
			}

			// solve back for the y posn
			yPos = h - slope*(xPos - l);
			xPos = l;
			while (yPos > 0)
			{
				startPoints.push_back(Point(xPos, yPos, z));
				float xInt = l - yPos / slope;
				if ((xInt > 0) && (xInt < l))
				{
					endPoints.push_back(Point(xInt, 0., z));
				}
				else
				{
					float yInt = yPos - l*slope;
					endPoints.push_back(Point(0., yInt, z));
				}
				yPos -= yIncr;
			}
		}
		else if (slope < 0)
		{
			while (xPos < l)
			{
				startPoints.push_back(Point(xPos, yPos, z));
				float xInt = xPos - h / slope;
				if ((xInt > 0) && (xInt < h))
				{
					endPoints.push_back(Point(xInt, 0., z));
				}
				else
				{
					float yInt = yPos  + slope * (h - xPos);
					endPoints.push_back(Point(l, yInt, z));
				}
				xPos += xIncr;
			}
			xPos = 0.;
			yPos = h - yIncr;
			while (yPos > 0)
			{
				startPoints.insert(startPoints.begin(), Point(0, yPos, z));
				float xInt = xPos - yPos / slope;
				if ((xInt > 0) && (xInt < l))
				{
					endPoints.insert(endPoints.begin(), Point(xInt, 0., z));
				}
				else
				{
					float yInt = yPos - slope * h;
					endPoints.insert(endPoints.begin(), Point(l, yInt, z));
				}
				yPos -= yIncr;
			}
		}
	}
	wireArray() : pitch(0.), angle(0.), height(0.), length(0.), z(0.) {};
};


// I haven't changed anything in wire array. I think it should be fine. What we care about is the intersecting stuff and solving for
// true charges.


vector <int> crossings(wireArray myArr, Point p1, Point p2)
{
    // figure out which lines are crossed in the wire array
	if (p1.x > p2.x)
	{
		Point p3 = p1;
		p1 = p2;
		p2 = p3;
	}
    vector <int> myCross;
	int firstWire = -1;
	int lastWire = -1;
	// cout << "my points to cross: \n\t";
	// p1.printPt();
	// cout << "\n\t";
	// p2.printPt();
	// cout << endl;
	for (int i = 0; i < (int)myArr.endPoints.size(); i++)
	{
		Point refPt1 = myArr.startPoints[i];
		Point refPt2 = myArr.endPoints[i];
		float slope = (refPt2.y - refPt1.y) / (refPt2.x - refPt1.x);
		// y at given x
		float refY_1 = refPt2.y + (slope * (p1.x - refPt2.x));
		float refY_2 = refPt2.y + (slope * (p2.x - refPt2.x));
		if ((slope > 0) && (p1.y >= refY_1) && (firstWire == -1))
		{
			firstWire = i;
		}
		else if ((slope < 0) && (p1.y <= refY_1) && (firstWire == -1))
		{
			firstWire = i;
		}
		if ((slope > 0) && (p2.y > refY_2) && (lastWire == -1))
		{
			lastWire = i;
		}
		else if ((slope < 0) && (p2.y < refY_2) && (lastWire == -1))
		{
			lastWire = i;
		}
		if ((firstWire != -1) && (lastWire != -1))
		{
			break;
		}
	}
	if (lastWire == -1) {
		lastWire = (int)myArr.endPoints.size();
	}
	// cout << "FIRST WIRE: " << firstWire << "\tLAST WIRE: " << lastWire << endl;
	for (int i = min(firstWire, lastWire); i < max(lastWire, firstWire); i++)
	{
		myCross.push_back(i);
		// cout << i << "\t";
	}
	// cout << endl;
    return myCross;
}


// I need to check that the indices are in the right order
vector <int> numCellPerWire(vector <float> mySig, vector <vector <int>> geoMat) {
    // cout << "Numcell per wire: " << mySig.size() << "\t" << geoMat.size() << "\t" << geoMat[0].size() << endl;
    if (mySig.empty()) {
        vector <int> toRet;
        return toRet;
    }
    assert (mySig.size() == geoMat[0].size());
    vector <int> numCells ((int)geoMat.size(), 0);
    for (int i = 0; i < (int)geoMat.size(); i++) {
        for (int j = 0; j < (int)geoMat[i].size(); j++) {
            if ((mySig[j] > 0) && (geoMat[i][j] == 1)) {
                numCells[i] += 1;
            }
        }
    }
    return numCells;
}

class wireLayer
{
public:
    wireArray vertical;
    wireArray lSlant;
    wireArray rSlant;
    float pitch;
    float length;
    float height;
    float z;
    vector <Point> grids;
    vector <Cell> cells; // Now we care about cells not grids, since the grids are meaningless
    vector <Cell> ambigCells;
    map<Cell, float> unambig;
    wireLayer(float h, float l, float z, float pitch) : pitch(pitch), length(l), height(h), z(z) {
        vertical = wireArray(pitch, 90., l, h, z);
        lSlant = wireArray(pitch, 54.3, l, h, z);
        rSlant = wireArray(pitch, -54.3, l, h, z);
    }
    wireLayer() {
        pitch = 0.;
        length = 0.;
        height = 0.;
        z = 0.;
    }
    void operator=(const wireLayer &w) {
        this->pitch = w.pitch;
		this->length = w.length;
		this->height = w.height;
		this->z = w.z;
		this->vertical = w.vertical;
		this->lSlant = w.lSlant;
		this->rSlant = w.rSlant;
        this->cells = w.cells;
        this->ambigCells = w.ambigCells;
        this->unambig = w.unambig;
    }

	vector <vector <int>> allCrossings(Point p1, Point p2) {
		vector <vector <int>> all3;
		all3.push_back(crossings(vertical, p1, p2));
		all3.push_back(crossings(lSlant, p1, p2));
		all3.push_back(crossings(rSlant, p1, p2));
		return all3;
	}

    bool inFrame(Point p) {
        if ((p.x < 0) || (p.x > this->length)) {return false;}
        if ((p.y < 0) || (p.y > this->length)) {return false;}
        return true;
    }

    /* geoMatrix: constructs the (reduced dimension) geometry matrix
     * Rows: correspond to the wires, ordered vertical, lSlant, rSlant
     * Columns: the relevant cells (based on the hit data)
     * Input: vector of floats indicating hits, this data is required to reduced
     * the dimensionality of the matrix so as to save on computation time.
     * Returns: the 2-D geometry matrix
     */

     // THIS IS GOING TO NEED SOME DEBUGGING TO CHECK THAT IT WORKS PROPERLY

     // I should change the hit thing so that it maps floats to cell indicators (i.e. 3 ints), this will make everything faster

    vector <vector <int>> geoMatrix(vector <float> hit) {
        int vSize = vertical.startPoints.size();
        int lSize = lSlant.startPoints.size();
        int rSize = rSlant.startPoints.size();
        int numRows = vSize + lSize + rSize - 3;
        vector <vector <int>> myMatrix;
        vector <vector <int>> transposeMat;
        cout << "vSize: " << vSize << endl;
        cout << "lSize: " << lSize << endl;
        cout << "rSize: " << rSize << endl;
        cells.clear();
        unambig.clear();
        ambigCells.clear();

        // generate all of the cells, then construct matrix accordingly
        // cout << "Entering the loop\n";
        for (int i = 0; i < vSize - 1; i++) {
        	if (hit[i] == 0) {
        		continue;
        	}
        	vector <vector <Point>> v;
        	vector <Point> tempV;
        	tempV.push_back(vertical.startPoints[i]);
        	tempV.push_back(vertical.endPoints[i]);
        	v.push_back(tempV);
        	tempV.clear();
        	tempV.push_back(vertical.startPoints[i+1]);
        	tempV.push_back(vertical.endPoints[i+1]);
        	v.push_back(tempV);
        	tempV.clear();
        	for (int j = 0; j < lSize - 1; j++) {
        		if (hit[vSize + j - 1] == 0) {
        			continue;
        		}
        		vector <vector <Point>> l;
        		vector <Point> tempL;
        		tempL.push_back(lSlant.startPoints[j]);
        		tempL.push_back(lSlant.endPoints[j]);
        		l.push_back(tempL);
        		tempL.clear();
        		tempL.push_back(lSlant.startPoints[j+1]);
        		tempL.push_back(lSlant.endPoints[j+1]);
        		l.push_back(tempL);
        		tempL.clear();
        		for (int k = 0; k < rSize - 1; k++) {
        			if (hit[vSize + lSize + k - 2] == 0) {
        				continue;
        			}
        			vector <vector <Point>> r;
        			vector <Point> tempR;
        			tempR.push_back(rSlant.startPoints[k]);
        			tempR.push_back(rSlant.endPoints[k]);
        			r.push_back(tempR);
        			tempR.clear();
        			tempR.push_back(rSlant.startPoints[k+1]);
        			tempR.push_back(rSlant.endPoints[k+1]);
        			r.push_back(tempR);
        			tempR.clear();
        			int myNums [3] = {i, j, k};
                    // cout << "Trying wire triplet: " << i << "\t" << j << "\t" << k << endl;
                    // cout << "\tLimits (v, l, r): " << vSize << ", " << lSize << ", " << rSize << endl;
        			Cell myCell = Cell(v, l, r, myNums);
                    // cout << "\t Made the cell\n";
                    if (myCell.vertices.size() >= 3) {
        				int flag = 0;
        				for (int c = 0; c < (int)myCell.vertices.size(); c++) {
        					Point myCheck = myCell.vertices[c];
                            if (!inFrame(myCheck)) {
                                flag = 1;
                                break;
                            }
        				}
        				if (flag == 0) {
        					this->cells.push_back(myCell);
                            // SET THE Transposed matrix Rows
                            // cout << "trying to make a matrix row\t";
                            vector <int> tempRow (numRows, 0);
                            tempRow[i] = 1;
                            tempRow[j + vSize - 1] = 1;
                            tempRow[k + vSize + lSize - 2] = 1;
                            transposeMat.push_back(tempRow);
                            // cout << "made the row \n";
        				}
        			}
        		}
        	}
        }
        // cout << "Out of the loop, trying to transpose" << endl;
        // Now transpose the matrix to get myMatrix
        int numCols = transposeMat.size();
        for (int i = 0; i < numRows; i++) {
            vector <int> tempRow(numCols, 0);
            for (int j = 0; j < numCols; j++) {
                tempRow[j] = transposeMat[j][i];
            }
            myMatrix.push_back(tempRow);
        }
        return myMatrix;
        // return transposeMat;
    }

    // pair<vector <vector <int>>, vector <float>> removeUnambig(vector <float> signals, vector <vector <int>> geoMat) {
    //     vector <float> fake(geoMat[0].size(), 1.);
    //     vector <int> numCells = numCellPerWire(fake, geoMat);
    //     vector <float> newSig = signals;
    //     vector <vector <int>> newGeoMat = geoMat;
    //     this->ambigCells = this->cells;
    //
    //     // cells to remove
    //     vector <int> cellIndices;
    //
    //     for (int i = 0; i < (int)geoMat.size(); i++) {
    //         if (numCells[i] == 1) {
    //             int myCellIndex = -1;
    //             for (int j = 0; j < (int)geoMat[i].size(); j++) {
    //                 if (geoMat[i][j] == 1) {
    //                     myCellIndex = j;
    //                     break;
    //                 }
    //             }
    //             if (myCellIndex == -1) {
    //                 continue;
    //             }
    //             for (int k = 0; k < (int)geoMat.size(); k++) {
    //                 if (geoMat[k][myCellIndex] == 1) {
    //                     newSig[k] -= signals[i];
    //                 }
    //             }
    //             cellIndices.push_back(myCellIndex);
    //             this->unambig[this->cells[myCellIndex]] = signals[i];
    //         }
    //     }
    //     cout << "Cell indices size: " << cellIndices.size() << endl;
    //     sort(cellIndices.begin(), cellIndices.end());
    //     for (int i = (int) newGeoMat.size() - 1; i >= 0; i--) {
    //         if (newSig[i] == 0) {
    //             newGeoMat.erase(newGeoMat.begin() + i);
    //             newSig.erase(newSig.begin() + i);
    //         }
    //         // else {
    //         //     for (int j = (int)cellIndices.size() - 1; j >= 0; j--) {
    //         //         newGeoMat[i].erase(newGeoMat[i].begin() + cellIndices[j]);
    //         //     }
    //         // }
    //     }
    //
    //     for (int i = (int)cellIndices.size() - 1; i >= 0; i--) {
    //         this->ambigCells.erase(this->ambigCells.begin() + cellIndices[i]);
    //         for (int j = 0; j < (int)newGeoMat.size(); j++) {
    //             vector <int> temp = newGeoMat[j];
    //             temp.erase(temp.begin() + cellIndices[i]);
    //             newGeoMat[j] = temp;
    //         }
    //     }
    //
    //     pair<vector <vector <int>>, vector <float>> toRet = make_pair(newGeoMat, newSig);
    //     return toRet;
    // }


    // pair <vector <vector <int>>, vector <float>> removeUnambig(vector <float> signals, vector <vector <int>> geoMat) {
    //
    //     vector <vector <int>> newGeoMat = geoMat;
    //     vector <float> newSig = signals;
    //     vector <float> fake (geoMat[0].size(), 1.);
    //     vector <int> numCells = numCellPerWire(fake, geoMat);
    //
    //     for (int i = 0; i < (int)numCells.size(); i++) {
    //         cout << numCells[i] << " ";
    //     }
    //     cout << endl;
    //
    //     int i = 0;
    //
    //     while (i != (int)newGeoMat.size()) {
    //         // find a wire that has unambiguous cell
    //         if (numCells[i] == 1) {
    //             cout << i << " ";
    //             float myCharge = newSig[i];
    //             int myCellIndex = -1;
    //             // identify which cell
    //             for (int j = 0; j < (int)newGeoMat[i].size(); j++) {
    //                 if (newGeoMat[i][j] == 1) {
    //                     myCellIndex = j;
    //                     break;
    //                 }
    //                 // for all the wires it touches, reduce the charge, and then delete that cell from the matrix
    //                 if (myCellIndex != -1) {
    //                     // write the cell into the unambig map and then remove from the regular cells vector
    //                     this->unambig[this->cells[myCellIndex]] = myCharge;
    //                     this->cells.erase(this->cells.begin() + myCellIndex);
    //                     for (int k = (int)newGeoMat.size() - 1; k >= 0; k--) {
    //                         if (newGeoMat[k][myCellIndex] == 1) {
    //                             // reduce charge
    //                             newSig[k] -= myCharge;
    //                         }
    //                         // remove the cell
    //                         vector <int> temp = newGeoMat[k];
    //                         temp.erase(temp.begin() + myCellIndex);
    //                         newGeoMat[k] = temp;
    //                         // check if the wire is now empty of charge
    //                         if (newSig[k] == 0.) {
    //                             // if so, remove the wire from geomat and signal
    //                             newGeoMat.erase(newGeoMat.begin() + k);
    //                             newSig.erase(newSig.begin() + k);
    //                         }
    //                     }
    //                 }
    //             }
    //             cout << endl;
    //         }
    //
    //         i++;
    //     }
    //     pair <vector <vector <int>>, vector <float>> toRet = make_pair(newGeoMat, newSig);
    //     return toRet;
    // }

    bool helpRemoveUnambig(vector <float> *sig, vector <vector <int>> *geoMat, vector <int> *numCells) {
        int i = 0;
	    bool toRet = false;
	    while ((i < (int)(*sig).size()) && (i < (int)(*numCells).size())) {
	       if (((*numCells)[i] == 0) || ((*sig)[i] == 0)) {
               toRet = true;
               (*geoMat).erase((*geoMat).begin() + i);
               (*sig).erase((*sig).begin() + i);
               (*numCells).erase((*numCells).begin() + i);
               // remove this wire entirely
               break;
           }
           if ((*numCells)[i] == 1) {
               toRet = true;
               // reduce the charge on the relevant wires and remove this particular wire
               int cellIndex = -1;
               for (int j = 0; j < (int)(*geoMat)[i].size(); j++) {
                   // temp
                   assert((*geoMat)[i].size() == this->cells.size());
                   if ((*geoMat)[i][j] == 1){
                       cellIndex = j;
                       break;
                   }
               }
               if (cellIndex == -1) {
                   continue;
               }
               float toDec = (*sig)[i];
               for (int j = 0; j < (int)(*geoMat).size(); j++) {
                   if ((*geoMat)[j][cellIndex] == 1) {
                       (*sig)[j] -= toDec;
                       (*numCells)[j]--;
                   }
                   // remove that column
                   vector <int> temp = (*geoMat)[j];
                   temp.erase(temp.begin() + cellIndex);
                   (*geoMat)[j] = temp;
                   // (*geoMat)[j].erase((*geoMat)[j].begin() + cellIndex);
               }
               // remove that cell
               cout << "This->cells.size() " << this->cells.size() << "\tCell Index: " << cellIndex << endl;
               Cell c = this->cells[cellIndex];
               this->unambig[c] = toDec;
               this->cells.erase(this->cells.begin() + cellIndex);
               break;
           }
           i++;
       }
       return toRet;
   }


    pair <vector <vector <int>>, vector <float>> removeUnambig(vector <float> signals, vector <vector <int>> geoMat) {
        vector <vector <int>> newGeoMat = geoMat;
        vector <float> newSig = signals;
        // cout << "Original Size of Sigs: " << signals.size() << endl;
        // cout << "Original Size of GeoMat: " << geoMat.size() << "\t" << geoMat[0].size() << endl;
        vector <float> fake (geoMat[0].size(), 1.);
        vector <int> numCells = numCellPerWire(fake, newGeoMat);

        cout << "Num Cells Vector: ";
        for (int i = 0; i < (int)numCells.size(); i++) {
            cout << numCells[i] << " ";
        }
        cout << "\n";

        bool hasChanged = true;

        while (hasChanged) {
            hasChanged = helpRemoveUnambig(&newSig, &newGeoMat, &numCells);
        }

        // cout << "Final Size of Sigs: " << newSig.size() << endl;
        // cout << "Final Size of GeoMat: " << newGeoMat.size() << "\t" << newGeoMat[0].size() << endl;

        pair <vector <vector <int>>, vector <float>> toRet = make_pair(newGeoMat, newSig);
        return toRet;
    }




    vector <float> signalVec(Path p) {
        int vDim = this->vertical.startPoints.size();
        int lDim = this->lSlant.startPoints.size();
        int rDim = this->rSlant.startPoints.size();
        vector <float> sigVec(vDim + lDim + rDim - 3, 0.);
        // go in short intervals?
        // break it up and find which cells I go through

        map <vector<int>, int> cellSeen;
        float stepSize = this->pitch / 10.;
        int numSteps = abs(int(p.length / stepSize));
        Point currentPt = p.vertex;
        Point step = p.path2vec().scalarMult(stepSize / p.length);
        for (int i = 0; i <= numSteps; i++) {
            if (i == numSteps) {
                currentPt = p.vertex + p.path2vec();
            }
            // find vertical boundaries
            int vPlace = -1;
            int lPlace = -1;
            int rPlace = -1;
            for (int j = 0; j < vDim; j++) {
                if (this->vertical.startPoints[j].x > currentPt.x) {
                    vPlace = j - 1;
                    break;
                }
            }
            for (int j = 0; j < lDim; j++) {
                float slope = tan(this->lSlant.angle * M_PI / 180.);
                float predY = this->lSlant.startPoints[j].y + slope * (currentPt.x - this->lSlant.startPoints[j].x);
                if (predY < currentPt.y) {
                    lPlace = j - 1;
                    break;
                }
            }

            for (int j = 0; j < rDim; j++) {
                float slope = tan(this->rSlant.angle * M_PI / 180.);
                float predY = this->rSlant.startPoints[j].y + slope * (currentPt.x - this->rSlant.startPoints[j].x);
                if (predY > currentPt.y) {
                    rPlace = j - 1;
                    break;
                }
            }
            vector <int> myNums;
            myNums.push_back(vPlace);
            myNums.push_back(lPlace);
            myNums.push_back(rPlace);
            cellSeen[myNums] = 1;
            currentPt += step;
        }
        for (map<vector <int>, int>::iterator j = cellSeen.begin(); j != cellSeen.end(); j++) {
            //add charge on the "wires" associated with the cell,
            // "wires" indexing has the same index as the bounding lines
            vector <int> myNums = j->first;
            sigVec[myNums[0]] += 10;
            sigVec[myNums[1] + vDim - 1] += 10;
            sigVec[myNums[2] + vDim + lDim - 2] += 10;
        }

        return sigVec;
    }

    vector <float> signalVec(vector <Path> myPaths) {
        int vDim = this->vertical.startPoints.size();
        int lDim = this->lSlant.startPoints.size();
        int rDim = this->rSlant.startPoints.size();
        vector <float> sigVec(vDim + lDim + rDim - 3, 0.);
        map <vector<int>, int> cellSeen;
        float stepSize = this->pitch / 10.;

        Path defPath = Path();

        for (int k = 0; k < (int)myPaths.size(); k++) {
            Path p = myPaths[k];
            if (p == defPath) {
                continue;
            }
            int numSteps = abs(int(p.length / stepSize));
            Point currentPt = p.vertex;
            Point step = p.path2vec().scalarMult(stepSize / p.length);
            for (int i = 0; i <= numSteps; i++) {
                if (i == numSteps) {
                    currentPt = p.vertex + p.path2vec();
                }
                // find vertical boundaries
                int vPlace = -1;
                int lPlace = -1;
                int rPlace = -1;
                for (int j = 0; j < vDim; j++) {
                    if (this->vertical.startPoints[j].x > currentPt.x) {
                        vPlace = j - 1;
                        break;
                    }
                }
                for (int j = 0; j < lDim; j++) {
                    float slope = tan(this->lSlant.angle * M_PI / 180.);
                    float predY = this->lSlant.startPoints[j].y + slope * (currentPt.x - this->lSlant.startPoints[j].x);
                    if (predY < currentPt.y) {
                        lPlace = j - 1;
                        break;
                    }
                }

                for (int j = 0; j < rDim; j++) {
                    float slope = tan(this->rSlant.angle * M_PI / 180.);
                    float predY = this->rSlant.startPoints[j].y + slope * (currentPt.x - this->rSlant.startPoints[j].x);
                    if (predY > currentPt.y) {
                        rPlace = j - 1;
                        break;
                    }
                }
                vector <int> myNums;
                myNums.push_back(vPlace);
                myNums.push_back(lPlace);
                myNums.push_back(rPlace);
                if (cellSeen.find(myNums) == cellSeen.end()) {
                    cellSeen[myNums]++;
                }
                else {
                    cellSeen[myNums] = 1;
                }
                currentPt += step;
            }
        }

        for (map<vector <int>, int>::iterator j = cellSeen.begin(); j != cellSeen.end(); j++) {
            //add charge on the "wires" associated with the cell,
            // "wires" indexing has the same index as the bounding lines
            vector <int> myNums = j->first;
            int numHits = j->second;
            sigVec[myNums[0]] += 10 * numHits;
            sigVec[myNums[1] + vDim - 1] += 10  * numHits;
            sigVec[myNums[2] + vDim + lDim - 2] += 10  * numHits;
        }
        return sigVec;
    }
};


pair <vector <vector <vector <int>>>, vector <vector <float>>> removeUnambigMulti(vector <vector <float>> allSigs, vector <vector <vector <int>>> allGeoMats, vector <wireLayer> *allLayers) {
    int numTrials = (int)allSigs.size();

    vector <vector <vector <int>>> newGeoMats(numTrials);
    vector <vector <float>> newSigs(numTrials);

    for (int i = 0; i < (int)allGeoMats.size(); i++) {
        if (allSigs[i].size() == 0) {
            continue;
        }
        pair <vector <vector <int>>, vector <float>> myNew = (*allLayers)[i].removeUnambig(allSigs[i], allGeoMats[i]);
        newGeoMats[i] = get<0>(myNew);
        newSigs[i] = get<1>(myNew);
    }

    pair <vector <vector <vector <int>>>, vector <vector <float>>> toRet = make_pair(newGeoMats, newSigs);
    return toRet;
}





vector <float> path2true(Path p, wireLayer myLayer) {
    int vDim = myLayer.vertical.startPoints.size();
    int lDim = myLayer.lSlant.startPoints.size();
    int rDim = myLayer.rSlant.startPoints.size();
    vector <float> sigVec(vDim + lDim + rDim - 3, 0.);
    // go in short intervals?
    // break it up and find which cells I go through
    map <vector<int>, int> cellSeen;
    float stepSize = myLayer.pitch / 10.;
    int numSteps = int(p.length / stepSize);
    Point currentPt = p.vertex;
    Point step = p.path2vec().scalarMult(stepSize / p.length);
    for (int i = 0; i <= numSteps; i++) {
        if (i == numSteps) {
            currentPt = p.vertex + p.path2vec();
        }
        // find vertical boundaries
        int vPlace = -1;
        int lPlace = -1;
        int rPlace = -1;
        for (int j = 0; j < vDim; j++) {
            if (myLayer.vertical.startPoints[j].x > currentPt.x) {
                vPlace = j - 1;
                break;
            }
        }
        for (int j = 0; j < lDim; j++) {
            float slope = tan(myLayer.lSlant.angle * M_PI / 180.);
            float predY = myLayer.lSlant.startPoints[j].y + slope * (currentPt.x - myLayer.lSlant.startPoints[j].x);
            if (predY < currentPt.y) {
                lPlace = j - 1;
                break;
            }
        }

        for (int j = 0; j < rDim; j++) {
            float slope = tan(myLayer.rSlant.angle * M_PI / 180.);
            float predY = myLayer.rSlant.startPoints[j].y + slope * (currentPt.x - myLayer.rSlant.startPoints[j].x);
            if (predY > currentPt.y) {
                rPlace = j - 1;
                break;
            }
        }
        vector <int> myNums;
        myNums.push_back(vPlace);
        myNums.push_back(lPlace);
        myNums.push_back(rPlace);
        cellSeen[myNums] = 1;
        currentPt += step;
    }
    cout << "NUM CELLS THAT SHOULD BE ACTIVE: " << cellSeen.size() << endl;
    vector <float> myCharge((int)myLayer.cells.size(), 0.);
    for (int i = 0; i < (int)myLayer.cells.size(); i++) {
        for (map<vector <int>, int>::iterator j = cellSeen.begin(); j != cellSeen.end(); j++) {
            vector <int> myNums = j->first;
            if ((myNums[0] == myLayer.cells[i].wireNums[0]) && (myNums[1] == myLayer.cells[i].wireNums[1]) && (myNums[2] == myLayer.cells[i].wireNums[2])) {
                myCharge[i] += 10.;
            }
        }
    }
    return myCharge;
}


vector <float> path2true(vector <Path> allPaths, wireLayer myLayer) {
    int vDim = myLayer.vertical.startPoints.size();
    int lDim = myLayer.lSlant.startPoints.size();
    int rDim = myLayer.rSlant.startPoints.size();
    vector <float> sigVec(vDim + lDim + rDim - 3, 0.);
    // go in short intervals?
    // break it up and find which cells I go through
    map <vector<int>, int> cellSeen;
    float stepSize = myLayer.pitch / 10.;

    Path defPath = Path();

    for (int k = 0; k < (int)allPaths.size(); k++) {
        map <vector <int>, int> thisPathSeen;
        Path p = allPaths[k];
        if (p == defPath) {
            continue;
        }
        int numSteps = int(p.length / stepSize);
        Point currentPt = p.vertex;
        Point step = p.path2vec().scalarMult(stepSize / p.length);
        for (int i = 0; i <= numSteps; i++) {
            if (i == numSteps) {
                currentPt = p.vertex + p.path2vec();
            }
            // find vertical boundaries
            int vPlace = -1;
            int lPlace = -1;
            int rPlace = -1;
            for (int j = 0; j < vDim; j++) {
                if (myLayer.vertical.startPoints[j].x > currentPt.x) {
                    vPlace = j - 1;
                    break;
                }
            }
            for (int j = 0; j < lDim; j++) {
                float slope = tan(myLayer.lSlant.angle * M_PI / 180.);
                float predY = myLayer.lSlant.startPoints[j].y + slope * (currentPt.x - myLayer.lSlant.startPoints[j].x);
                if (predY < currentPt.y) {
                    lPlace = j - 1;
                    break;
                }
            }

            for (int j = 0; j < rDim; j++) {
                float slope = tan(myLayer.rSlant.angle * M_PI / 180.);
                float predY = myLayer.rSlant.startPoints[j].y + slope * (currentPt.x - myLayer.rSlant.startPoints[j].x);
                if (predY > currentPt.y) {
                    rPlace = j - 1;
                    break;
                }
            }
            vector <int> myNums;
            myNums.push_back(vPlace);
            myNums.push_back(lPlace);
            myNums.push_back(rPlace);
            if (cellSeen.find(myNums) == cellSeen.end()) {
                cellSeen[myNums] = 1;
                thisPathSeen[myNums] = 1;
            }
            else {
                if (thisPathSeen.find(myNums) == thisPathSeen.end()) {
                    thisPathSeen[myNums] = 1;
                    cellSeen[myNums]++;
                }
            }
            currentPt += step;
        }
    }
    vector <float> myCharge((int)myLayer.cells.size(), 0.);
    for (int i = 0; i < (int)myLayer.cells.size(); i++) {
        for (map<vector <int>, int>::iterator j = cellSeen.begin(); j != cellSeen.end(); j++) {
            vector <int> myNums = j->first;
            if ((myNums[0] == myLayer.cells[i].wireNums[0]) && (myNums[1] == myLayer.cells[i].wireNums[1]) && (myNums[2] == myLayer.cells[i].wireNums[2])) {
                myCharge[i] += 10. * j->second;
            }
        }
    }
    return myCharge;
}










vector <vector <float>> randChargeDistrib(vector <float> mySig, vector <vector <int>> geoMat, wireLayer myLayer) {
    vector <float> remain = mySig;
    if (myLayer.cells.empty() || mySig.empty()) {
        vector <float> temp1;
        vector <float> temp2;
        vector <vector <float>> toRet;
        toRet.push_back(temp1);
        toRet.push_back(temp2);
        return toRet;
    }

    vector <float> myCharge ((int)geoMat[0].size(), 0.);
    // int vDim = myLayer.vertical.startPoints.size();
    // int lDim = myLayer.lSlant.startPoints.size();
    // int rDim = myLayer.rSlant.startPoints.size();
    vector <int> myOrder((int)myLayer.cells.size(), 0);
    for (int i = 0; i < (int)myOrder.size(); i++) {
        myOrder[i] = i;
    }
    random_shuffle(myOrder.begin(), myOrder.end());
    for (int i = 0; i < (int)myLayer.cells.size(); i++) {
        vector <int> relevantWires;
        for (int j = 0; j < (int)geoMat.size(); j++) {
            if (geoMat[j][myOrder[i]] == 1) {
                relevantWires.push_back(j);
            }
        }


        float minRem = numeric_limits<float>::infinity();

        for (int j = 0; j < (int)relevantWires.size(); j++) {
            minRem = min(minRem, remain[relevantWires[j]]);
        }

        // float f1 = remain[myLayer.cells[myOrder[i]].wireNums[0]];
        // float f2 = remain[myLayer.cells[myOrder[i]].wireNums[1] + vDim - 1];
        // float f3 = remain[myLayer.cells[myOrder[i]].wireNums[2] + vDim + lDim - 2];
        // minRem = min(minRem, f1);
        // minRem = min(minRem, f2);
        // minRem = min(minRem, f3);
        // vector <int> vDec;
        // vector <int> lDec;
        // vector <int> rDec;
        // if (myLayer.cells[myOrder[i]].touches[0][0]) {
        //     minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[0]]);
        //     vDec.push_back(myLayer.cells[myOrder[i]].wireNums[0]);
        // }
        // if (myLayer.cells[myOrder[i]].touches[0][1]) {
        //     minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[0] + 1]);
        //     vDec.push_back(myLayer.cells[myOrder[i]].wireNums[0] + 1);
        // }
        // if (myLayer.cells[myOrder[i]].touches[1][0]) {
        //     minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[1] + vDim]);
        //     lDec.push_back(myLayer.cells[myOrder[i]].wireNums[1] + vDim);
        // }
        // if (myLayer.cells[myOrder[i]].touches[1][1]) {
        //     minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[1] + vDim + 1]);
        //     lDec.push_back(myLayer.cells[myOrder[i]].wireNums[1] + vDim + 1);
        // }
        // if (myLayer.cells[myOrder[i]].touches[2][0]) {
        //     minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[2] + vDim + lDim]);
        //     rDec.push_back(myLayer.cells[myOrder[i]].wireNums[2] + vDim + lDim);
        // }
        // if (myLayer.cells[myOrder[i]].touches[2][1]) {
        //     minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[2] + vDim + lDim + 1]);
        //     rDec.push_back(myLayer.cells[myOrder[i]].wireNums[2] + vDim + lDim + 1);
        // }
        // if (minRem > 0) {
        //     // put charge in the cells
        //     myCharge[myOrder[i]] = 10.;
        //     // decrement the wires
        //     for (int j = 0; j < (int)vDec.size(); j++) {
        //         remain[vDec[j]] -= 10. / float(vDec.size());
        //     }
        //     for (int j = 0; j < (int)lDec.size(); j++) {
        //         remain[lDec[j]] -= 10. / float(lDec.size());
        //     }
        //     for (int j = 0; j < (int)rDec.size(); j++) {
        //         remain[rDec[j]] -= 10. / float(rDec.size());
        //     }
        // }
        // minRem *= float(rand()) / float(RAND_MAX);
        myCharge[myOrder[i]] += minRem;

        for (int j = 0; j < (int)relevantWires.size(); j++) {
            remain[relevantWires[j]] -= minRem;
        }

        // remain[myLayer.cells[myOrder[i]].wireNums[0]] -= minRem;
        // remain[myLayer.cells[myOrder[i]].wireNums[1] + vDim - 1] -= minRem;
        // remain[myLayer.cells[myOrder[i]].wireNums[2] + vDim + lDim - 2] -= minRem;
    }
    vector <vector <float>> toReturn(2);
    toReturn[0] = myCharge;
    toReturn[1] = remain;

    //toReturn.push_back(myCharge);
    //toReturn.push_back(remain);
    return toReturn;
}





/* This method triese to just invert the matrix, sort of like wirecell */
vector <float> solveTrue(vector <float> mySig, vector <vector <int>> geoMat) {
	mat geo(geoMat.size(), geoMat[0].size());
	vec sig(mySig.size());
	for (int i = 0; i < (int)geoMat.size(); i++) {
		for (int j = 0; j < (int)geoMat[i].size(); j++) {
			geo(i, j) = geoMat[i][j];
		}
	}
	for (int i = 0; i < (int)mySig.size(); i++) {
		sig[i] = mySig[i];
	}
	vec trueSig = solve(geo, sig);
	vector <float> trueSignal;
	for (int i = 0; i < (int)geoMat[0].size(); i++) {
		trueSignal.push_back(trueSig[i]);
	}
	return trueSignal;
}

vector <vector <float>> solveTrue_multi(vector <vector <float>> allSigs, vector <vector <vector <int>>> allGeoMats) {
    vector <vector <float>> toReturn(allSigs.size());
    for (int i = 0; i < (int)allSigs.size(); i++) {
        if (allSigs[i].size() == 0) {
            continue;
        }
        toReturn[i] = solveTrue(allSigs[i], allGeoMats[i]);
    }
    return toReturn;
}



// Defines chi-squared computation for two vectors of floats.
// When the expected value is zero, just set to slightly
// greater than zero so that you don't throw nan but still
// get a very large penalty
float chiSqr(vector <float> expect, vector <float> observed) {
    assert (expect.size() == observed.size());
    float mySum = 0.;
    for (int i = 0; i < (int)expect.size(); i++) {
        if (expect[i] == 0) {
            expect[i] += 0.00001;
        }
        mySum += pow(expect[i] - observed[i], 2) / expect[i];
    }
    return mySum;
}

float chiSqr(vector <int> expect, vector <int> observed) {
    assert (expect.size() == observed.size());
    float mySum = 0.;
    for (int i = 0; i < (int)expect.size(); i++) {
        float e = (float)expect[i];
        if (e == 0.) {
            e += 0.00001;
        }
        mySum += pow(e - (float)observed[i], 2) / e;
    }
    return mySum;
}




float computeCost(vector <vector <float>> trial, vector <float> mySig, vector<vector <int>> geoMat, wireLayer myLayer) {
    float myCost = 0.;
    assert (mySig.size() == trial[1].size());
	//cout << "trying compute Cost\n";
    for (int i = 0; i < (int)mySig.size(); i++) {
        if (mySig[i] == 0.) {
            mySig[i] += 0.00001;
        }
        myCost += pow(trial[1][i],2)/mySig[i]; //1000 prefactor works nicely
        // myCost += pow(trial[1][i], 2);
    }
	//out << "computed the chi-squared\n";

	for (int i = 0; i < (int)myLayer.cells.size(); i++) {
		for (int j = i+1; j < (int)myLayer.cells.size(); j++) {
            // for now just use the average values of the vertices
            float x1 = 0.;
            float x2 = 0.;
            float y1 = 0.;
            float y2 = 0.;
            int s1 = myLayer.cells[i].vertices.size();
            int s2 = myLayer.cells[j].vertices.size();
            for (int k = 0; k < s1; k++) {
                x1 += myLayer.cells[i].vertices[k].x / float(s1);
                y1 += myLayer.cells[i].vertices[k].y / float(s1);
            }
            for (int k = 0; k < s2; k++) {
                x2 += myLayer.cells[j].vertices[k].x / float(s2);
                y2 += myLayer.cells[j].vertices[k].y / float(s2);
            }
			float xDiff = x1 - x2;
			float yDiff = y1 - y2;
			myCost -=  0.1 * trial[0][i] * trial[0][j] / (pow(xDiff,2) + pow(yDiff, 2));
		}
	}
	// cout << "computed packing component\n";

    // Chi squared component on the number of cells actived per wire
    // cout << "trial sizes [0] [1] " << trial[0].size() << "\t" << trial[1].size() << "\tGeomat sizes" << geoMat.size() << " " << geoMat[0].size() << endl;
    vector <int> oCells = numCellPerWire(trial[0], geoMat);
    vector <int> eCells;
    for (int i = 0; i < (int)mySig.size(); i++) {
        eCells.push_back((int)(mySig[i] / 10));
    }
    float cellCost = chiSqr(eCells, oCells);

    myCost += cellCost;

    return myCost;
}

// A naive random iteration function to minimize the cost.
// Repeats until 1000 trials have passed without the score improving
// At each iteration, it generates a new random charge distribution
// then computes the cost and tries to update the minimum cost distribution.
vector <float> solveCharge(vector <float> mySig, vector <vector <int>> geoMat, wireLayer myLayer) {
    int tryCount = 0;
    float minCost = numeric_limits<float>::infinity();
    vector <float> bestDistrib ((int)geoMat[0].size(), -1.);
    while (tryCount < 1000) {
		vector <vector <float>> myTrial = randChargeDistrib(mySig, geoMat, myLayer);
        float trialCost = computeCost(myTrial, mySig, geoMat, myLayer);
        if (trialCost < minCost) {
            minCost = trialCost;
            bestDistrib = myTrial[0];
            tryCount = 0;
        }
        tryCount++;
    }
    // cout << "Best solve charge cost: " << minCost << endl;
    return bestDistrib;
}

vector <vector <vector <float>>> randChargeDistrib_multi(vector <vector <float>> allSigs, vector <vector <vector <int>>> allGeoMats, vector <wireLayer> allLayers) {
    vector <vector <vector <float>>> allDistribs;
    for (int i = 0; i < (int)allSigs.size(); i++) {
        allDistribs.push_back(randChargeDistrib(allSigs[i], allGeoMats[i], allLayers[i]));
    }
    return allDistribs;
}

float computeCostMulti(vector <vector <vector <float>>> allTrials, vector <vector <float>> allSigs, vector <vector <vector <int>>> allGeoMats, vector <wireLayer> allLayers) {
    // cout << "All Trials sizes: " << allTrials.size() << "\t" << allTrials[0].size() << endl;
    // cout << "Sigs sizes: " << allSigs.size() << "\t" << allSigs[0].size() << endl;
    // cout << "geoMats sizes: " << allGeoMats.size() << "\t" << allGeoMats[0].size() << "\t" << allGeoMats[0][0].size() << endl;
    // cout << "Layers size: " << allLayers.size() << endl;
    assert(allTrials.size() == allSigs.size());
    assert(allSigs.size() == allGeoMats.size());
    assert(allGeoMats.size() == allLayers.size());

    float myCost = 0.;
    if (allTrials.empty()) {
        return myCost;
    }
    for (int i = 0; i < (int)allTrials.size(); i++) {
        if (allGeoMats[i].empty() || allSigs[i].empty() || allLayers[i].cells.empty()) {
            continue;
        }
        else {
            myCost += computeCost(allTrials[i], allSigs[i], allGeoMats[i], allLayers[i]);
        }
    }

    for (int i = 0; i < (int)allTrials.size() - 1; i++) {
        if (allLayers[i].cells.empty()) {
            continue;
        }
        for (int j = i+1; j < (int)allTrials.size(); j++) {
            // Do the spacial computation between time slices
            if (allLayers[j].cells.empty()) {
                continue;
            }
            for (int a = 0; a < (int)allLayers[i].cells.size(); a++) {
                if (allLayers[i].cells[a].vertices.empty()) {
                    continue;
                }
                for (int b = 0; b < (int)allLayers[j].cells.size(); b++) {
                    if (allLayers[j].cells[b].vertices.empty()) {
                        continue;
                    }
                    float x1 = 0.;
                    float x2 = 0.;
                    float y1 = 0.;
                    float y2 = 0.;
                    int s1 = (int)allLayers[i].cells[a].vertices.size();
                    int s2 = (int)allLayers[j].cells[b].vertices.size();
                    for (int k = 0; k < s1; k++) {
                        x1 += allLayers[i].cells[a].vertices[k].x / float(s1);
                        y1 += allLayers[i].cells[a].vertices[k].y / float(s1);
                    }
                    for (int k = 0; k < s2; k++) {
                        x2 += allLayers[j].cells[b].vertices[k].x / float(s2);
                        y2 += allLayers[j].cells[b].vertices[k].y / float(s2);
                    }
                    float xDiff = x1 - x2;
                    float yDiff = y1 - y2;
                    if ((xDiff != 0.) || (yDiff != 0.)) {
                        myCost -= 0.1 * allTrials[i][0][a] * allTrials[j][0][b] / (pow(xDiff, 2) + pow(yDiff, 2));
                    }
                }
            }
        }
    }
    return myCost;
}

vector <vector <float>> solveChargeMulti(vector <vector <float>> allSigs, vector <vector <vector <int>>> allGeoMats, vector <wireLayer> allLayers)
{
    int tryCount = 0;
    float minCost = numeric_limits<float>::infinity();
    // vector <vector <float>> bestDistrib;
    vector <vector <vector <float>>> initDistrib = randChargeDistrib_multi(allSigs, allGeoMats, allLayers);
    vector <vector <float>> bestDistrib;
    for (int i = 0; i < (int)initDistrib.size(); i++) {
        bestDistrib.push_back(initDistrib[i][0]);
    }
    minCost = computeCostMulti(initDistrib, allSigs, allGeoMats, allLayers);
    while (tryCount < 1000) {
        if (tryCount % 50 == 0) {
            cout << tryCount << endl;
        }
        vector <vector <vector <float>>> allTrials = randChargeDistrib_multi(allSigs, allGeoMats, allLayers);
        float trialCost = computeCostMulti(allTrials, allSigs, allGeoMats, allLayers);
        if (trialCost < minCost) {
            minCost = trialCost;
            bestDistrib.clear();
            for (int i = 0; i < (int)allTrials.size(); i++) {
                bestDistrib.push_back(allTrials[i][0]);
            }
            tryCount = 0;
        }
        tryCount++;
    }
    cout << "Solve charge multi, best soln cost: " << minCost << endl;
    return bestDistrib;
}




// Given a charge distribution, creates a new one by mutating the given
// one with some mutation frequency. Does not allow negative charges.
vector <vector <float>> mutateCharge(vector <vector <float>> myTrial, vector <vector <int>> geoMat) {
    vector <vector <float>> mutated = myTrial;
    for (int i = 0; i < (int) myTrial[0].size(); i++) {
        float myRand = float(rand()) / float (RAND_MAX);
        if (myRand < 0.1) {
            float randDir = float(rand()) / float (RAND_MAX);
            int change = 10;
            if ((randDir < 0.5) && (mutated[0][i] >= 10)) {
                change = -10;
            }
            mutated[0][i] += change;
            for (int j = 0; j < (int) geoMat.size(); j++) {
                if (geoMat[j][i] == 1) {
                    mutated[1][j] += change;
                }
            }
        }
    }
    return mutated;
}

vector <vector <vector <float>>> mutateChargeMulti(vector <vector <vector <float>>> myTrial, vector <vector <vector <int>>> allGeoMats) {
    vector <vector <vector <float>>> mutated;
    for (int i = 0; i < (int)myTrial.size(); i++) {
        mutated.push_back(mutateCharge(myTrial[i], allGeoMats[i]));
    }
    return mutated;
}

// Takes a vector of charges on cells, as well as the wire layer and geomatrix
// Uses these to construct a vector of "remaining charge" for each wire
vector <float> charge2remain(vector <float> myCharge, vector <float> mySig, vector <vector <int>> geoMat, wireLayer myLayer) {
    vector <float> remain = mySig;
    for (int i = 0; i < (int)myCharge.size(); i++) {
        if (myCharge[i] == 0) {
            continue;
        }
        for (int j = 0; j < (int) geoMat.size(); j++) {
            if (geoMat[j][i] == 1) {
                remain[j] -= myCharge[i];
            }
        }
    }
    return remain;
}


vector <vector <float>> charge2remain(vector <vector <float>> allCharge, vector <vector <float>> allSig, vector <vector <vector <int>>> allGeoMat, vector <wireLayer> allLayers) {
    vector <vector <float>> remain;
    for (int i = 0; i < (int)allSig.size(); i++) {
        remain.push_back(charge2remain(allCharge[i], allSig[i], allGeoMat[i], allLayers[i]));
    }
    return remain;
}

vector <vector <float>> zip1D(vector <float> a, vector <float> b) {
    assert(a.size() == b.size());
    vector <vector <float>> toRet;
    for (int i = 0; i < (int)a.size(); i++) {
        vector <float> temp;
        temp.push_back(a[i]);
        temp.push_back(b[i]);
        toRet.push_back(temp);
    }
    return toRet;
}

vector <vector <Path>> zip1D(vector <Path> a, vector <Path> b) {
    assert(a.size() == b.size());
    vector <vector <Path>> toRet;
    for (int i = 0; i < (int)a.size(); i++) {
        vector <Path> temp;
        temp.push_back(a[i]);
        temp.push_back(b[i]);
        toRet.push_back(temp);
    }
    return toRet;
}

vector <vector <vector <float>>> zip2D(vector <vector <float>> a, vector <vector <float>> b) {
    assert(a.size() == b.size());
    vector <vector <vector <float>>> toRet;
    for (int i = 0; i < (int)a.size(); i++) {
        vector <vector <float>> temp;
        temp.push_back(a[i]);
        temp.push_back(b[i]);
        toRet.push_back(temp);
    }
    return toRet;
}


// Randomly select a region in the charge vector to swap between the two
// solutions. Modify the remain vectors appropriately and then return
// the crossover reslts.
vector <vector <vector <float>>> crossoverCharge(vector <float> mySig, vector <vector <float>> t1, vector <vector <float>> t2, vector <vector <int>> geoMat) {
    assert (t1.size() == t2.size());
    assert (t1[0].size() == t2[0].size());
    if (t1[0].size() == 0) {
        vector <vector <vector <float>>> toReturn;
        toReturn.push_back(t1);
        toReturn.push_back(t2);
        return toReturn;
    }
    int startPosn = rand() % t1[0].size();
    int endPosn = (rand() % (t1[0].size() - startPosn)) + startPosn;
    for (int i = startPosn; i <= endPosn; i++) {
        float diff = t1[0][i] - t2[0][i];
        t2[0][i] += diff;
        t1[0][i] -= diff;
        for (int j = 0; j < (int)geoMat.size(); j++) {
            if (geoMat[j][i] == 1) {
                t1[1][j] -= diff;
                t2[1][j] += diff;
            }
        }
    }



    //
    // vector <float> r1 = mySig;
    // vector <float> r2 = mySig;
    // for (int i = 0; i < (int)mySig.size(); i++) {
    //     for (int j = 0; j < (int)t1[0].size(); j++) {
    //         if (geoMat[i][j] == 1) {
    //             r1[i] -= t1[0][j];
    //             r2[i] -= t2[0][j];
    //         }
    //     }
    // }
    // t1[1] = r1;
    // t2[1] = r2;
    vector <vector <vector <float>>> toReturn;
    toReturn.push_back(t1);
    toReturn.push_back(t2);
    return toReturn;
}

vector <vector <vector <vector <float>>>> crossoverChargeMulti(vector <vector <float>> allSigs, vector <vector <vector <float>>> t1, vector <vector <vector <float>>> t2, vector <vector <vector <int>>> allGeoMats) {
    assert (t1.size() == t2.size());
    assert (t1.size() == allSigs.size());
    vector <vector <vector <float>>> t1_new;
    vector <vector <vector <float>>> t2_new;
    //vector <vector <vector <vector <float>>>> toReturn((int)allSigs.size());
    for (int i = 0; i < (int)allSigs.size(); i++) {
        vector <vector <vector <float>>> temp = crossoverCharge(allSigs[i], t1[i], t2[i], allGeoMats[i]);
        t1_new.push_back(temp[0]);
        t2_new.push_back(temp[1]);
    }
    vector <vector <vector <vector <float>>>> toReturn;
    toReturn.push_back(t1);
    toReturn.push_back(t2);
    return toReturn;
}






/*
 * solveCharge_genetic: tries to find the best charge distribution via a genetic
 * algorithm approach. Seeds the population with 200 random distributions. Then
 * takes the top 40. Each iteration mutates with some probability and then re-evaluates.
 * After a fixed number of generations, the top scoring disribution is returned
 */
vector <float> solveCharge_genetic(vector <float> mySig, vector <vector <int>> geoMat, wireLayer myLayer) {
    int genCount = 0;
    map <vector <vector <float>> , float> genePool;
    for (int i = 0; i < 1000; i++) {
        vector <vector <float>> temp = randChargeDistrib(mySig, geoMat, myLayer);
        genePool[temp] = computeCost(temp, mySig, geoMat, myLayer);
    }

    vector<pair<vector<vector<float>>, float>> allPairs;
    for (auto i = genePool.begin(); i != genePool.end(); i++) {
        allPairs.push_back(*i);
    }
    sort(allPairs.begin(), allPairs.end(), [=](const pair<vector<vector<float>>, float>& a, const pair<vector <vector <float>>, float>& b) {
        return a.second < b.second;
    });
    genePool.clear();
    for (int i = 0; i < min(40, (int)allPairs.size()); i++) {
        genePool.insert(allPairs[i]);
    }

    while (genCount < 200) {
        vector<pair<vector<vector<float>>, float>> allPairs;
        for (auto i = genePool.begin(); i != genePool.end(); i++) {
            allPairs.push_back(*i);
            for (int j = 0; j < 5; j++) {
                vector<vector <float>> tempMutated = mutateCharge(i->first, geoMat);
                float myCost = computeCost(tempMutated, mySig, geoMat, myLayer);
                pair<vector <vector <float>>, float> tempPair = make_pair(tempMutated, myCost);
                allPairs.push_back(tempPair);
                map<vector<vector<float>>, float>::iterator k = genePool.begin();
                for (int a = 0; a < (int)(rand() % genePool.size()); a++) {
                    k++;
                }
                vector <vector <vector <float>>> tempCrossed = crossoverCharge(mySig, i->first, k->first, geoMat);
                float c1 = computeCost(tempCrossed[0], mySig, geoMat, myLayer);
                float c2 = computeCost(tempCrossed[1], mySig, geoMat, myLayer);
                pair<vector <vector <float>>, float> p1 = make_pair(tempCrossed[0], c1);
                pair<vector <vector <float>>, float> p2 = make_pair(tempCrossed[1], c2);
                allPairs.push_back(p1);
                allPairs.push_back(p2);
            }
        }
        sort(allPairs.begin(), allPairs.end(), [=](const pair<vector<vector<float>>, float>& a, const pair<vector <vector <float>>, float>& b) {
            return a.second < b.second;
        });
        genePool.clear();
        for (int i = 0; i < min(40, (int)allPairs.size()); i++) {
            genePool.insert(allPairs[i]);
        }
        genCount++;
    }
    allPairs.clear();
    for (auto i = genePool.begin(); i != genePool.end(); i++) {
        allPairs.push_back(*i);
    }
    sort(allPairs.begin(), allPairs.end(), [=](const pair<vector<vector<float>>, float>& a, const pair<vector <vector <float>>, float>& b) {
        return a.second < b.second;
    });
    return get<0>(allPairs[0])[0];
}

vector <vector <float>> solveChargeMulti_genetic(vector <vector <float>> allSigs, vector <vector <vector <int>>> allGeoMats, vector <wireLayer> allLayers) {
    int genCount = 0;
    map <vector <vector <vector <float>>>, float> genePool;
    for (int i = 0; i < 1000; i++) {
        if (i % 50 == 0) {
            cout << "Seeding: " << i << endl;
        }
        vector <vector <vector <float>>> temp = randChargeDistrib_multi(allSigs, allGeoMats, allLayers);
        genePool[temp] = computeCostMulti(temp, allSigs, allGeoMats, allLayers);
    }
    vector <pair<vector <vector <vector <float>>>, float>> allPairs;
    for (auto i = genePool.begin(); i != genePool.end(); i++) {
        allPairs.push_back(*i);
    }
    sort(allPairs.begin(), allPairs.end(), [=](const pair<vector <vector <vector <float>>>, float>&a, const pair <vector <vector <vector <float>>>, float> &b) {
        return a.second < b.second;
    });
    genePool.clear();
    for (int i = 0; i < min(40, (int)allPairs.size()); i++) {
        genePool.insert(allPairs[i]);
    }

    while (genCount < 10) {
        cout << "Generation: " << genCount << endl;
        vector<pair<vector<vector <vector <float>>>, float>> allPairs;
        for (auto i = genePool.begin(); i != genePool.end(); i++) {
            allPairs.push_back(*i);
            for (int j = 0; j < 5; j++) {
                vector<vector<vector <float>>> tempMutated = mutateChargeMulti(i->first, allGeoMats);
                float myCost = computeCostMulti(tempMutated, allSigs, allGeoMats, allLayers);
                pair<vector <vector <vector <float>>>, float> tempPair = make_pair(tempMutated, myCost);
                allPairs.push_back(tempPair);
                map<vector<vector<vector<float>>>, float>::iterator k = genePool.begin();
                for (int a = 0; a < (int)(rand() % genePool.size()); a++) {
                    k++;
                }
                vector <vector <vector <vector <float>>>> tempCrossed = crossoverChargeMulti(allSigs, i->first, k->first, allGeoMats);
                float c1 = computeCostMulti(tempCrossed[0], allSigs, allGeoMats, allLayers);
                float c2 = computeCostMulti(tempCrossed[1], allSigs, allGeoMats, allLayers);
                pair <vector <vector <vector <float>>>, float> p1 = make_pair(tempCrossed[0], c1);
                pair <vector <vector <vector <float>>>, float> p2 = make_pair(tempCrossed[1], c2);
                allPairs.push_back(p1);
                allPairs.push_back(p2);
            }
        }
        sort(allPairs.begin(), allPairs.end(), [=](const pair<vector <vector<vector<float>>>, float>& a, const pair<vector <vector <vector <float>>>, float>& b) {
            return a.second < b.second;
        });
        genePool.clear();

        for (int i = 0; i < min(40, (int)allPairs.size()); i++) {
            genePool.insert(allPairs[i]);
        }
        genCount++;
    }
    allPairs.clear();
    for (auto i = genePool.begin(); i != genePool.end(); i++) {
        allPairs.push_back(*i);
    }
    sort(allPairs.begin(), allPairs.end(), [=](const pair<vector <vector<vector<float>>>, float>& a, const pair<vector <vector <vector <float>>>, float>& b) {
        return a.second < b.second;
    });
    vector <vector <vector <float>>> bestTrial = get<0>(allPairs[0]);
    vector <vector <float>> toReturn;
    for (int i = 0; i < (int)bestTrial.size(); i++) {
        toReturn.push_back(bestTrial[i][0]);
    }
    cout << "Solve genetic multi, best soln cost: " << get<1>(allPairs[0]) << endl;

    return toReturn;
}



// Returns the average of the square difference between corresponding
// elements in two equally sized float vectors
float sumSqrDiff(vector <float> a, vector <float> b) {
    assert (a.size() == b.size());
    float mySum = 0.;
    for (int i = 0; i < (int)a.size(); i++) {
        mySum += pow(a[i] - b[i], 2);
    }
    return mySum / float(a.size());
}

// Sasme as the above, but defined for 2-D float vectors
float sumSqrDiff(vector <vector <float>> a, vector <vector <float>> b) {
    assert(a.size() == b.size());
    float mySum = 0.;
    for (int i = 0; i < (int)a.size(); i++) {
        mySum += sumSqrDiff(a[i], b[i]);
    }
    return mySum / float(a.size());
}


vector <vector <int>> reduceGeoMatrix(vector <vector <int>> myGeoMat, vector <float> signals) {
    vector <vector <int>> reduced;
    for (int i = 0; i < (int)signals.size(); i++) {
        if (signals[i] > 0.) {
            reduced.push_back(myGeoMat[i]);
        }
    }
    return reduced;
}
