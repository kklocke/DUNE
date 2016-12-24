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
    wireLayer(float h, float l, float z, float pitch) : pitch(pitch), length(l), height(h), z(z) {
        vertical = wireArray(pitch, 90., l, h, z);
        lSlant = wireArray(pitch, 54.3, l, h, z);
        rSlant = wireArray(pitch, -54.3, l, h, z);
    }
    void operator=(const wireLayer &w) {
        this->pitch = w.pitch;
		this->length = w.length;
		this->height = w.height;
		this->z = w.z;
		this->vertical = w.vertical;
		this->lSlant = w.lSlant;
		this->rSlant = w.rSlant;
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
        cout << "STEP LENGTH: " << stepSize << endl;
        cout << "STEP NUMS: " << numSteps << endl;
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
        // for (map<vector <int>, int>::iterator j = cellSeen.begin(); j != cellSeen.end(); j++) {
        //     // build this cell
        //     vector <int> myNums = j->first;
        //     int nums[3] = {myNums[0], myNums[1], myNums[2]};
        //     vector <vector <Point>> v;
        //     vector <Point> tempV1;
        //     vector <Point> tempV2;
        //     tempV1.push_back(this->vertical.startPoints[myNums[0]]);
        //     tempV1.push_back(this->vertical.endPoints[myNums[0]]);
        //     tempV2.push_back(this->vertical.startPoints[myNums[0] + 1]);
        //     tempV2.push_back(this->vertical.endPoints[myNums[0] + 1]);
        //     v.push_back(tempV1);
        //     v.push_back(tempV2);
        //     vector <vector <Point>> l;
        //     vector <Point> tempL1;
        //     vector <Point> tempL2;
        //     tempL1.push_back(this->lSlant.startPoints[myNums[1]]);
        //     tempL1.push_back(this->lSlant.endPoints[myNums[1]]);
        //     tempL2.push_back(this->lSlant.startPoints[myNums[1] + 1]);
        //     tempL2.push_back(this->lSlant.endPoints[myNums[1] + 1]);
        //     l.push_back(tempL1);
        //     l.push_back(tempL2);
        //     vector <vector <Point>> r;
        //     vector <Point> tempR1;
        //     vector <Point> tempR2;
        //     tempR1.push_back(this->rSlant.startPoints[myNums[2]]);
        //     tempR1.push_back(this->rSlant.endPoints[myNums[2]]);
        //     tempR2.push_back(this->rSlant.startPoints[myNums[2] + 1]);
        //     tempR2.push_back(this->rSlant.endPoints[myNums[2] + 1]);
        //     r.push_back(tempR1);
        //     r.push_back(tempR2);
        //     Cell c(v, l , r, nums);
        //     if (c.touches[0][0]) {
        //         if (c.touches[0][1]) {
        //             sigVec[myNums[0]] += 5.;
        //             sigVec[myNums[0] + 1] += 5.;
        //         }
        //         else {
        //             sigVec[myNums[0]] += 10.;
        //         }
        //     }
        //     else if (c.touches[0][1]) {
        //         sigVec[myNums[0] + 1] += 10.;
        //     }
        //     if (c.touches[1][0]) {
        //         if (c.touches[1][1]) {
        //             sigVec[vDim + myNums[1]] += 5.;
        //             sigVec[vDim + myNums[1] + 1] += 5.;
        //         }
        //         else {
        //             sigVec[vDim + myNums[1]] += 10.;
        //         }
        //     }
        //     else if (c.touches[1][1]) {
        //         sigVec[vDim + myNums[1] + 1] += 10.;
        //     }
        //     if (c.touches[2][0]) {
        //         if (c.touches[2][1]) {
        //             sigVec[vDim + lDim + myNums[2]] += 5.;
        //             sigVec[vDim + lDim + myNums[2] + 1] += 5.;
        //         }
        //         else {
        //             sigVec[vDim + lDim + myNums[2]] += 10.;
        //         }
        //     }
        //     else if (c.touches[2][1]) {
        //         sigVec[vDim + lDim + myNums[2] + 1] += 10.;
        //     }
        // }
        return sigVec;
    }

    // vector <float> signalVec(Path p) {
    // 	int vDim = this->vertical.startPoints.size();
    // 	int lDim = this->lSlant.startPoints.size();
    // 	int rDim = this->rSlant.startPoints.size();
    // 	vector <int> crossV = crossings(this->vertical, p.vertex, (p.vertex + p.path2vec()));
    // 	vector <int> crossL = crossings(this->lSlant, p.vertex, (p.vertex + p.path2vec()));
    // 	vector <int> crossR = crossings(this->rSlant, p.vertex, (p.vertex + p.path2vec()));
    //     // cout << "vDim: " << vDim << endl;
    //     // cout << "lDim: " << lDim << endl;
    //     // cout << "rDim: " << rDim << endl;
    //     // cout << "crossV dim: " << crossV.size() << endl;
    //     // cout << "crossL dim: " << crossL.size() << endl;
    //     // cout << "crossR dim: " << crossR.size() << endl;
    //     vector <float> sigVec(vDim + lDim + rDim - 3, 0.);
    //     // cout << "test 1" << endl;
    //
    //     if (crossV.empty()) {
    //         // cout << "V EMPTY\n";
    //         int vPlace = -1;
    //         for (int i = 0; i < vDim; i++) {
    //             if (p.vertex.x < this->vertical.startPoints[i].x) {
    //                 vPlace = i - 1;
    //             }
    //         }
    //         if (vPlace != -1) {
    //             sigVec[vPlace] += 10;
    //         }
    //     }
    //     else {
    //         // if (crossV[0] == -1) {
    //         //     cout << "v neg\n";
    //         // }
    //         for (int i = 0; i < (int)crossV.size(); i++) {
    //     		sigVec[crossV[i]] += 10.;
    //     	}
    //         sigVec[crossV[0] - 1] += 10;
    //         sigVec[crossV[crossV.size() - 1] + 1] += 10;
    //         // cout << "wrote vert parts" << endl;
    //     }
    //
    //     if (crossL.empty()) {
    //         // cout << "L EMPTY\n";
    //         int lPlace = -1;
    //         for (int i = 0; i < lDim; i++) {
    //             float predY = this->lSlant.startPoints[i].y + tan(this->lSlant.angle)*(p.vertex.x - this->lSlant.startPoints[i].x);
    //             if (predY < p.vertex.y) {
    //                 lPlace = i - 1;
    //             }
    //         }
    //         if (lPlace != -1) {
    //             sigVec[lPlace + vDim - 1] += 10;
    //         }
    //     }
    //     else {
    //         // if (crossL[0] == -1) {
    //         //     cout << "l neg\n";
    //         // }
    //         for (int i = 0; i < (int)crossL.size(); i++) {
    //             sigVec[crossL[i] + vDim - 1] += 10.;
    //         }
    //         sigVec[crossL[0] + vDim - 1 - 1] += 10;
    //         sigVec[crossL[crossL.size() - 1] + vDim - 1 + 1] += 10;
    //         // cout << "wrote lSlant parts" << endl;
    //     }
    //
    //     if (crossR.empty()) {
    //         // cout << "R EMPTY\n";
    //         int rPlace = -1;
    //         for (int i = 0; i < rDim; i++) {
    //             float predY = this->rSlant.startPoints[i].y + tan(this->rSlant.angle) * (p.vertex.x - this->rSlant.startPoints[i].x);
    //             if (predY > p.vertex.y) {
    //                 rPlace = i - 1;
    //             }
    //         }
    //         if (rPlace != -1) {
    //             sigVec[rPlace + vDim + lDim - 2] += 10;
    //         }
    //     }
    //     else {
    //         // if (crossR[0] == -1) {
    //         //     cout << "r neg\n";
    //         // }
    //         for (int i = 0; i < (int)crossR.size(); i++) {
    //     		sigVec[crossR[i] + vDim + lDim - 2] += 10.;
    //     	}
    //         sigVec[crossR[0] + vDim + lDim - 2 - 1] += 10;
    //         sigVec[crossR[crossR.size() - 1] + vDim + lDim - 2 + 1] += 10;
    //         // cout << "wrote rSlant parts" << endl;
    //     }
    //
    //     // cout << "all done with sigvec" << endl;
    // 	return sigVec;
    // }
};

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

vector <vector <float>> randChargeDistrib(vector <float> mySig, vector <vector <int>> geoMat, wireLayer myLayer) {
    vector <float> remain = mySig;
    vector <float> myCharge ((int)geoMat[0].size(), 0.);
    int vDim = myLayer.vertical.startPoints.size();
    int lDim = myLayer.lSlant.startPoints.size();
    int rDim = myLayer.rSlant.startPoints.size();
    vector <int> myOrder((int)myLayer.cells.size(), 0);
    for (int i = 0; i < (int)myOrder.size(); i++) {
        myOrder[i] = i;
    }
    random_shuffle(myOrder.begin(), myOrder.end());
    for (int i = 0; i < (int)myLayer.cells.size(); i++) {
        float minRem = numeric_limits<float>::infinity();
        minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[0]]);
        minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[1] + vDim - 1]);
        minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[2] + vDim + lDim - 2]);
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
        remain[myLayer.cells[myOrder[i]].wireNums[0]] -= minRem;
        remain[myLayer.cells[myOrder[i]].wireNums[1] + vDim - 1] -= minRem;
        remain[myLayer.cells[myOrder[i]].wireNums[2] + vDim + lDim - 2] -= minRem;
    }
    vector <vector <float>> toReturn;
    toReturn.push_back(myCharge);
    toReturn.push_back(remain);
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
    vector <vector <float>> toReturn;
    for (int i = 0; i < (int)allSigs.size(); i++) {
        toReturn.push_back(solveTrue(allSigs[i], allGeoMats[i]));
    }
    return toReturn;
}

float computeCost(vector <vector <float>> trial, vector <float> mySig, vector<vector <int>> geoMat, wireLayer myLayer) {
    float myCost = 0.;
	//cout << "trying compute Cost\n";
    for (int i = 0; i < (int)mySig.size(); i++) {
        if (mySig[i] == 0.) {
            mySig[i] += 0.0001;
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

    return myCost;
}

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
    float myCost = 0.;
    if (allTrials.empty()) {
        return myCost;
    }
    for (int i = 0; i < (int)allTrials.size(); i++) {
        myCost += computeCost(allTrials[i], allSigs[i], allGeoMats[i], allLayers[i]);
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


vector <vector <vector <float>>> mutateChargeMulti(vector <vector <vector <float>>> myTrial, vector <vector <vector <int>>> allGeoMats) {
    vector <vector <vector <float>>> mutated;
    for (int i = 0; i < (int)myTrial.size(); i++) {
        mutated.push_back(mutateCharge(myTrial[i], allGeoMats[i]));
    }
    return mutated;
}


vector <vector <vector <float>>> crossoverCharge(vector <float> mySig, vector <vector <float>> t1, vector <vector <float>> t2, vector <vector <int>> geoMat) {
    int startPosn = rand() % t1[0].size();
    int endPosn = (rand() % (t1[0].size() - startPosn)) + startPosn;
    for (int i = startPosn; i <= endPosn; i++) {
        float temp = t1[0][i];
        t1[0][i] = t2[0][i];
        t2[0][i] = temp;
    }
    vector <float> r1 = mySig;
    vector <float> r2 = mySig;
    for (int i = 0; i < (int)mySig.size(); i++) {
        for (int j = 0; j < (int)t1[0].size(); j++) {
            if (geoMat[i][j] == 1) {
                r1[i] -= t1[0][j];
                r2[i] -= t2[0][j];
            }
        }
    }
    t1[1] = r1;
    t2[1] = r2;
    vector <vector <vector <float>>> toReturn;
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
    for (int i = 0; i < 40; i++) {
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
        for (int i = 0; i < 40; i++) {
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
    for (int i = 0; i < 40; i++) {
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
            }
        }
        sort(allPairs.begin(), allPairs.end(), [=](const pair<vector <vector<vector<float>>>, float>& a, const pair<vector <vector <vector <float>>>, float>& b) {
            return a.second < b.second;
        });
        genePool.clear();
        for (int i = 0; i < 40; i++) {
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

float sumSqrDiff(vector <float> a, vector <float> b) {
    assert (a.size() == b.size());
    float mySum = 0.;
    for (int i = 0; i < (int)a.size(); i++) {
        mySum += pow(a[i] - b[i], 2);
    }
    return mySum / float(a.size());
}

float sumSqrDiff(vector <vector <float>> a, vector <vector <float>> b) {
    assert(a.size() == b.size());
    float mySum = 0.;
    for (int i = 0; i < (int)a.size(); i++) {
        mySum += sumSqrDiff(a[i], b[i]);
    }
    return mySum / float(a.size());
}
