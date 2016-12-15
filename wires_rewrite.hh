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

        // generate all of the cells, then construct matrix accordingly

        for (int i = 1; i < vSize - 1; i++) {
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
        	for (int j = 1; j < lSize - 1; j++) {
        		if (hit[vSize + j] == 0) {
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
        		for (int k = 1; k < rSize - 1; k++) {
        			if (hit[vSize + lSize + k] == 0) {
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
        			Cell myCell = Cell(v, l, r, myNums);
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
                            vector <int> tempRow (numRows, 0);
                            tempRow[i - 1] = 1;
                            tempRow[j + vSize - 2] = 1;
                            tempRow[k + vSize + lSize - 3] = 1;
                            transposeMat.push_back(tempRow);
        				}
        			}
        		}
        	}
        }
        // Now transpose the matrix to get myMatrix
        int numCols = transposeMat.size();
        for (int i = 0; i < numRows; i++) {
            vector <int> tempRow;
            for (int j = 0; j < numCols; j++) {
                tempRow[j] = transposeMat[j][i];
            }
            myMatrix.push_back(tempRow);
        }
        return myMatrix;
    }

    vector <float> signalVec(Path p) {
        int vDim = this->vertical.startPoints.size();
        int lDim = this->lSlant.startPoints.size();
        int rDim = this->rSlant.startPoints.size();
        vector <float> sigVec(vDim + lDim + rDim, 0.);
        // go in short intervals?
        // break it up and find which cells I go through
        map <vector<int>, int> cellSeen;
        float stepSize = this->pitch / 3.;
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
                if (this->vertical.startPoints[j].x < currentPt.x) {
                    vPlace = j;
                    break;
                }
            }
            for (int j = 0; j < lDim; j++) {
                float slope = tan(this->lSlant.angle);
                float predY = this->lSlant.startPoints[j].y + slope * (currentPt.x - this->lSlant.startPoints[j].x);
                if (predY < currentPt.y) {
                    lPlace = j - 1;
                    break;
                }
            }

            for (int j = 0; j < rDim; j++) {
                float slope = tan(this->rSlant.angle);
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
            // build this cell
            vector <int> myNums = j->first;
            int nums[3] = {myNums[0], myNums[1], myNums[2]};
            vector <vector <Point>> v;
            vector <Point> tempV1;
            vector <Point> tempV2;
            tempV1.push_back(this->vertical.startPoints[myNums[0]]);
            tempV1.push_back(this->vertical.endPoints[myNums[0]]);
            tempV2.push_back(this->vertical.startPoints[myNums[0] + 1]);
            tempV2.push_back(this->vertical.endPoints[myNums[0] + 1]);
            v.push_back(tempV1);
            v.push_back(tempV2);
            vector <vector <Point>> l;
            vector <Point> tempL1;
            vector <Point> tempL2;
            tempL1.push_back(this->lSlant.startPoints[myNums[1]]);
            tempL1.push_back(this->lSlant.endPoints[myNums[1]]);
            tempL2.push_back(this->lSlant.startPoints[myNums[1] + 1]);
            tempL2.push_back(this->lSlant.endPoints[myNums[1] + 1]);
            l.push_back(tempL1);
            l.push_back(tempL2);
            vector <vector <Point>> r;
            vector <Point> tempR1;
            vector <Point> tempR2;
            tempR1.push_back(this->rSlant.startPoints[myNums[2]]);
            tempR1.push_back(this->rSlant.endPoints[myNums[2]]);
            tempR2.push_back(this->rSlant.startPoints[myNums[2] + 1]);
            tempR2.push_back(this->rSlant.endPoints[myNums[2] + 1]);
            r.push_back(tempR1);
            r.push_back(tempR2);
            Cell c(v, l , r, nums);
            if (c.touches[0][0]) {
                if (c.touches[0][1]) {
                    sigVec[myNums[0]] += 5.;
                    sigVec[myNums[0] + 1] += 5.;
                }
                else {
                    sigVec[myNums[0]] += 10.;
                }
            }
            else if (c.touches[0][1]) {
                sigVec[myNums[0] + 1] += 10.;
            }
            if (c.touches[1][0]) {
                if (c.touches[1][1]) {
                    sigVec[vDim + myNums[1]] += 5.;
                    sigVec[vDim + myNums[1] + 1] += 5.;
                }
                else {
                    sigVec[vDim + myNums[1]] += 10.;
                }
            }
            else if (c.touches[1][1]) {
                sigVec[vDim + myNums[1] + 1] += 10.;
            }
            if (c.touches[2][0]) {
                if (c.touches[2][1]) {
                    sigVec[vDim + lDim + myNums[2]] += 5.;
                    sigVec[vDim + lDim + myNums[2] + 1] += 5.;
                }
                else {
                    sigVec[vDim + lDim + myNums[2]] += 10.;
                }
            }
            else if (c.touches[2][1]) {
                sigVec[vDim + lDim + myNums[2] + 1] += 10.;
            }
        }
        return sigVec;
    }
};

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
        vector <int> vDec;
        vector <int> lDec;
        vector <int> rDec;
        if (myLayer.cells[myOrder[i]].touches[0][0]) {
            minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[0]]);
            vDec.push_back(myLayer.cells[myOrder[i]].wireNums[0]);
        }
        if (myLayer.cells[myOrder[i]].touches[0][1]) {
            minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[0] + 1]);
            vDec.push_back(myLayer.cells[myOrder[i]].wireNums[0] + 1);
        }
        if (myLayer.cells[myOrder[i]].touches[1][0]) {
            minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[1] + vDim]);
            lDec.push_back(myLayer.cells[myOrder[i]].wireNums[1] + vDim);
        }
        if (myLayer.cells[myOrder[i]].touches[1][1]) {
            minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[1] + vDim + 1]);
            lDec.push_back(myLayer.cells[myOrder[i]].wireNums[1] + vDim + 1);
        }
        if (myLayer.cells[myOrder[i]].touches[2][0]) {
            minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[2] + vDim + lDim]);
            rDec.push_back(myLayer.cells[myOrder[i]].wireNums[2] + vDim + lDim);
        }
        if (myLayer.cells[myOrder[i]].touches[2][1]) {
            minRem = min(minRem, remain[myLayer.cells[myOrder[i]].wireNums[2] + vDim + lDim + 1]);
            rDec.push_back(myLayer.cells[myOrder[i]].wireNums[2] + vDim + lDim + 1);
        }
        if (minRem > 0) {
            // put charge in the cells
            myCharge[myOrder[i]] = 10.;
            // decrement the wires
            for (int j = 0; j < (int)vDec.size(); j++) {
                remain[vDec[j]] -= 10. / float(vDec.size());
            }
            for (int j = 0; j < (int)lDec.size(); j++) {
                remain[lDec[j]] -= 10. / float(lDec.size());
            }
            for (int j = 0; j < (int)rDec.size(); j++) {
                remain[rDec[j]] -= 10. / float(rDec.size());
            }
        }
    }
    vector <vector <float>> toReturn;
    toReturn.push_back(myCharge);
    toReturn.push_back(remain);
    return toReturn;
}
