#include "events.hh"

using namespace std;

class wireArray
{
public:
	float pitch; // mm
	float angle; // angle with respect to x axis
	float height; // y dim
	float length; // x dim
	float z; // z posn
    vector<Point> startPoints;
    vector<Point> endPoints;
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
    wireArray vertical; // theta = Pi/2
    wireArray lSlant; // theta = Pi/4
    wireArray rSlant; // theta = -Pi/4
	float pitch;
	float length;
	float height;
	float z;
	vector <Point> grids;
	wireLayer(float h, float l, float z, float pitch) : pitch(pitch), length(l), height(h), z(z) {
		vertical = wireArray(pitch, 90., l, h, z);
		lSlant = wireArray(pitch, 45., l, h, z);
		rSlant = wireArray(pitch, -45., l, h, z);
	}
	vector <vector <int>> allCrossings(Point p1, Point p2) {
		vector <vector <int>> all3;
		all3.push_back(crossings(vertical, p1, p2));
		all3.push_back(crossings(lSlant, p1, p2));
		all3.push_back(crossings(rSlant, p1, p2));
		return all3;
	}
	vector <vector <int>> geoMatrix()
	{
		int vSize = vertical.startPoints.size();
		int lSize = lSlant.startPoints.size();
		int rSize = rSlant.startPoints.size();
		int numRows = vSize + lSize + rSize;
		vector <vector <int>> myMatrix;
		map<Point, vector<int>> gridPts;

		for (int i = 0; i < (int)lSlant.startPoints.size(); i++)
		{
			Point p1 = lSlant.startPoints[i];
			Point p2 = lSlant.endPoints[i];
			vector <int> wires = crossings(vertical, p1, p2);
			for (int j = 0; j < (int)wires.size(); j++) {
				Point refp1 = vertical.startPoints[wires[j]];
				Point refp2 = vertical.endPoints[wires[j]];
				Point intPoint = intersection(p1, p2, refp1, refp2);
				if (intPoint != Point(-1., -1., -1.)) {
					intPoint.round(2);
					// cout << "Point found: ";
					// intPoint.printPt();
					// cout << endl;
					if (gridPts.find(intPoint) == gridPts.end())
					{
						int wireRefs [3] = {wires[j], i, -1};
						vector<int> newEntry;
						newEntry.assign(wireRefs, wireRefs+3);
						gridPts[intPoint] = newEntry;
					}
					else {
						if (gridPts[intPoint][0] == -1) {
							gridPts[intPoint][0] = wires[j];
						}
						if (gridPts[intPoint][1] == -1) {
							gridPts[intPoint][1] = i;
						}
					}
				}
			}
		}
		for (int i = 0; i < (int)rSlant.startPoints.size(); i++)
		{
			Point p1 = rSlant.startPoints[i];
			Point p2 = rSlant.endPoints[i];
			vector <int> wires = crossings(vertical, p1, p2);
			for (int j = 0; j < (int)wires.size(); j++) {
				Point refp1 = vertical.startPoints[wires[j]];
				Point refp2 = vertical.endPoints[wires[j]];
				Point intPoint = intersection(p1, p2, refp1, refp2);
				if (intPoint != Point(-1., -1., -1.)) {
					intPoint.round(2);
					// cout << "Point found: ";
					// intPoint.printPt();
					// cout << endl;
					if (gridPts.find(intPoint) == gridPts.end())
					{
						int wireRefs [3] = {wires[j], -1, i};
						vector<int> newEntry;
						newEntry.assign(wireRefs, wireRefs+3);
						gridPts[intPoint] = newEntry;
					}
					else {
						if (gridPts[intPoint][0] == -1) {
							gridPts[intPoint][0] = wires[j];
						}
						if (gridPts[intPoint][2] == -1) {
							gridPts[intPoint][2] = i;
						}
					}
				}
			}
		}
		for (int i = 0; i < (int)rSlant.startPoints.size(); i++)
		{
			Point p1 = rSlant.startPoints[i];
			Point p2 = rSlant.endPoints[i];
			vector <int> wires = crossings(lSlant, p1, p2);
			for (int j = 0; j < (int)wires.size(); j++) {
				Point refp1 = lSlant.startPoints[wires[j]];
				Point refp2 = lSlant.endPoints[wires[j]];
				Point intPoint = intersection(p1, p2, refp1, refp2);
				if (intPoint != Point(-1., -1., -1.)) {
					intPoint.round(2);
					// cout << "Point found: ";
					// intPoint.printPt();
					// cout << endl;
					if ((intPoint.x <= 0) || (intPoint.y <= 0)) {
						continue;
					}
					if (gridPts.find(intPoint) == gridPts.end())
					{
						int wireRefs [3] = {-1, wires[j], i};
						vector<int> newEntry;
						newEntry.assign(wireRefs, wireRefs+3);
						gridPts[intPoint] = newEntry;
					}
					else {
						if (gridPts[intPoint][1] == -1) {
							gridPts[intPoint][1] = wires[j];
						}
						if (gridPts[intPoint][2] == -1) {
							gridPts[intPoint][2] = i;
						}
					}
				}
			}
		}
		int numCols = gridPts.size();
		for (int i = 0; i < numRows; i++) {
			vector <int> temp(numCols, 0);
			myMatrix.push_back(temp);
		}
		int count = 0;
		this->grids.clear();
		for (map<Point, vector<int>>::iterator j = gridPts.begin(); j != gridPts.end(); j++) {
			//j->first = key
			//j->second = vertical
			Point myPt = j->first;
			// cout << "Intersection Point: ";
			// myPt.printPt();
			// cout <<"\n";
			this->grids.push_back(myPt);
			vector<int> my3 = j->second;
			if (my3[0] != -1) {
				myMatrix[my3[0]][count] = 1;
			}
			if (my3[1] != -1) {
				myMatrix[my3[1] + vSize][count] = 1;
			}
			if (my3[0] != -1) {
				myMatrix[my3[2] + vSize + lSize][count] = 1;
			}
			count++;
		}
		return myMatrix;
	}

	bool inFrame(Point p) {
    	if ((p.x < 0) || (p.x > this->length)) {return false;}
    	if ((p.y < 0) || (p.y > this->height)) {return false;}
    	//if ((p.z < 0) || (p.z > zMax)) {return false;}
    	return true;
	}

vector <vector <int>> geoMatrix_test(vector <float> mySig) {
	vector<int> vPosn;
	vector<int> lPosn;
	vector<int> rPosn;
	int vSize = (int) this->vertical.startPoints.size();
	int lSize = (int) this->lSlant.startPoints.size();
	int rSize = (int) this->rSlant.startPoints.size();
	int i = 0;
	while (i < vSize) {
    	if (mySig[i] >= 1.) { // this is temporarily set to 1
        	// we care about this wire
        	vPosn.push_back(i);
    	}
    	i++;
	}
	while (i < vSize + lSize) {
    	if (mySig[i] >= 1.) { // this is temporarily set to 1
        	// we care about this wire
			lPosn.push_back(i-vSize);
    	}
    	i++;
	}
	while (i < vSize + lSize + rSize) {
    	if (mySig[i] >= 1.) { // this is temporarily set to 1
        	// we care about this wire
        	rPosn.push_back(i - vSize - lSize);
    	}
    	i++;
	}
	map<Point, vector<int>> gridPts;
	vector <vector <int>> myMatrix;
	for (int i = 0; i < (int)vPosn.size(); i++) {
    	// intersect with all of the lSlant ones
    	// intersect with all of the rSlant ones
    	// add intersection points to a map after checking within bounds and not (-1, -1, -1)
		for (int j = 0; j < (int)lPosn.size(); j++) {
        	Point pv_s = this->vertical.startPoints[vPosn[i]];
        	Point pv_e = this->vertical.endPoints[vPosn[i]];
        	Point pl_s = this->lSlant.startPoints[lPosn[j]];
        	Point pl_e = this->lSlant.endPoints[lPosn[j]];
        	Point intPt = intersection(pv_s, pv_e, pl_s, pl_e);
        	if ((intPt != Point(-1., -1., -1.)) && (inFrame(intPt))) {
            	intPt.round(2);
            	if (gridPts.find(intPt) == gridPts.end()) {
                	int wireRefs [3] = {i, j, -1};
                	vector <int> newEntry;
                	newEntry.assign(wireRefs, wireRefs+3);
                	gridPts[intPt] = newEntry;
            	}
            	else {
                	if (gridPts[intPt][0] == -1) {gridPts[intPt][0] = i;}
                	if (gridPts[intPt][1] == -1) {gridPts[intPt][1] = j;}
            	}
        	}
    	}
		for (int j = 0; j < (int)rPosn.size(); j++) {
        	Point pv_s = this->vertical.startPoints[vPosn[i]];
        	Point pv_e = this->vertical.endPoints[vPosn[i]];
        	Point pl_s = this->rSlant.startPoints[rPosn[j]];
        	Point pl_e = this->rSlant.endPoints[rPosn[j]];
        	Point intPt = intersection(pv_s, pv_e, pl_s, pl_e);
        	if ((intPt != Point(-1., -1., -1.)) && (inFrame(intPt))) {
            	intPt.round(2);
            	if (gridPts.find(intPt) == gridPts.end()) {
                	int wireRefs [3] = {i, -1, j};
                	vector <int> newEntry;
                	newEntry.assign(wireRefs, wireRefs+3);
                	gridPts[intPt] = newEntry;
            	}
            	else {
                	if (gridPts[intPt][0] == -1) {gridPts[intPt][0] = i;}
                	if (gridPts[intPt][2] == -1) {gridPts[intPt][2] = j;}
            	}
        	}
    	}
	}
	for (int i = 0; i < (int)lPosn.size(); i++) {
    	// intersect with all of the lSlant ones
    	// intersect with all of the rSlant ones
    	// add intersection points to a map after checking within bounds and not (-1, -1, -1)
		for (int j = 0; j < (int)rPosn.size(); j++) {
        	Point pv_s = this->lSlant.startPoints[lPosn[i]];
        	Point pv_e = this->lSlant.endPoints[lPosn[i]];
        	Point pl_s = this->rSlant.startPoints[rPosn[j]];
        	Point pl_e = this->rSlant.endPoints[rPosn[j]];
        	Point intPt = intersection(pv_s, pv_e, pl_s, pl_e);
        	if ((intPt != Point(-1., -1., -1.)) && (inFrame(intPt))) {
            	intPt.round(2);
            	if (gridPts.find(intPt) == gridPts.end()) {
                	int wireRefs [3] = {-1, i, j};
                	vector <int> newEntry;
                	newEntry.assign(wireRefs, wireRefs+3);
                	gridPts[intPt] = newEntry;
            	}
            	else {
                	if (gridPts[intPt][1] == -1) {gridPts[intPt][1] = i;}
                	if (gridPts[intPt][2] == -1) {gridPts[intPt][2] = j;}
            	}
        	}
    	}
    }
    int numCols = gridPts.size();
    for (int i = 0; i < (int)(rPosn.size() + lPosn.size() + vPosn.size()); i++) {
    	vector <int> tempVec(numCols, 0);
    	myMatrix.push_back(tempVec);
    }
    int count = 0;
    this->grids.clear();
    for (map<Point, vector<int>>::iterator j = gridPts.begin(); j != gridPts.end(); j++) {
    	Point myPt = j->first;
    	this->grids.push_back(myPt);
    	vector <int> my3 = j->second;
    	if (my3[0] != -1) {myMatrix[my3[0]][count] = 1;}
		if (my3[1] != -1) {myMatrix[my3[1] + vPosn.size()][count] = 1;}
    	if (my3[2] != -1) {myMatrix[my3[2] + vPosn.size() + lPosn.size()][count] = 1;}
    	count++;
    }
    return myMatrix;
}






	vector <float> signalVec(Path p) {
		int vDim = this->vertical.startPoints.size();
		int lDim = this->lSlant.startPoints.size();
		int rDim = this->rSlant.startPoints.size();
		vector <int> crossV = crossings(this->vertical, p.vertex, (p.vertex + p.path2vec()));
		vector <int> crossL = crossings(this->lSlant, p.vertex, (p.vertex + p.path2vec()));
		vector <int> crossR = crossings(this->rSlant, p.vertex, (p.vertex + p.path2vec()));
		vector <float> sigVec(vDim + lDim + rDim, 0.);
		for (int i = 0; i < (int)crossV.size(); i++) {
			sigVec[crossV[i]] += 10.;
		}
		for (int i = 0; i < (int)crossL.size(); i++) {
			sigVec[crossL[i] + vDim] += 10.;
		}
		for (int i = 0; i < (int)crossR.size(); i++) {
			sigVec[crossR[i] + vDim + lDim] += 10.;
		}
		return sigVec;
	}
};






vector <float> solveTrue(vector <float> mySig, vector <vector <int>> geoMat) {
	cout << geoMat.size();
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


vector <Cell> tiling(wireLayer w, float xMax, float yMax) {
	// for all the combos of 2 adjacent v, 2 adjacent l, and 2 adjacent r
	// check if the cell generated has at least 3 vertices
	// if so, append to the to return vector
	vector <Cell> allCells;
	// NOTE: THIS IS CERTAINLY AN INEFFICIENT VERSION, WILL IMPROVE LATER
	// NOTE: Can make a helper function to extract the desired two lines
	for (int i = 1; i < (int)w.vertical.startPoints.size() - 1; i++) {
		vector <vector <Point>> v;
		vector <Point> tempV;
		tempV.push_back(w.vertical.startPoints[i]);
		tempV.push_back(w.vertical.endPoints[i]);
		v.push_back(tempV);
		tempV.clear();
		tempV.push_back(w.vertical.startPoints[i+1]);
		tempV.push_back(w.vertical.endPoints[i+1]);
		v.push_back(tempV);
		tempV.clear();
		for (int j = 1; j < (int)w.lSlant.startPoints.size() - 1; j++) {
			vector <vector <Point>> l;
			vector <Point> tempL;
			tempL.push_back(w.lSlant.startPoints[j]);
			tempL.push_back(w.lSlant.endPoints[j]);
			l.push_back(tempL);
			tempL.clear();
			tempL.push_back(w.lSlant.startPoints[j+1]);
			tempL.push_back(w.lSlant.endPoints[j+1]);
			l.push_back(tempL);
			tempL.clear();
			for (int k = 1; k < (int)w.rSlant.startPoints.size() - 1; k++) {
				vector <vector <Point>> r;
				vector <Point> tempR;
				tempR.push_back(w.rSlant.startPoints[k]);
				tempR.push_back(w.rSlant.endPoints[k]);
				r.push_back(tempR);
				tempR.clear();
				tempR.push_back(w.rSlant.startPoints[k+1]);
				tempR.push_back(w.rSlant.endPoints[k+1]);
				r.push_back(tempR);
				tempR.clear();
				Cell myCell = Cell(v, l, r);
				if (myCell.vertices.size() >= 3) {
					int flag = 0;
					for (int c = 0; c < (int)myCell.vertices.size(); c++) {
						Point myCheck = myCell.vertices[c];
						if ((myCheck.x < 0) || (myCheck.x > xMax)) {
							flag = 1;
							break;
						}
						if ((myCheck.y < 0) || (myCheck.y > yMax)) {
							flag = 1;
							break;
						}
					}
					if (flag == 0) {
						allCells.push_back(myCell);
					}
				}
			}
		}
	}
	return allCells;
}
