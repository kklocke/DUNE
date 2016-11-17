#include "point.hh"

using namespace std;

bool betweenLines(Point p, vector <vector <Point>> lines);

class Cell{
public:
    vector <vector <Point>> vLines;
    vector <vector <Point>> lLines;
    vector <vector <Point>> rLines;
    vector <Point> vertices;
    int wireNums [3];



    Cell(vector <vector <Point>> v, vector <vector <Point>> l, vector <vector <Point>> r, int w[3]) : vLines(v), lLines(l), rLines(r){
        // assert the size of vLines, lLines, rLines
        assert(vLines.size() == 2);
        assert(lLines.size() == 2);
        assert(rLines.size() == 2);

        wireNums[0] = w[0];
        wireNums[1] = w[1];
        wireNums[2] = w[2];
        // ordering of the vertices is really important for tile coloring for display

        Point p1 = intersection(lLines[0][0], lLines[0][1], rLines[1][0], rLines[1][1]);
        if (betweenLines(p1, vLines)) {
            vertices.push_back(p1);
        }
        Point p2 = intersection(vLines[0][0], vLines[0][1], rLines[1][0], rLines[1][1]);
        if (betweenLines(p2, lLines)) {
            vertices.push_back(p2);
        }
        Point p3 = intersection(vLines[0][0], vLines[0][1], lLines[0][0], lLines[0][1]);
        if (betweenLines(p3, rLines)) {
            vertices.push_back(p3);
        }
        Point p4 = intersection(lLines[0][0], lLines[0][1], rLines[0][0], rLines[0][1]);
        if (betweenLines(p4, vLines)) {
            vertices.push_back(p4);
        }
        Point p5 = intersection(vLines[0][0], vLines[0][1], rLines[0][0], rLines[0][1]);
        if (betweenLines(p5, lLines)) {
            vertices.push_back(p5);
        }
        Point p6 = intersection(lLines[1][0], lLines[1][1], vLines[0][0], vLines[0][1]);
        if (betweenLines(p6, rLines)) {
            vertices.push_back(p6);
        }
        Point p7 = intersection(rLines[0][0], rLines[0][1], lLines[1][0], lLines[1][1]);
        if (betweenLines(p7, vLines)) {
            vertices.push_back(p7);
        }
        Point p8 = intersection(rLines[0][0], rLines[0][1], vLines[1][0], vLines[1][1]);
        if (betweenLines(p8, lLines)) {
            vertices.push_back(p8);
        }
        Point p9 = intersection(lLines[1][0], lLines[1][1], vLines[1][0], vLines[1][1]);
        if (betweenLines(p9, rLines)) {
            vertices.push_back(p9);
        }
        Point p10 = intersection(lLines[1][0], lLines[1][1], rLines[1][0], rLines[1][1]);
        if (betweenLines(p10, vLines)) {
            vertices.push_back(p10);
        }
        Point p11 = intersection(vLines[1][0], vLines[1][1], rLines[1][0], rLines[1][1]);
        if (betweenLines(p11, lLines)) {
            vertices.push_back(p11);
        }
        Point p12 = intersection(lLines[0][0], lLines[0][1], vLines[1][0], vLines[1][1]);
        if (betweenLines(p12, rLines)) {
            vertices.push_back(p12);
        }

        // for (int i = 0; i < 2; i++) {
        //     for (int j = 0; j < 2; j++) {
        //         Point p1 = intersection(vLines[i][0], vLines[i][1], lLines[1-j][0], lLines[1-j][1]);
        //         Point p2 = intersection(vLines[i][0], vLines[i][1], rLines[j][0], rLines[j][1]);
        //         Point p3 = intersection(lLines[1-i][0], lLines[1-i][1], rLines[j][0], rLines[j][1]);
        //         if (betweenLines(p1, rLines)) {
        //             vertices.push_back(p1);
        //         }
        //         if (betweenLines(p2, lLines)) {
        //             vertices.push_back(p2);
        //         }
        //         if (betweenLines(p3, vLines)) {
        //             vertices.push_back(p3);
        //         }
        //     }
        // }
    }
    // bool inCell(Point p, int dir) {
    //     // dir = 0 means we want to check
    //     bool tryV = betweenLines(p, vLines);
    //     bool tryL = betweenLines(p, lLines);
    //     bool tryR = betweenLines(p, rLines);
    //     cout <<  "In Cell: " << int(tryV) << " " << int(tryL) << " " << int(tryR) << endl;
    //     return (tryV && tryL && tryR);
    // }

};

bool betweenLines(Point p, vector< vector <Point>> lines) {
    assert(lines.size() == 2);
    // if any of the slopes are infinite, then we assume (v1) is left and (v2) is right so just check the x coord
    // if any slopes are 0, then just check the y coord
    float xDiff = lines[0][0].x - lines[0][1].x;
    float yDiff = lines[0][0].y - lines[0][1].y;
    if (xDiff == 0.) {
        // cout << "Infinite slope, left v, right v, my x" << "\n\t" << lines[0][0].x << "\t" << lines[1][0].x << "\t" << p.x << endl;
        // cout << "\tEvals: " << int((p.x < lines[0][0].x) || (p.x > lines[1][0].x)) << endl;
        if ((p.x < min(lines[0][0].x, lines[1][0].x)) || (p.x > max(lines[1][0].x, lines[0][1].x))) {
            return false;
        }
    }
    else if (yDiff == 0.) {
        if ((p.y < min(lines[0][0].y, lines[1][0].y)) || (p.y > max(lines[1][0].y, lines[0][0].y))) {
            return false;
        }
    }
    else {
        float slope = yDiff / xDiff;
        if (slope > 0) {
            if (p.y < lines[1][0].y + slope * (p.x - lines[1][0].x) - 0.001) {
                return false;
            }
            if (p.y > lines[0][0].y + slope * (p.x - lines[0][0].x) + 0.001) {
                return false;
            }
        }
        else if (slope < 0) {
            if (p.y < lines[0][0].y + slope * (p.x - lines[0][0].x) - 0.001) {
                return false;
            }
            if (p.y > lines[1][0].y + slope * (p.x - lines[1][0].x) + 0.001) {
                return false;
            }
        }

    }
    return true;
}
