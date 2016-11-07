#include "point.hh"

class cell{
public:
    vector <vector <Point>> vLines;
    vector <vector <Point>> lLines;
    vector <vector <Point>> rLines;
    vector <Point> vertices;

    bool inCell(Point p) {
        return (betweenLines(p, vLines) && betweenLines(p, lLines) && betweenLines(p, rLines));
    }

    cell(int v [], int l [], int r []) : vLines(v), lLines(l), rLines(r){
        // assert the size of vLines, lLines, rLines
        assert(vLines.size() == 2);
        assert(lLines.size() == 2);
        assert(rLines.size() == 2);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                Point p1 = intersection(vLines[i][0], vLines[i][1], lLines[j][0], lLines[j][1]);
                Point p2 = intersection(vLines[i][0], vLines[i][1], rLines[j][0], rLines[j][1]);
                Point p3 = intersection(lLines[i][0], lLines[i][1], rLines[j][0], rLines[j][1]);
                if (inCell(p1)) {
                    vertices.push_back(p1);
                }
                if (inCell(p2)) {
                    vertices.push_back(p2);
                }
                if (inCell(p3)) {
                    vertices.push_back(p3);
                }
            }
        }
    }

};

bool betweenLines(Point p, vector< <vector Point>> lines) {
    assert(lines.size() == 2);
    // if any of the slopes are infinite, then we assume (v1) is left and (v2) is right so just check the x coord
    // if any slopes are 0, then just check the y coord

    float xDiff = lines[0][0].x - lines[0][1].x;
    float yDiff = lines[0][0].y - lines[0][1].y;
    float xDiff2 = lines[1][0].x - lines[1][1].x;
    if (xDiff == 0.) {
        if ((p.x < lines[0][0].x) || (p.x > lines[1][0].x)) {
            return false;
        }
    }
    if (yDiff == 0.) {
        if ((p.y < lines[0][0].y) || (p.y > lines[1][0].y)) {
            return false;
        }
    }
    else {
        float slope = yDiff / xDiff;
        if (p.y < lines[0][0].y - slope * xDiff) {
            return false;
        }
        if (p.y > lines[1][0].y - slope * xDiff2) {
            return false;
        }
    }
    return true;
}
