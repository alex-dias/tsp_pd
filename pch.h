#ifndef PCH_H
#define PCH_H

class point2D {
public:
	float x, y;
	int id = 0;
	void printPoint();
	void printID();
};

struct distFloat {
	float minP, minD;
};

struct distTotal {
	float minTotal;
	int clusterIndex;
};

float pointDist(point2D, point2D);

void pairCostPrint(point2D, point2D, point2D, point2D);

#endif //PCH_H
