#include "pch.h"
#include <iostream>
#include <istream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>	
#include <math.h>
#include <time.h>

using namespace std;

void point2D::printPoint() {
	cout << "X: " << x << " " << "Y: " << y;
};

void point2D::printID() {
	cout << "ID: " << id;
};

float pointDist(point2D point1, point2D point2) {
	return sqrt(pow((point1.x - point2.x), 2) + pow((point1.y - point2.y), 2));
};

void pairCostPrint(point2D point1, point2D point2, point2D point3, point2D point4) {
	float case1 = pointDist(point1, point3) + pointDist(point2, point4)
		- pointDist(point1, point2) - pointDist(point3, point4); // sem cruzar
	float case2 = pointDist(point1, point4) + pointDist(point2, point3)
		- pointDist(point1, point2) - pointDist(point3, point4); // cruzar

	cout << "Case 1: " << case1 << " Case 2: " << case2 << endl;
};

// case1 eh a comparacao direta e case 2 eh a comparacao cruzada
float pairCost(point2D point1, point2D point2, point2D point3, point2D point4, bool &isDirect) {
	float case1 = pointDist(point1, point3) + pointDist(point2, point4)
		- pointDist(point1, point2) - pointDist(point3, point4); // sem cruzar
	float case2 = pointDist(point1, point4) + pointDist(point2, point3)
		- pointDist(point1, point2) - pointDist(point3, point4); // cruzar

	// retorna o menor caso e muda o valor de isDirect
	if (case1 < case2) {
		isDirect = 1;
		return case1;
	}
	else {
		isDirect = 0;
		return case2;
	}
};

int returnID(point2D point1) {
	return point1.id;
};

// funcao para manipular e retorna o vetor que sera inserido
// a sequencia dos valores no vetor representa a ordem dos pontos de um grafo
// begin e end sao necessarios para simular a situacao de uma lista circular
vector<point2D> createNewVector(vector<point2D>::iterator it1, vector<point2D>::iterator itv1, vector<point2D>::iterator begin, vector<point2D>::iterator end, bool isDirect) {

	vector<point2D> newVector;

	// se a ligacao eh direta o novo vetor eh criado a partir it1, seguidos de seus vizinhos a ESQUERDA ateh itv1
	if (isDirect) {
		while (it1 != itv1) {
			newVector.push_back(*it1);
			if (it1 == begin)
				it1 = end;
			it1--;
		}
		newVector.push_back(*itv1);
	}
	// se a ligacao eh cruzada o novo vetor eh criado a partir itv1, seguidos de seus vizinhos a DIREITA ateh it1
	else {
		while (itv1 != it1) {
			newVector.push_back(*itv1);
			if (itv1 == end - 1) {
				itv1 = begin;
				break;
			}
			itv1++;
		}
		newVector.push_back(*it1);
	}

	return newVector;
};

double returnpairsPDCost(vector<vector<point2D>> pairsPD) {
	double cost = 0;

	for (int i = 0; i < pairsPD.size(); i++) {
		for (std::vector<point2D>::iterator it = pairsPD[i].begin(); it != pairsPD[i].end(); ++it) {
			std::vector<point2D>::iterator vizinho;
			if (it == pairsPD[i].end() - 1) {
				vizinho = pairsPD[i].begin();
			}
			else {
				vizinho = it;
				vizinho++;
			}
			cost += pointDist(*it, *vizinho);
		}
	}

	return cost;
};

void printpairsPD(vector<vector<point2D>> pairsPD) {
	for (int i = 0; i < pairsPD.size(); i++) {
		for (std::vector<point2D>::iterator it = pairsPD[i].begin(); it != pairsPD[i].end(); ++it) {
			cout << returnID(*it) << '>';
		}
		cout << endl;
	}
}

void insertBestPosition(vector<point2D>& pairsPD, vector<point2D> trivial) { // TO BE TESTED!!!!!!!!!!!!!!1

	float costP, costD, costPD, lowP = std::numeric_limits<float>::infinity(), lowD = std::numeric_limits<float>::infinity(), lowPD = std::numeric_limits<float>::infinity();
	int vizinho, insertP, insertD, insertPD;

	for (int i = 0; i < pairsPD.size(); i++) {

		if (i == pairsPD.size() - 1) {
			vizinho = 0;
		}
		else {
			vizinho = i;
			vizinho++;
		}

		costP = pointDist(trivial[0], pairsPD[i]) + pointDist(trivial[0], pairsPD[vizinho]) - pointDist(pairsPD[i], pairsPD[vizinho]);
		costD = pointDist(trivial[1], pairsPD[i]) + pointDist(trivial[1], pairsPD[vizinho]) - pointDist(pairsPD[i], pairsPD[vizinho]);
		costPD = pointDist(trivial[0], pairsPD[i]) + pointDist(trivial[1], pairsPD[vizinho]) + pointDist(trivial[0], trivial[1]) - pointDist(pairsPD[i], pairsPD[vizinho]);

		if (costP < lowP) {
			lowP = costP;
			insertP = i;
		}
		if (costD < lowD) {
			lowD = costD;
			insertD = i;
		}
		if (costPD < lowPD) {
			lowPD = costPD;
			insertPD = i;
		}
	}

	if (lowPD < lowP + lowD || insertP == insertD) {
		vector<point2D> vector2insert;
		vector2insert.push_back(trivial[0]);
		vector2insert.push_back(trivial[1]);
		pairsPD.insert(pairsPD.begin() + insertP + 1, vector2insert.begin(), vector2insert.end());
	}
	else {
		if (insertP < insertD) {
			pairsPD.insert(pairsPD.begin() + insertP + 1, trivial[0]);
			pairsPD.insert(pairsPD.begin() + insertD + 2 , trivial[1]);
		}
		else {
			pairsPD.insert(pairsPD.begin() + insertD + 1, trivial[1]);
			pairsPD.insert(pairsPD.begin() + insertP + 2, trivial[0]);
		}
	}
	
}

point2D twoPointsCentroid(point2D point1, point2D point2) {
	point2D centroid;

	centroid.x = (point1.x + point2.x) / 2;
	centroid.y = (point1.y + point2.y) / 2;

	return centroid;
}

void lowestCostMatrix(vector<vector<float>> costs, int &i, int &j) {
	float lowestCost = std::numeric_limits<float>::infinity();

	for (int k = 0; k < costs.size(); k++) {
		for (int l = 0; l < costs[k].size(); l++) {
			if (costs[k][l] < lowestCost) {
				lowestCost = costs[k][l];
				i = k;
				j = l;
			}
		}
	}
}

void updateCentroid(vector<vector<point2D>> &clusters, vector<point2D> &centroids, int index) {
	centroids[index] = twoPointsCentroid(centroids[index], *(clusters[index].end() - 1));
}

void updatePairsCosts(vector<vector<point2D>> &pairsPD, vector<vector<float>> &costs, vector<point2D> &centroids, int clusterIndex) {
	for (int i = 0; i < costs.size(); i++) {
		costs[i][clusterIndex] = pointDist(twoPointsCentroid(pairsPD[i][0], pairsPD[i][1]), centroids[clusterIndex]);
	}
}

int main()
{
	// variavel de contagem de tempo
	clock_t tStart = clock();

	// variaveis responsaveis pela leitura do arquivo e atribuicao das variaveis
	ifstream inFile;
	string line;
	string searchD = "DIMENSION: ";
	string searchN = "NODE_COORD_SECTION";
	string DIM;
	string NOD;
	int intDIM = 0;
	istringstream iss;
	int clustersN = 3;

	// vetor dinamico de pontos, usada para receber os pontos sem seus respectivos pares
	vector<point2D> node;

	// estrutura dinamica de vetor de vetores
	// representa uma colecao de grafos, onde cada vetor eh um grafo
	vector<vector<point2D>> pairsPD;

	// clusters
	vector<vector<point2D>> clusters;

	// vetor que guarda os centroides dos clusters respectivamente
	vector<point2D> centroids;

	// leitura do arquivo
	inFile.open("rg-256-q-2x2-W0.4.0001.ccpdp");
	if (!inFile) {
		cout << "Unable to open file" << endl;
		exit(1);
	}
	else {
		cout << "File opened" << endl;
	}

	// coletando a dimencao do grafo
	while (getline(inFile, line)) {
		if (line.find(searchD) != string::npos) {
			DIM = line.substr(11);
			cout << searchD << DIM << endl;
			intDIM = std::stoi(DIM);
			break;
		}
	}

	// com o valor da dimencao eh possivel saber que o tamanho da colecao de grafos (vetor de vetores)
	pairsPD.resize(intDIM / 2);

	// alimentacao da estrutura de pontos com os valores do arquivo
	while (getline(inFile, line)) {
		if (line.find(searchN) != string::npos) {
			for (int i = 0; i < intDIM; i++) {
				getline(inFile, line);
				iss.str(line);

				point2D point;
				iss >> point.id >> point.x >> point.y; // iss retira os espacos entre os valores

				node.push_back(point);

				iss.clear();
			}
		}
	}

	// associacao de cada ponto com seu respectivo par (pickup, delivery), que passam a ser considerados como grafos
	for (int i = 0; i < intDIM / 2; i++) {
		pairsPD[i].push_back(node[i]);
		pairsPD[i].push_back(node[i + intDIM / 2]);
	}

	for (int i = 0; i < clustersN; i++) {
		clusters.push_back(pairsPD[0]);
		centroids.push_back(twoPointsCentroid(clusters[i][0], clusters[i][1]));
		pairsPD.erase(pairsPD.begin());
	}

	// 
	vector<vector<float>> pairsPDCosts;
	pairsPDCosts.resize(pairsPD.size());

	for (int i = 0; i < pairsPD.size(); i++) {
		for (int j = 0; j < centroids.size(); j++) {
			pairsPDCosts[i].push_back(pointDist(twoPointsCentroid(pairsPD[i][0], pairsPD[i][1]), centroids[j]));
		}

	}

	while (pairsPD.size() > 0) {

		int pointIndex, clusterIndex;

		lowestCostMatrix(pairsPDCosts, pointIndex, clusterIndex);

		insertBestPosition(clusters[clusterIndex], pairsPD[pointIndex]);

		pairsPD.erase(pairsPD.begin() + pointIndex);
		pairsPDCosts.erase(pairsPDCosts.begin() + pointIndex);

		updateCentroid(clusters, centroids, clusterIndex);

		updatePairsCosts(pairsPD, pairsPDCosts, centroids, clusterIndex);
	}

	for (int i = 0; i < clusters.size(); i++) {
		for (int j = 0; j < clusters[i].size(); j++) {
			cout << returnID(clusters[i][j]) << ">";
		}
		cout << endl << endl;
	}

	// print do tempo de execucao
	cout << endl << "Time taken: " <<  (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;

	inFile.close();
}
