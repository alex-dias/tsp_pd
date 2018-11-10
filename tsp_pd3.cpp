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
#include <algorithm>
#include <stdlib.h>

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

void printpairsPD(vector<vector<point2D>> pairsPD) {
	for (int i = 0; i < pairsPD.size(); i++) {
		for (std::vector<point2D>::iterator it = pairsPD[i].begin(); it != pairsPD[i].end(); ++it) {
			cout << returnID(*it) << '>';
		}
		cout << endl;
	}
}

void insertBestPosition(vector<point2D>& pairsPD, vector<point2D> trivial, vector<point2D> lastPairCluster) {

	float costP, costD, costPD, costPD2, 
		lowP = std::numeric_limits<float>::infinity(), 
		lowD = std::numeric_limits<float>::infinity(), 
		lowPD = std::numeric_limits<float>::infinity(),
		lowPD2 = std::numeric_limits<float>::infinity();
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
		costPD2 = pointDist(trivial[0], pairsPD[vizinho]) + pointDist(trivial[1], pairsPD[i]) + pointDist(trivial[0], trivial[1]) - pointDist(pairsPD[i], pairsPD[vizinho]);

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
		if (costPD2 < lowPD2) {
			lowPD2 = costPD2;
			insertPD = i;
		}
	}

	if (lowPD < lowP + lowD || insertP == insertD) {
		vector<point2D> vector2insert;
		vector2insert.push_back(trivial[0]);
		vector2insert.push_back(trivial[1]);
		pairsPD.insert(pairsPD.begin() + insertP + 1, vector2insert.begin(), vector2insert.end());
		lastPairCluster = vector2insert;
	}
	else {
		if (lowPD2 < lowP + lowD) {
			vector<point2D> vector2insert;
			vector2insert.push_back(trivial[1]);
			vector2insert.push_back(trivial[0]);
			pairsPD.insert(pairsPD.begin() + insertP + 1, vector2insert.begin(), vector2insert.end());
			lastPairCluster = vector2insert;
		}
		else {
			if (insertP < insertD) {
				pairsPD.insert(pairsPD.begin() + insertP + 1, trivial[0]);
				pairsPD.insert(pairsPD.begin() + insertD + 2, trivial[1]);
				lastPairCluster[0] = trivial[0];
				lastPairCluster[1] = trivial[1];
			}
			else {
				pairsPD.insert(pairsPD.begin() + insertD + 1, trivial[1]);
				pairsPD.insert(pairsPD.begin() + insertP + 2, trivial[0]);
				lastPairCluster[1] = trivial[1];
				lastPairCluster[0] = trivial[0];
			}
		}
	}

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

float twoPairsCost(vector<point2D> pair1, vector<point2D> pair2) {
	vector<float> costs; // costPD, costDP, costPP, costDD;

	costs.push_back(pointDist(pair1[0], pair2[1])); // costPD
	costs.push_back(pointDist(pair1[1], pair2[0])); // costDP
	costs.push_back(pointDist(pair1[0], pair2[0])); // costPP
	costs.push_back(pointDist(pair1[1], pair2[1])); // costDD

	return *min_element(std::begin(costs), std::end(costs));
}

void updatePairsCosts(vector<vector<point2D>> &pairsPD, vector<vector<float>> &costs, vector<vector<point2D>> &lastPairCluster, int clusterIndex) {
	for (int i = 0; i < costs.size(); i++) {
		costs[i][clusterIndex] = twoPairsCost(pairsPD[i], lastPairCluster[clusterIndex]);
	}
}

void clusterCalculator(vector<vector<point2D>> &clusters, vector<vector<point2D>> &lastPairCluster, vector<vector<point2D>> &pairsPD, int k, int c) {
	clusters.push_back(pairsPD[c]);
	pairsPD.erase(pairsPD.begin() + c);

	vector<float> minDist;

	for (int i = 0; i < pairsPD.size(); i++) {
		minDist.push_back(twoPairsCost(clusters[0], pairsPD[i]));
	}

	for (int i = 1; i < k; i++) {
		float maxDist = minDist[0];
		int iMax = 0;

		for (int j = 1; j < pairsPD.size(); j++) {
			if (minDist[j] > maxDist) {
				maxDist = minDist[j];
				iMax = j;
			}
		}

		clusters.push_back(*(pairsPD.begin() + iMax));
		lastPairCluster.push_back(*(pairsPD.begin() + iMax));
		pairsPD.erase(pairsPD.begin() + iMax);
		minDist.erase(minDist.begin() + iMax);

		for (int j = 0; j < pairsPD.size(); j++) {
			float newCost = twoPairsCost(*(clusters.end() - 1), pairsPD[j]);
			if (newCost < minDist[j])
				minDist[j] = newCost;
		}
	}
}

double returnCycleCost(vector<vector<point2D>> cycle) {
	double cost = 0;
	int x = 0;

	for (int i = 0; i < cycle.size(); i++) {
		for (std::vector<point2D>::iterator it = cycle[i].begin(); it != cycle[i].end(); ++it) {
			std::vector<point2D>::iterator vizinho;
			if (it == cycle[i].end() - 1) {
				vizinho = cycle[i].begin();
			}
			else {
				vizinho = it;
				vizinho++;
			}
			cost += pointDist(*it, *vizinho);
			x++;
		}
	}

	//cout << "Nodes: " << x << endl;

	return cost;
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

	// vetor dinamico de pontos, usada para receber os pontos sem seus respectivos pares
	vector<point2D> node;

	// estrutura dinamica de vetor de vetores
	// representa uma colecao de grafos, onde cada vetor eh um grafo
	vector<vector<point2D>> pairsPD;

	// clusters
	vector<vector<point2D>> clusters;

	// vetor que guarda os centroides dos clusters respectivamente
	vector<vector<point2D>> lastPairCluster;

	// leitura do arquivo
	inFile.open("rg-032-q-2x2-W0.4.0001.ccpdp");
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

	vector<vector<point2D>> bestCaseClusters;
	double bestCaseCost = std::numeric_limits<double>::infinity();

	vector<double> allCosts;

	for (int k = 2; k < intDIM / 2; k++) {

		//cout << "K: " << k << endl << endl;

		for (int c = 0; c < intDIM / 2; c++) {

			//cout << "C: " << c << endl;

			// estrutura dinamica de vetor de vetores
			// representa uma colecao de grafos, onde cada vetor eh um grafo
			vector<vector<point2D>> pairsPD;
			pairsPD.resize(intDIM / 2);

			// clusters
			vector<vector<point2D>> clusters;

			// vetor que guarda os centroides dos clusters respectivamente
			vector<vector<point2D>> lastPairCluster;

			// associacao de cada ponto com seu respectivo par (pickup, delivery), que passam a ser considerados como grafos
			for (int i = 0; i < intDIM / 2; i++) {
				pairsPD[i].push_back(node[i]);
				pairsPD[i].push_back(node[i + intDIM / 2]);
			}

			clusterCalculator(clusters, lastPairCluster, pairsPD, k, c);

			// 
			vector<vector<float>> pairsPDCosts;
			pairsPDCosts.resize(pairsPD.size());

			for (int i = 0; i < pairsPD.size(); i++) {
				for (int j = 0; j < lastPairCluster.size(); j++) {
					pairsPDCosts[i].push_back(twoPairsCost(pairsPD[i], lastPairCluster[j]));
				}
			}

			while (pairsPD.size() > 0) {

				int pairIndex, clusterIndex;

				lowestCostMatrix(pairsPDCosts, pairIndex, clusterIndex);

				insertBestPosition(clusters[clusterIndex], pairsPD[pairIndex], lastPairCluster[clusterIndex]);

				pairsPD.erase(pairsPD.begin() + pairIndex);
				pairsPDCosts.erase(pairsPDCosts.begin() + pairIndex);

				updatePairsCosts(pairsPD, pairsPDCosts, lastPairCluster, clusterIndex);
			}

			double cost = returnCycleCost(clusters);

			/*
			if (cost < bestCaseCost) {
				bestCaseClusters = clusters;
				bestCaseCost = cost;
			} */

			allCosts.push_back(cost);
		}

		//cout << endl;
	}

	/*
	for (int i = 0; i < bestCaseClusters.size(); i++) {
		for (int j = 0; j < bestCaseClusters[i].size(); j++) {
			cout << returnID(bestCaseClusters[i][j]) << ">";
		}
		cout << endl << endl;
	} */

	

	cout  << "Custo: " << *min_element(allCosts.begin(), allCosts.end()) << endl;

	// print do tempo de execucao
	cout << endl << "Time taken: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;

	inFile.close();
}
