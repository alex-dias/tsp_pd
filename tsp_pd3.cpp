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
#include <Windows.h>

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

void printpairsPD(vector<vector<point2D>> pairsPD) {
	for (int i = 0; i < pairsPD.size(); i++) {
		for (std::vector<point2D>::iterator it = pairsPD[i].begin(); it != pairsPD[i].end(); ++it) {
			cout << returnID(*it) << '>';
		}
		cout << endl;
	}
}

// Dado um clusters e um ciclo trivial, essa funcao verifica qual a melhor forma e posicao de inserir esse ciclo trivial no cluster
// Alem disso, tambem atualiza o lastPairCluster com o ciclo recem inserido
void insertBestPosition(vector<point2D>& cluster, vector<point2D> trivial, vector<point2D> lastPairCluster) {

	float costP, costD, costPD, costPD2, 
		lowP = std::numeric_limits<float>::infinity(), 
		lowD = std::numeric_limits<float>::infinity(), 
		lowPD = std::numeric_limits<float>::infinity(),
		lowPD2 = std::numeric_limits<float>::infinity();
	int vizinho, insertP, insertD, insertPD, insertPD2;

	for (int i = 0; i < cluster.size(); i++) {

		if (i == cluster.size() - 1) {
			vizinho = 0;
		}
		else {
			vizinho = i;
			vizinho++;
		}

		costP = pointDist(trivial[0], cluster[i]) + pointDist(trivial[0], cluster[vizinho]) - pointDist(cluster[i], cluster[vizinho]);
		costD = pointDist(trivial[1], cluster[i]) + pointDist(trivial[1], cluster[vizinho]) - pointDist(cluster[i], cluster[vizinho]);
		costPD = pointDist(trivial[0], cluster[i]) + pointDist(trivial[1], cluster[vizinho]) + pointDist(trivial[0], trivial[1]) - pointDist(cluster[i], cluster[vizinho]);
		costPD2 = pointDist(trivial[0], cluster[vizinho]) + pointDist(trivial[1], cluster[i]) + pointDist(trivial[0], trivial[1]) - pointDist(cluster[i], cluster[vizinho]);

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
			insertPD2 = i;
		}
	}

	if (lowPD < lowP + lowD || insertP == insertD) {
		vector<point2D> vector2insert;
		vector2insert.push_back(trivial[0]);
		vector2insert.push_back(trivial[1]);
		cluster.insert(cluster.begin() + insertPD + 1, vector2insert.begin(), vector2insert.end());
		lastPairCluster = vector2insert;
	}
	else {
		if (lowPD2 < lowP + lowD) {
			vector<point2D> vector2insert;
			vector2insert.push_back(trivial[1]);
			vector2insert.push_back(trivial[0]);
			cluster.insert(cluster.begin() + insertPD2 + 1, vector2insert.begin(), vector2insert.end());
			lastPairCluster = vector2insert;
		}
		else {
			if (insertP < insertD) {
				cluster.insert(cluster.begin() + insertP + 1, trivial[0]);
				cluster.insert(cluster.begin() + insertD + 2, trivial[1]);
				lastPairCluster[0] = trivial[0];
				lastPairCluster[1] = trivial[1];
			}
			else {
				cluster.insert(cluster.begin() + insertD + 1, trivial[1]);
				cluster.insert(cluster.begin() + insertP + 2, trivial[0]);
				lastPairCluster[1] = trivial[1];
				lastPairCluster[0] = trivial[0];
			}
		}
	}
}

// Percorre a matriz de custo dos pares nao indexados e encontra o menor custo
void lowestCostMatrix(vector<distTotal> costs, int &i, int &j) {
	float lowestCost = std::numeric_limits<float>::infinity();

	for (int k = 0; k < costs.size(); k++) {
		if (costs[k].minTotal < lowestCost) {
			lowestCost = costs[k].minTotal;
			i = k;
			j = costs[k].clusterIndex;
		}
	}
}

// Retorna a distancia entre dois pares, que a menor distancia entre os quatro pontos que compoem os dois pares
float twoPairsCost(vector<point2D> pair1, vector<point2D> pair2) {
	vector<float> costs; // costPD, costDP, costPP, costDD;

	costs.push_back(pointDist(pair1[0], pair2[1])); // costPD
	costs.push_back(pointDist(pair1[1], pair2[0])); // costDP
	costs.push_back(pointDist(pair1[0], pair2[0])); // costPP
	costs.push_back(pointDist(pair1[1], pair2[1])); // costDD

	return *min_element(std::begin(costs), std::end(costs));
}

distFloat pointToPair(vector<point2D> pairsPD, vector<point2D> lastPairCluster) {
	distFloat dist;
	dist.minP = pointDist(pairsPD[0], lastPairCluster[0]); //PP
	float cost = pointDist(pairsPD[0], lastPairCluster[1]); //PD
	if (cost < dist.minP)
		dist.minP = cost;
	dist.minD = pointDist(pairsPD[1], lastPairCluster[0]); //DP
	cost = pointDist(pairsPD[1], lastPairCluster[1]); //DD
	if (cost < dist.minD)
		dist.minD = cost;
	
	return dist;
}

// Atualiza a distancia de todos os pares nao indexados em relacao a um cluster
void updatePairsCosts(vector<vector<point2D>> &pairsPD, vector<vector<distFloat>> &costs, vector<vector<point2D>> &lastPairCluster, vector<distTotal> &total, int clusterIndex) {
	for (int i = 0; i < costs.size(); i++) {
		distFloat cost = pointToPair(pairsPD[i], lastPairCluster[clusterIndex]);
		if (cost.minP < costs[i][clusterIndex].minP)
			costs[i][clusterIndex].minP = cost.minP;
		if (cost.minD < costs[i][clusterIndex].minD)
			costs[i][clusterIndex].minD = cost.minD;
		if (costs[i][clusterIndex].minP + costs[i][clusterIndex].minD < total[i].minTotal) {
			total[i].minTotal = costs[i][clusterIndex].minP + costs[i][clusterIndex].minD;
			total[i].clusterIndex = clusterIndex;
		}
	}
}

// Obtem os K clusters, segundo a heuristica
void clusterCalculator(vector<vector<point2D>> &clusters, vector<vector<point2D>> &lastPairCluster, vector<vector<point2D>> &pairsPD, int k, int c) {
	clusters.push_back(pairsPD[c]);
	lastPairCluster.push_back(pairsPD[c]);
	pairsPD.erase(pairsPD.begin() + c);

	vector<float> minDist;

	for (int i = 0; i < pairsPD.size(); i++) {
		minDist.push_back(twoPairsCost(clusters[0], pairsPD[i]));
	}

	for (int i = 1; i < k; i++) {
		float maxDist = minDist[0];
		int iMax = 0;

		for (int j = 0; j < pairsPD.size(); j++) {
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
	ofstream myfile;
	myfile.open("tests.csv");

	myfile << "NAME|DIMENSION|RESULT|CLUSTERS|TIME\n";

	WIN32_FIND_DATA data;
	HANDLE hFind = FindFirstFile("testInst\\*", &data);      // DIRECTORY

	FindNextFile(hFind, &data);
	FindNextFile(hFind, &data);

	int inst = 0;

	if (hFind != INVALID_HANDLE_VALUE) {
		do {

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
			inFile.open("testInst\\" + string(data.cFileName));
			if (!inFile) {
				cout << "Unable to open file" << endl;
				exit(1);
			}
			else {
				//cout << "File opened" << endl;
			}

			// coletando a dimencao do grafo
			while (getline(inFile, line)) {
				if (line.find(searchD) != string::npos) {
					DIM = line.substr(11);
					//cout << searchD << DIM << endl;
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

			// variavel de contagem de tempo
			clock_t tStart = clock();

			vector<vector<point2D>> bestCaseClusters;
			double bestCaseCost = std::numeric_limits<double>::infinity();

			int nClusters;

			int x = 0;

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

					// vetor que guarda o ultimo par inserido em um cluster
					vector<vector<point2D>> lastPairCluster;

					// associacao de cada ponto com seu respectivo par (pickup, delivery), que passam a ser considerados como grafos
					for (int i = 0; i < intDIM / 2; i++) {
						pairsPD[i].push_back(node[i]);
						pairsPD[i].push_back(node[i + intDIM / 2]);
					}

					clusterCalculator(clusters, lastPairCluster, pairsPD, k, c);

					// matriz de custo dos pares em relacao aos clusters
					vector<vector<distFloat>> pairsPDCosts;
					pairsPDCosts.resize(pairsPD.size());

					vector<distTotal> minDistTotal;
					minDistTotal.resize(pairsPD.size());
					for (int i = 0; i < pairsPD.size(); i++)
						minDistTotal[i].minTotal = std::numeric_limits<float>::infinity();

					for (int i = 0; i < pairsPD.size(); i++) {
						for (int j = 0; j < lastPairCluster.size(); j++) {
							//pointToPair(pairsPD[i], pairsPDCosts[i][j], lastPairCluster[j]);
							pairsPDCosts[i].push_back(pointToPair(pairsPD[i], lastPairCluster[j]));
							float cost = pairsPDCosts[i][j].minP + pairsPDCosts[i][j].minD;
							if (cost < minDistTotal[i].minTotal) {
								minDistTotal[i].minTotal = cost;
								minDistTotal[i].clusterIndex = j;
							}
						}
					}

					while (pairsPD.size() > 0) {

						int pairIndex, clusterIndex;

						lowestCostMatrix(minDistTotal, pairIndex, clusterIndex);

						insertBestPosition(clusters[clusterIndex], pairsPD[pairIndex], lastPairCluster[clusterIndex]);

						pairsPD.erase(pairsPD.begin() + pairIndex);
						pairsPDCosts.erase(pairsPDCosts.begin() + pairIndex);
						minDistTotal.erase(minDistTotal.begin() + pairIndex);

						updatePairsCosts(pairsPD, pairsPDCosts, lastPairCluster, minDistTotal, clusterIndex);
					}

					double caseCost = returnCycleCost(clusters);

					if (caseCost < bestCaseCost) {
						bestCaseCost = caseCost;
						nClusters = k;
					}
					//cout << x++ << endl;
				}

				//cout << endl;
			}

			cout << "Custo: " << fixed << bestCaseCost << endl;

			cout << inst << endl;

			inst++;

			// print do tempo de execucao
			cout  << "Time taken: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl << endl;

			myfile << string(data.cFileName) << "|" << intDIM / 2 << "|" << fixed << bestCaseCost << "|" << nClusters << "|" << ("%.5fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC) << "\n";

			inFile.close();
		} while (FindNextFile(hFind, &data));
		FindClose(hFind);
		myfile.close();
	}
}
