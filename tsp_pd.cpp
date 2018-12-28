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
#include <windows.h>

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
				if (itv1 == it1)
					break;
				newVector.push_back(*itv1);
			}
			itv1++;
		}
		newVector.push_back(*it1);
	}

	return newVector;
};

double returnCycleCost(vector<vector<point2D>> cycle) {
	double cost = 0;

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
		}
	}

	return cost;
};

void printCycle(vector<vector<point2D>> cycle) {
	int x = 0;
	for (int i = 0; i < cycle.size(); i++) {
		for (std::vector<point2D>::iterator it = cycle[i].begin(); it != cycle[i].end(); ++it) {
			cout << returnID(*it) << '>';
			x++;
		}
		cout << endl;
	}
	cout << endl << x << endl;
}

int main() {

	ofstream myfile;
	myfile.open("tests.csv");

	myfile << "NAME|DIMENSION|RESULT|TIME\n";

	WIN32_FIND_DATA data;
	HANDLE hFind = FindFirstFile("testInst\\*", &data);      // DIRECTORY

	FindNextFile(hFind, &data);
	FindNextFile(hFind, &data);

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
			vector<vector<point2D>> cycle;

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
			cycle.resize(intDIM / 2);

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
				cycle[i].push_back(node[i]);
				cycle[i].push_back(node[i + intDIM / 2]);
			}

			vector<vector<point2D>> bestCase = cycle;
			double bestCaseCost = returnCycleCost(bestCase);

			// variavel de contagem de tempo
			clock_t tStart = clock();

			// aqui comeca a manipulacao em si
			// o objetivo eh tornar todos os grafos em um soh, por isso a condicao de parada eh quando a colecao de grafos tiver soh um grafo
			while (cycle.size() > 1) {

				/*
				   vizinho1 e vizinho 2 sao os proximos pontos em relacao aos iteradores it1 e it2 respectivamente
				   lowest1 e lowest2 guardam os dois iteradores (it1, it2) respectivos a comparacao de menor custo
				   lowest2V eh o proximo ponto em relacao a lowest2, lowest1V nao eh salvo pq o vetor resultante sempre eh inserido a direita de lowest1, indenpendente de lowest1V
				   eraseP eh uma variavel auxiliar para remover um vetor da colecao de vetores
				   lowestCost guarda o menor custo encontrado depois das comparacoes, por isso eh iniciado com o valor 'infinito'
				   cost armazena o custo da comparacao
				   isDirect armazena se o custo retornado foi em relacao a um comparacao cruzada ou nao, iDirectLow armazena a direcao da comparacao de menor custo
				   insertI eh uma variavel auxiliar que guardar a posicao do vetor no qual sera realizado a insercao
				*/
				vector<point2D>::iterator vizinho1, vizinho2, lowest1, lowest2, lowest2V, lowBegin, lowEnd;
				vector<vector<point2D>>::iterator eraseP;
				float lowestCost = std::numeric_limits<float>::infinity();
				float cost;
				bool isDirect, isDirectLow;
				int insertI;
				double cycleCost;

				for (int i = 0; i < cycle.size() - 1; ++i) // - 1 pq o ultimo grafo nao compara com ninguem
				{
					//cout << "------ EXTERNO " << i + 1 << " -------" << endl;
					for (std::vector<point2D>::iterator it1 = cycle[i].begin(); it1 != cycle[i].end(); ++it1) { // it1 eh um ponteiro que percorre o vetor de vetores
						// eh preciso tratar caso it1 seja a ultima posicao do vetor, de modo que seu vizinho devera ser a primeira posicao do vetor
						if (it1 == cycle[i].end() - 1) { // -1 pq a funcao end() retornar a posicao seguinte a ultima posicao do vetor
							vizinho1 = cycle[i].begin();
						}
						else {
							vizinho1 = it1;
							vizinho1++;
						}
						for (int j = i + 1; j < cycle.size(); ++j) { // j sempre diminui em relacao a i, para evitar comparacoes ja realizadas
							//cout << "------ INTERNO " << j << " -------" << endl;
							for (std::vector<point2D>::iterator it2 = cycle[j].begin(); it2 != cycle[j].end(); ++it2) { // it2 eh um ponteiro que percorre cada vetor da colecao
								vizinho2 = it2;
								// mesmo caso de it1
								if (vizinho2 == cycle[j].end() - 1) {
									vizinho2 = cycle[j].begin();
								}
								else {
									vizinho2++;
								}
								// calculo do custo entre os pares de pontos
								cost = pairCost(*it1, *vizinho1, *it2, *vizinho2, isDirect);

								// prints para auxiliar o debug
								//cout << returnID(*it1) << " - " << returnID(*vizinho1) << " X " << returnID(*it2) << " - " << returnID(*vizinho2) << endl;
								//pairCostPrint(*it1, *vizinho1, *it2, *vizinho2);
								//cout << endl << endl;

								// atribuicoes caso o custo encontrado seja o menor ateh agora
								if (cost < lowestCost) {
									lowestCost = cost;
									lowest1 = it1;
									lowest2 = it2;
									lowest2V = vizinho2;
									insertI = i;
									isDirectLow = isDirect;
									lowBegin = cycle[j].begin();
									lowEnd = cycle[j].end();
									eraseP = cycle.begin() + j;
								}
							}
						}
					}
				}

				// prints para auxiliar o debug
				//cout << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
				//cout << endl << returnID(*lowest1) << " X " << returnID(*lowest2) << " - " << returnID(*lowest2V) << endl;
				//cout << "Custo: " << lowestCost << endl;
				//cout << "Direto? " << isDirectLow << endl;

				// criacao de um vetor baseado no menor custo
				vector<point2D> vectorToInsert = createNewVector(lowest2, lowest2V, lowBegin, lowEnd, isDirectLow);

				// insercao do vetor criado no vetor alvo (lowest1)
				// +1 pq a funcao insert() insere na posicao anterior ao ponteiro dado
				cycle[insertI].insert(lowest1 + 1, vectorToInsert.begin(), vectorToInsert.end());

				// o vetor que foi inserido eh apagado
				cycle.erase(eraseP);

				cycleCost = returnCycleCost(cycle);

				if (cycleCost < bestCaseCost) {
					bestCase = cycle;
					bestCaseCost = cycleCost;
				}

				// prints para auxiliar o debug

				/*
				for (std::vector<point2D>::iterator it1 = cycle[insertI].begin(); it1 != cycle[insertI].end(); ++it1) {
					cout << returnID(*it1) << ">";
				}

				cout << endl << endl;
				*/

			}

			//cout << "Best case : " << endl << endl;

			//printCycle(bestCase);

			cout << fixed << endl << "Custo: " << bestCaseCost << endl;

			// print do tempo de execucao
			printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

			cout << endl << endl;

			myfile << string(data.cFileName) << "|" << intDIM/2 << "|" << fixed << bestCaseCost << "|" << ("%.5fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC) << "\n";

			inFile.close();

		} while (FindNextFile(hFind, &data));
		FindClose(hFind);
	}
}
