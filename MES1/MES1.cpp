#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <future>
#include <windows.h>

//Control panel
#define N_POINTS 2
#define FILE_NAME "Test1_4_4.txt"
//#define FILE_NAME "Test2_4_4_MixGrid.txt"
//#define FILE_NAME "Test3_31_31_kwadrat.txt"
#define CALCULATE_HBC true
#define PRINT_HBC false
#define PRINT_P_LOCAL false
#define PRINT_P_GLOBAL false
#define PRINT_H_LOCAL false
#define PRINT_GLOBAL_H false
#define PRINT_GLOBAL_C false
#define PRINT_DATA false
#define PRINT_LOCAL_C false
#define PRINT_SIMULATION false
#define PRINT_SIMULATION_SHORT true

using namespace std;

//Print vector NxM
void coutVector2D(vector<vector<double>> vec, int size) {
	for (auto& x : vec)
	{
		for (auto& y : x) {
			if (abs(y) < 0.000001) cout << 0;
			else cout << ceil(y * 100.0) / 100.0;
			for (int i = 0; i < size; i++)
				cout << "\t";
		}
		cout << "\n";
	}
	cout << "\n";
}

//Print vector 1xN
void coutVector1D(vector<double> vec, int size) {
	for (auto& x : vec)
	{
		if (abs(x) < 0.000001) cout << 0;
		else cout << ceil(x * 100.0) / 100.0;
		for (int i = 0; i < size; i++)
			cout << "\t";
		cout << "\n";
	}
	cout << "\n";
}

//Info about point
struct Node
{
	double x, y, t = 0;
	short int BC = 0;
};

struct Element
{
	int id[4];
};

struct Grid
{
	int NodesNumber;
	int ElementsNumber;
	std::vector <Node> nodes;
	std::vector <Element> elements;
};

struct GlobalData
{
	int SimulationTime; // Czas symulacji
	int SimulationStepTime; // Przyrost czasu
	int Conductivity; // Współczynnik przewodności cieplnej
	int Alfa;
	int Tot; // Tempreratura otoczenia
	int InitialTemp; // Temperatura otoczenia
	int Density; // Gęstość
	int SpecificHeat; // Ciepło właściwe
};

class parametry
{
public:
	Grid grid;
	GlobalData data;
};

void readFile(parametry& obj)
{
	std::ifstream infile(FILE_NAME);

	std::string keyword;
	while (infile >> keyword)
	{
		if (keyword == "SimulationTime")
			infile >> obj.data.SimulationTime;
		else if (keyword == "SimulationStepTime")
			infile >> obj.data.SimulationStepTime;
		else if (keyword == "Conductivity")
			infile >> obj.data.Conductivity;
		else if (keyword == "Alfa")
			infile >> obj.data.Alfa;
		else if (keyword == "Tot")
			infile >> obj.data.Tot;
		else if (keyword == "InitialTemp")
			infile >> obj.data.InitialTemp;
		else if (keyword == "Density")
			infile >> obj.data.Density;
		else if (keyword == "SpecificHeat")
			infile >> obj.data.SpecificHeat;
		else if (keyword == "Nodes")
		{
			std::string dumb;
			infile >> dumb;
			infile >> obj.grid.NodesNumber;
		}
		else if (keyword == "Elements")
		{
			std::string dumb;
			infile >> dumb;
			infile >> obj.grid.ElementsNumber;
		}
		else if (keyword == "*Node")
		{
			for (int i = 0; i < obj.grid.NodesNumber; i++)
			{
				std::string word;
				Node temp;
				infile >> word;
				infile >> word;
				temp.x = stod(word.substr(0, word.length() - 1));
				infile >> word;
				temp.y = stod(word);
				obj.grid.nodes.push_back(temp);
			}
		}
		else if (keyword == "*BC")
		{
			std::string word;
			while (infile >> word)
			{
				if (word[word.length()] == ',') word = word.substr(0, word.length() - 1);
				obj.grid.nodes[stoi(word) - 1].BC = 1;
			}
		}
		else if (keyword == "*Element,")
		{
			std::string word;
			infile >> word;
			for (int i = 0; i < obj.grid.ElementsNumber; i++)
			{
				Element temp;
				infile >> word;
				for (int j = 0; j < 3; j++)
				{
					infile >> word;
					temp.id[j] = stoi(word.substr(0, word.length() - 1));
				}
				infile >> word;
				temp.id[3] = stoi(word);
				obj.grid.elements.push_back(temp);
			}
		}
		else
		{
			std::cout << "Niepoprawny parametr: '" << keyword << "', pomijam." << std::endl;
		}
	}
}
//Divides dobules
double divf(double a, double b) { return (double)a / (double)b; }

struct e2v {
	vector<vector<double>> vectDnDksi;
	vector<vector<double>> vectDnDeta;
	vector<vector<double>> weights = {
		{2},
		{ 1,1 },
		{divf(5,9), divf(8,9), divf(5,9)},
		{divf(18 + sqrt(30),36), divf(18 - sqrt(30),36), divf(18 - sqrt(30),36), divf(18 + sqrt(30),36)}
	};
	vector<vector<double>> points = {
		/*n=1*/{0},
		/*n=2*/{(-1 / sqrt(3)), (1 / sqrt(3))},
		/*n=3*/{(-sqrt(divf(3, 5))), 0, sqrt(divf(3, 5))},
		/*n=4*/{-sqrt(3.0f / 7.0f + 2.0f / 7.0f * sqrt(6.0f / 5.0f)),-sqrt(3.0f / 7.0f - 2.0f / 7.0f * sqrt(6.0f / 5.0f)),sqrt(3.0f / 7.0f - 2.0f / 7.0f * sqrt(6.0f / 5.0f)),sqrt(3.0f / 7.0f + 2.0f / 7.0f * sqrt(6.0f / 5.0f))}
	};
	double Dn1Deta(double eta) {
		return -0.25 * (1 - eta);
	}
	double Dn2Deta(double eta) {
		return -0.25 * (1 + eta);
	}
	double Dn3Deta(double eta) {
		return 0.25 * (1 + eta);
	}
	double Dn4Deta(double eta) {
		return 0.25 * (1 - eta);
	}
	double Dn1Dksi(double ksi) {
		return -0.25 * (1 - ksi);
	}
	double Dn2Dksi(double ksi) {
		return 0.25 * (1 - ksi);
	}
	double Dn3Dksi(double ksi) {
		return 0.25 * (1 + ksi);
	}
	double Dn4Dksi(double ksi) {
		return -0.25 * (1 + ksi);
	}

	e2v(int option) {
		for (int i = 0; i < option; i++) {
			for (int j = 0; j < option; j++)
			{
				vectDnDksi.push_back(vector<double> {Dn1Dksi(points[option - 1][i]), Dn2Dksi(points[option - 1][i]), Dn3Dksi(points[option - 1][i]), Dn4Dksi(points[option - 1][i])});
				vectDnDeta.push_back(vector<double> {Dn1Deta(points[option - 1][j]), Dn2Deta(points[option - 1][j]), Dn3Deta(points[option - 1][j]), Dn4Deta(points[option - 1][j])});
			}
		}
	}
};

class calculate_mes {
private:
	vector<vector<double>> global_H;
	vector<vector<double>> global_C;
	vector<double> global_P;
	vector<double> weigths_line; //All weigths combinations in 1D array
	e2v element = e2v(N_POINTS);
	vector<double> detJ_H; //Changes in every agregation iteration in local_H, needed to calculate local_C
	vector<vector<double>> walls_cords{ {0,1}, {1,2}, {2,3}, {3,0} }; //Stepping by walls in local_P and HBC
	vector<vector<double>> cords_ksi; //Wall points needed to calculate local_P and HBC
	vector<vector<double>> cords_eta; //As above

public:
	parametry obj;

	int nPoints;
	int nPoints2;

	calculate_mes() {
		nPoints = N_POINTS;
		nPoints2 = nPoints * nPoints;
		readFile(obj);
		global_H = vector<vector<double>>(obj.grid.NodesNumber, vector<double>(obj.grid.NodesNumber, 0));
		global_C = vector<vector<double>>(obj.grid.NodesNumber, vector<double>(obj.grid.NodesNumber, 0));
		global_P = vector<double>(obj.grid.NodesNumber, 0);
		//Varaibles to read from file in future
		//Now static to develop
		//int obj.data.Conductivity = 25; // TO CHANGE
		//int obj.data.Alfa = 300; // TO CHANGE
		//int Tot = 1200; //TEMP
		if (PRINT_DATA) {
			cout << "kT: " << obj.data.Conductivity << endl;
			cout << "Alpha: " << obj.data.Alfa << endl;
			cout << "Tot: " << obj.data.Tot << endl;
			cout << "Density: " << obj.data.Density << endl;
			cout << "Specific Heat: " << obj.data.SpecificHeat << endl;
			cout << "Simulation Time: " << obj.data.SimulationTime << endl;
			cout << "Simulation Step Time: " << obj.data.SimulationStepTime << endl;
			cout << "Initial Temp: " << obj.data.InitialTemp << endl;
		}
		for (int i = 0; i < nPoints; i++) {
			vector<double> tempV;
			for (int j = 0; j < nPoints; j++) {
				weigths_line.push_back(element.weights[nPoints - 1][i] * element.weights[nPoints - 1][j]);
			}
		}
		for (int i = 0; i < 4; i++) {
			vector<double> temp;
			for (int j = 0; j < nPoints; j++) {
				temp.push_back(element.points[nPoints - 1][j]);
			}
			cords_ksi.push_back(temp);
			cords_eta.push_back(temp);
		}
		for (int j = 0; j < nPoints; j++) {
			cords_eta[0][j] = -1;
			cords_ksi[1][j] = 1;
			cords_eta[2][j] = 1;
			cords_ksi[3][j] = -1;
		}
	}

	vector<double> ElimGauss(std::vector<std::vector<double>> H, std::vector<double> P)
	{ //Gauss elimination method
		int size = obj.grid.NodesNumber;
		vector<double> result(obj.grid.NodesNumber, 0);
		for (int i = 0; i < size - 1; i++)
		{
			for (int j = i; j < size - 1; j++)
			{
				double temp = H[j + 1][i];
				for (int k = i; k < size; k++)
				{
					H[j + 1][k] = H[j + 1][k] - (temp / H[i][i]) * H[i][k];
				}
				P[j + 1] = P[j + 1] - (temp / H[i][i]) * P[i];
			}
		}
		//Wynik
		for (int i = 0; i < size; i++)
		{
			result[i] = P[i];
		}
		result[size - 1] = P[size - 1] / H[size - 1][size - 1];
		for (int i = size - 2; i >= 0; i--)
		{
			int temp = 0;
			for (int j = size - 1; j > i; j--)
			{
				temp = j;
				result[i] = result[i] - H[i][j] * result[j];
			}
			result[i] = result[i] / H[i][temp - 1];
		}
		return result;
	}

	double get_N(int number, double ksi, double eta) {
		if (number == 1)
			return 0.25 * (1 - ksi) * (1 - eta);
		else if (number == 2)
			return 0.25 * (1 + ksi) * (1 - eta);
		else if (number == 3)
			return 0.25 * (1 + ksi) * (1 + eta);
		else if (number == 4)
			return 0.25 * (1 - ksi) * (1 + eta);

		cout << "ERROR, getting unknown N" << endl; //error case, common mixtake in develop

		return 0.0;
	}

	vector<double> get_nodes_x(int numberOfElement) {
		vector<double> result;
		for (int i = 0; i < 4; i++) {
			result.push_back(obj.grid.nodes[obj.grid.elements[numberOfElement].id[i] - 1].x);
		}
		return result;
	}

	vector<double> get_nodes_y(int numberOfElement) {
		vector<double> result;
		for (int i = 0; i < 4; i++) {
			result.push_back(obj.grid.nodes[obj.grid.elements[numberOfElement].id[i] - 1].y);
		}
		return result;
	}

	vector<int> get_nodes_bc(int numberOfElement) {
		vector<int> result;
		for (int i = 0; i < 4; i++) {
			result.push_back(obj.grid.nodes[obj.grid.elements[numberOfElement].id[i] - 1].BC);
		}
		return result;
	}

	vector<vector<double>> calculate_local_H(int numberOfElement) {
		vector<double> dxDksi;
		vector<double> dyDksi;
		vector<double> dxDeta;
		vector<double> dyDeta;
		vector<double> overDetJ;
		vector<vector<double>> dnDx;
		vector<vector<double>> dnDy;
		vector<vector<vector<double>>> matrixH;
		vector<vector<double>> matrixHTotal;
		vector<double> x = get_nodes_x(numberOfElement);
		vector<double> y = get_nodes_y(numberOfElement);

		//Calculate (dx/dksi,dy/dksi,dx/deta,dy/deta)
		for (int i = 0; i < nPoints2; i++) {
			double temp1 = 0;
			double temp2 = 0;
			double temp3 = 0;
			double temp4 = 0;
			for (int j = 0; j < 4; j++) {
				temp1 += element.vectDnDksi[i][j] * x[j];
				temp2 += element.vectDnDksi[i][j] * y[j];
				temp3 += element.vectDnDeta[i][j] * x[j];
				temp4 += element.vectDnDeta[i][j] * y[j];
			}
			dxDksi.push_back(temp1);
			dyDksi.push_back(temp2);
			dxDeta.push_back(temp3);
			dyDeta.push_back(temp4);
		}

		//calculate detJ
		detJ_H.clear();
		for (int i = 0; i < nPoints2; i++) {
			double temp = 0;
			temp = dxDksi[i] * dyDeta[i] - dxDeta[i] * dyDksi[i];
			detJ_H.push_back(temp);
			overDetJ.push_back(1 / temp);
		}

		//elemets * 1/detJ
		for (int i = 0; i < nPoints2; i++) {
			dxDksi[i] *= overDetJ[i];
			dyDksi[i] *= overDetJ[i];
			dxDeta[i] *= overDetJ[i];
			dyDeta[i] *= overDetJ[i];
		}
		for (int i = 0; i < nPoints2; i++) {
			vector<double> tempV;
			vector<double> tempV2;
			for (int j = 0; j < 4; j++) {
				tempV.push_back(dyDksi[i] * element.vectDnDeta[i][j] - dyDeta[i] * element.vectDnDksi[i][j]);
				tempV2.push_back(-dxDksi[i] * element.vectDnDeta[i][j] + dxDeta[i] * element.vectDnDksi[i][j]);
			}
			dnDx.push_back(tempV);
			dnDy.push_back(tempV2);
		}

		for (int i = 0; i < nPoints2; i++) {
			vector<vector<double>> tempV2;
			for (int j = 0; j < 4; j++) {
				vector<double> tempV;
				for (int k = 0; k < 4; k++) {
					tempV.push_back((dnDx[i][k] * dnDx[i][j] + dnDy[i][k] * dnDy[i][j]) * obj.data.Conductivity * detJ_H[i]);
				}
				tempV2.push_back(tempV);
			}
			matrixH.push_back(tempV2);
		}

		// Calculating HBC
		vector<vector<double>> final_hbc(4, vector<double>(4));
		if (CALCULATE_HBC)
			final_hbc = calculate_hbc(numberOfElement);

		for (int i = 0; i < 4; i++) {
			vector<double> tempV;
			for (int j = 0; j < 4; j++) {
				double temp = 0;
				temp += final_hbc[i][j];
				for (int k = 0; k < nPoints2; k++) {
					temp += matrixH[k][i][j] * weigths_line[k];
				}
				tempV.push_back(temp);
			}
			matrixHTotal.push_back(tempV);
		}

		if (PRINT_H_LOCAL) {
			cout << "LOCAL_H FOR " << numberOfElement << " ELEMENT\n";
			coutVector2D(matrixHTotal, 1);
		}

		return matrixHTotal;
	}

	vector<vector<double>> calculate_hbc(int numberOfElement) {
		vector<double> wagi = element.weights[nPoints - 1];
		vector<int> nodes_bc = get_nodes_bc(numberOfElement);
		//vector<int> nodes_bc = {1,1,0,0}; //change to auto select  //TEMP

		vector<double> x = get_nodes_x(numberOfElement);
		vector<double> y = get_nodes_y(numberOfElement);
		//vector<double> x = { 0, 0.025, 0.025, 0 }; //TEMP
		//vector<double> y = { 0, 0, 0.025, 0.025 }; //TEMP

		vector<int> walls_bc(4, 0);
		vector<vector<double>> final_hbc(4, vector<double>(4));
		vector<double> detJ;

		for (int i = 0; i < 4; i++)
			if (nodes_bc[walls_cords[i][0]] == 1 && nodes_bc[walls_cords[i][1]] == 1) walls_bc[i] = 1;

		//oblicznie detJ
		for (int i = 0; i < 4; i++) {
			detJ.push_back(sqrt(pow(x[walls_cords[i][0]] - x[walls_cords[i][1]], 2) + pow(y[walls_cords[i][0]] - y[walls_cords[i][1]], 2)) / 2);
		}

		for (int i = 0; i < 4; i++) { //walls
			if (walls_bc[i] == 0) continue;

			for (int j = 0; j < nPoints; j++) { //points on wall
				vector<double> part_pc;
				for (int k = 0; k < 4; k++) part_pc.push_back(get_N(k + 1, cords_ksi[i][j], cords_eta[i][j])); //Get N1-4
				//mnożenie przez wage i mnozenie transponowanej macierzy potem zsumowac party i pomnozyc przez detJ
				vector<vector<double>> part2_pc; //hbc for wall (part2_pc)
				for (int o = 0; o < 4; o++) { //(part2_pc.column)
					vector<double> tempV;
					for (int p = 0; p < 4; p++) {//(part2_pc.row)
						tempV.push_back(wagi[j] * obj.data.Alfa * part_pc[p] * part_pc[o]);
					}
					part2_pc.push_back(tempV);
				}
				//suma
				for (int o = 0; o < 4; o++) { //sum all walls
					for (int p = 0; p < 4; p++) {
						final_hbc[o][p] += part2_pc[o][p] * detJ[i];
					}
				}
			}
		}

		if (PRINT_HBC) {
			cout << "LOCAL_HBC FOR " << numberOfElement << " ELEMENT\n";
			coutVector2D(final_hbc, 1);
		}
		return final_hbc;
	}

	vector<double> calculate_local_P(int numberOfElement) {
		vector<int> nodes_bc = get_nodes_bc(numberOfElement);
		//vector<int> nodes_bc = { 1,0,0,1 }; //change to auto select  //TEMP

		vector<double> x = get_nodes_x(numberOfElement);
		vector<double> y = get_nodes_y(numberOfElement);
		//vector<double> x = { 0, 0.025, 0.025, 0 }; //TEMP
		//vector<double> y = { 0, 0, 0.025, 0.025 }; //TEMP

		vector<double> wagi = element.weights[nPoints - 1];
		vector<int> walls_bc(4, 0);
		vector<double> final_P(4, 0);
		vector<double> detJ;

		for (int i = 0; i < 4; i++)
			if (nodes_bc[walls_cords[i][0]] == 1 && nodes_bc[walls_cords[i][1]] == 1) walls_bc[i] = 1;

		//oblicznie detJ
		for (int i = 0; i < 4; i++) {
			detJ.push_back(sqrt(pow(x[walls_cords[i][0]] - x[walls_cords[i][1]], 2) + pow(y[walls_cords[i][0]] - y[walls_cords[i][1]], 2)) / 2);
		}

		//here calculate vect P
		for (int i = 0; i < 4; i++) { //walls
			if (walls_bc[i] == 0) continue;

			for (int j = 0; j < nPoints; j++) { //points on wall
				vector<double> part_pc;
				for (int k = 0; k < 4; k++) part_pc.push_back(get_N(k + 1, cords_ksi[i][j], cords_eta[i][j])); //Get N1-4
				vector<double> part2_pc; //vectP for wall (part2_pc)
				for (int o = 0; o < 4; o++) { //(part2_pc.column)
					double tempV;
					tempV = wagi[j] * obj.data.Alfa * part_pc[o] * obj.data.Tot;
					part2_pc.push_back(tempV);
				}
				//suma
				for (int o = 0; o < 4; o++) { //sum all walls
					final_P[o] += part2_pc[o] * detJ[i];
				}
			}
		}
		if (PRINT_P_LOCAL) {
			cout << "LOCAL_P FOR " << numberOfElement << " ELEMENT\n";
			coutVector1D(final_P, 1);
		}
		return final_P;
	}

	vector<vector<double>> calculate_local_C(int numberOfElement) {
		vector<vector<double>> matrixCTotal;
		vector<vector<double>> matrix_N;

		for (int i = 0; i < nPoints; i++) { //Calculate N1-4 in every combination ksi eta
			for (int j = 0; j < nPoints; j++) {
				vector<double> temp;
				for (int k = 0; k < 4; k++)
					temp.push_back(get_N(k + 1, element.points[nPoints - 1][i], element.points[nPoints - 1][j]));
				matrix_N.push_back(temp);
			}
		}

		for (int i = 0; i < 4; i++) { //Calculate matrix C
			vector<double> tempV;
			for (int j = 0; j < 4; j++) {
				double temp = 0;
				for (int k = 0; k < nPoints2; k++) {
					temp += matrix_N[k][i] * matrix_N[k][j] * detJ_H[k] * weigths_line[k] * obj.data.Density * obj.data.SpecificHeat;
				}
				tempV.push_back(temp);
			}
			matrixCTotal.push_back(tempV);
		}

		if (PRINT_LOCAL_C) {
			cout << "LOCAL_C FOR " << numberOfElement << " ELEMENT\n";
			coutVector2D(matrixCTotal, 1);
		}
		return matrixCTotal;
	}

	void agregation() {
		for (int k = 0; k < obj.grid.ElementsNumber; k++) {
			vector<vector<double>> local_H = calculate_local_H(k); //MUST BE FIRST! (FILLS DETJ_H AND OTHERS)
			vector<double> local_P = calculate_local_P(k);
			vector<vector<double>> local_C = calculate_local_C(k);
			for (int i = 0; i < 4; i++) {
				global_P[obj.grid.elements[k].id[i] - 1] += local_P[i];
				for (int j = 0; j < 4; j++) {
					global_H[obj.grid.elements[k].id[i] - 1][obj.grid.elements[k].id[j] - 1] += local_H[i][j];
					global_C[obj.grid.elements[k].id[i] - 1][obj.grid.elements[k].id[j] - 1] += local_C[i][j];
				}
			}
		}

		if (PRINT_P_GLOBAL) {
			cout << "GLOBAL_P:\n";
			coutVector1D(global_P, 1);
		}
		if (PRINT_GLOBAL_C) {
			cout << "GLOBAL_C:\n";
			coutVector2D(global_C, 1);
		}
		if (PRINT_GLOBAL_H) {
			cout << "GLOBAL_MATRIX:\n";
			coutVector2D(global_H, 1);
		}
	}

	void time_loop() {
		agregation();
		vector<vector<double>>H_C_dT(obj.grid.NodesNumber, vector<double>(obj.grid.NodesNumber, 0));
		vector<double>P_nowe(obj.grid.NodesNumber, 0);
		vector<double>t0(obj.grid.NodesNumber, obj.data.InitialTemp);
		for (int i = 0; i < obj.grid.NodesNumber; i++) {
			for (int j = 0; j < obj.grid.NodesNumber; j++) {
				H_C_dT[i][j] += global_H[i][j] + (global_C[i][j] / obj.data.SimulationStepTime);
			}
		}

		for (int i = 0; i < obj.data.SimulationTime; i += obj.data.SimulationStepTime) {
			for (int i = 0; i < obj.grid.NodesNumber; i++) {
				for (int j = 0; j < obj.grid.NodesNumber; j++) {
					P_nowe[i] += (global_C[i][j] / obj.data.SimulationStepTime) * t0[j];
				}
				P_nowe[i] += global_P[i];
			}
			t0 = ElimGauss(H_C_dT, P_nowe);
			P_nowe = vector<double>(obj.grid.NodesNumber, 0);

			if (PRINT_SIMULATION) {
				cout << "TIME: " << i + obj.data.SimulationStepTime << endl;
				coutVector1D(t0, 1);
			}
			if (PRINT_SIMULATION_SHORT) {
				cout << "TIME: " << i + obj.data.SimulationStepTime << endl;
				cout << "MIN: " << *min_element(t0.begin(), t0.end()) << "\nMAX: " << *max_element(t0.begin(), t0.end()) << "\n\n";
			}
		}
	}
};

//not using now
double fun1d(double x) {
	//return 2 * pow(x, 2) + 3 * x - 8;
	return pow(x, 2) + 3 * x + 2;
}
//not using now
double fun2d(double ksi, double eta) {
	return -5 * pow(ksi, 2) * eta + 2 * ksi * pow(eta, 2) + 10;
}
//not using now
double integral(int k, std::vector<double> pc, std::vector<double> a, int dim)
{
	double result = 0;
	if (dim == 1)
	{
		for (int i = 0; i < k; i++)
			result += fun1d(pc[i]) * a[i];
	}
	else if (dim == 2)
	{
		for (int i = 0; i < k; i++)
			for (int j = 0; j < k; j++)
				result += fun2d(pc[i], pc[j]) * a[i] * a[j];
	}
	return result;
}

int main()
{
	calculate_mes mes;
	mes.time_loop();

	//cout << "Hello!";
	//std::cout << integral(2, std::vector<double>{ -divf(1, sqrtf(3)), divf(1, sqrtf(3)) }, { 1,1 }, 1);
	//std::cout << integral(3, std::vector<double>{ -sqrtf(.6), 0, sqrtf(.6) }, { divf(5,9),divf(8,9),divf(5,9) }, 1);
	//std::cout << integral(3, std::vector<double>{ -sqrtf(.6), 0, sqrtf(.6) }, { divf(5,9),divf(8,9),divf(5,9) }, 2);
}

////////////////////////////////////////////////////////////////////////
//
// When my code doesn't compile
// I recall what Winston Churchill said
//
// "Underidoderidoderiododeridoo" ~ Winston Churchill
//
// And that makes my code work well :)
//
////////////////////////////////////////////////////////////////////////