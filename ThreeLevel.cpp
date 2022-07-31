#include "stdafx.h"
#include "ThreeLevel.h"

#include "Structure.h"
#include "Fields.h"
#include "Constants.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <time.h> 
using namespace std;

ThreeLevel::ThreeLevel(Structure &sS, char *fileName)
{
	cout << "Reading four level gain material properties from file \"" << fileName << "\" ";
	ifstream InputFile;
	string sline;
	istringstream iss;
	InputFile.open(fileName);
	getline(InputFile, sline);
	iss.str(sline);
	iss >> distri;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> NumSources;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> i1;
	iss >> i2;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> j1;
	iss >> j2;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> lambda_a;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> lambda_b;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> gamma_a;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> gamma_b;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> t20;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> t21;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> t10;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> N_0;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> N_1;
	cout << "..." << endl << "Constructing gain material";

	omega_a = 2.0*pi*C / lambda_a;
	omega_b = 2.0*pi*C / lambda_b;

	CPa1 = 2 / ((sS.dt*sS.dt + omega_a*omega_a) * (1 / (sS.dt*sS.dt) + gamma_a / (2 * sS.dt)));
	CPa2 = (2 / (sS.dt*sS.dt) + 1/(2*sS.dt)) / (1 / (sS.dt*sS.dt) + gamma_a / (2 * sS.dt));
	CPa3 = (qe*qe / me) / (1 / (sS.dt*sS.dt) + gamma_a / (2 * sS.dt)) / N_0;

	CPb1 = 2 / ((sS.dt*sS.dt + omega_b*omega_b) * (1 / (sS.dt*sS.dt) + gamma_b / (2 * sS.dt)));
	CPb2 = (2 / (sS.dt*sS.dt) + 1 / (2 * sS.dt)) / (1 / (sS.dt*sS.dt) + gamma_b / (2 * sS.dt));
	CPb3 = (qe*qe / me) / (1 / (sS.dt*sS.dt) + gamma_b / (2 * sS.dt)) / N_0;

	CN2N2 = (1.0 / sS.dt - 1.0 / (2.0 * t20) - 1.0 / (2.0 * t21)) / (1.0 / sS.dt + 1.0 / (2.0 * t20) + 1.0 / (2.0 * t21));
	CN2E = (1.0 / (2.0 * sS.dt*hPlanck*omega_a / 2.0 / pi)) / (1.0 / sS.dt + 1.0 / (2.0 * t20) + 1.0 / (2.0 * t21));
	CN1N1 = (1.0 / sS.dt - 1.0 / (2.0 * t10)) / (1.0 / sS.dt + 1.0 / (2.0 * t10));
	CN1N2 = (1.0 / (2.0 * t21)) / (1.0 / sS.dt + 1.0 / (2.0 * t10));
	CN1E = (1.0 / (2.0 * sS.dt*hPlanck*omega_b/2.0/pi)) / (1.0 / sS.dt + 1.0 / (2.0 * t10));
	CN0N1 = sS.dt / (2.0 * t10);
	CN0N2 = sS.dt / (2.0 * t20);
	CN0Ea = 1.0 / (2.0*hPlanck*omega_a / 2.0 / pi);
	CN0Eb = 1.0 / (2.0*hPlanck*omega_b / 2.0 / pi);

	//Locating
	idif = i2 - i1;
	jdif = j2 - j1;
	if (distri == "random"){
		srand(time(NULL));
		location = new int*[NumSources];
		for (i = 0; i < NumSources; ++i){
			location[i] = new int[2];
			location[i][0] = rand() % idif + i1;
			location[i][1] = rand() % jdif + j1;
		}
	}
	else{
		ii = 0;
		for (i = 0; i < sS.Nx; ++i){
			for (j = 0; j < sS.Ny; ++j){
				if (sS.ID[i][j] == 5){
					ii = ii + 1;
				}
			}
		}
		NumSources = ii + 1;//
		location = new int*[NumSources];
		ii = 0;
		for (i = 0; i < sS.Nx; ++i){
			for (j = 0; j < sS.Ny; ++j){
				if (sS.ID[i][j] == 5){
					location[ii] = new int[2];
					location[ii][0] = i;
					location[ii][1] = j;
					ii = ii + 1;
				}
			}
		}
		location[ii] = new int[2];
		location[ii][0] = i;
		location[ii][1] = j;
	}
	Pax = new double[NumSources];
	Pay = new double[NumSources];
	//Paz = new double[NumSources];
	Pbx = new double[NumSources];
	Pby = new double[NumSources];
	//Pbz = new double[NumSources];

	Pax_1 = new double[NumSources];
	Pay_1 = new double[NumSources];
	//Paz_1 = new double[NumSources];
	Pbx_1 = new double[NumSources];
	Pby_1 = new double[NumSources];
	//Pbz_1 = new double[NumSources];

	Pax_2 = new double[NumSources];
	Pay_2 = new double[NumSources];
	//Paz_2 = new double[NumSources];
	Pbx_2 = new double[NumSources];
	Pby_2 = new double[NumSources];
	//Pbz_2 = new double[NumSources];

	Ex_1 = new double[NumSources];
	Ey_1 = new double[NumSources];
	//Ez_1 = new double[NumSources];

	N0 = new double[NumSources];
	N1 = new double[NumSources];
	N2 = new double[NumSources];
	N1_1 = new double[NumSources];
	N2_1 = new double[NumSources];

	for (i = 0; i < NumSources; ++i){
		Pax[i] = 0.0;
		Pay[i] = 0.0;
		//Paz[i] = 0.0;
		Pbx[i] = 0.0;
		Pby[i] = 0.0;
		//Pbz[i] = 0.0;
		Pax_1[i] = 0.0;
		Pay_1[i] = 0.0;
		//Paz_1[i] = 0.0;
		Pbx_1[i] = 0.0;
		Pby_1[i] = 0.0;
		//Pbz_1[i] = 0.0;
		Pax_2[i] = 0.0;
		Pay_2[i] = 0.0;
		//Paz_2[i] = 0.0;
		Pbx_2[i] = 0.0;
		Pby_2[i] = 0.0;
		//Pbz_2[i] = 0.0;
		N0[i] = N_0; //0.0;
		N1[i] = 0.0; // N_1;
		N2[i] = 0.0;
		N1_1[i] = 0.0;// N_1;
		N2_1[i] = 0.0;
		Ex_1[i] = 0.0;
		Ey_1[i] = 0.0;
		//Ez_1[i] = 0.0;
	}
	CEP = 1.0 / sS.dt / eps0;
	cout << "..." << endl;
}


ThreeLevel::~ThreeLevel()
{
}


void ThreeLevel::UpdatePolarizationDensities(Field &fF)
{
	for (i = 0; i < NumSources; ++i)
	{
		Pax_2[i] = Pax_1[i];
		Pax_1[i] = Pax[i];
		Pax[i] = CPa1*Pax[i] + CPa2*Pax_2[i] + fF.Ex[location[i][0]][location[i][1]] * CPa3 * (N2[i] - N0[i]);
		Pay_2[i] = Pay_1[i];
		Pay_1[i] = Pay[i];
		Pay[i] = CPa1*Pay[i] + CPa2*Pay_2[i] + fF.Ey[location[i][0]][location[i][1]] * CPa3 * (N2[i] - N0[i]);
		//		Paz_2[i] = Paz_1[i];
		//		Paz_1[i] = Paz[i];
		//		Paz[i] = CPa1*Paz[i] + CPa2*Paz_2[i] + CPa3*(N2[i] - N0[i])*fF.Ez[location[i][0]][location[i][1]];

		Pbx_2[i] = Pbx_1[i];
		Pbx_1[i] = Pbx[i];
		Pbx[i] = CPb1*Pbx[i] + CPb2*Pbx_2[i] + fF.Ex[location[i][0]][location[i][1]] * CPb3 * (N1[i] - N0[i]);
		Pby_2[i] = Pby_1[i];
		Pby_1[i] = Pby[i];
		Pby[i] = CPb1*Pby[i] + CPb2*Pby_2[i] + fF.Ey[location[i][0]][location[i][1]] * CPb3*(N1[i] - N0[i]);
		//		Pbz_2[i] = Pbz_1[i];
		//		Pbz_1[i] = Pbz[i];
		//		Pbz[i] = CPb1*Pbz[i] + CPb2*Pbz_2[i] + CPb3*(N1[i] - N0[i])*fF.Ez[location[i][0]][location[i][1]];

		Ex_1[i] = fF.Ex[location[i][0]][location[i][1]];
		Ey_1[i] = fF.Ey[location[i][0]][location[i][1]];
		//		Ez_1[i] = fF.Ez[location[i][0]][location[i][1]];
	}
}

void ThreeLevel::UpdateElectricField(Field &fF)
{
	for (i = 0; i < NumSources; ++i){
		fF.Ex[location[i][0]][location[i][1]] = fF.Ex[location[i][0]][location[i][1]] - CEP*(Pax[i] - Pax_1[i] + Pbx[i] - Pbx_1[i]);
		fF.Ey[location[i][0]][location[i][1]] = fF.Ey[location[i][0]][location[i][1]] - CEP*(Pay[i] - Pay_1[i] + Pby[i] - Pby_1[i]);
		//		fF.Ez[location[i][0]][location[i][1]] = fF.Ez[location[i][0]][location[i][1]] - CEP*(Paz[i] - Paz_1[i] + Pbz[i] - Pbz_1[i]);
	}
}

void ThreeLevel::UpdateStatesPopulations(Field &fF)
{
	for (i = 0; i < NumSources; ++i){
		N2_1[i] = N2[i];
		N2[i] = CN2N2 * N2[i] + CN2E * (fF.Ex[location[i][0]][location[i][1]] + Ex_1[i]) * (Pax[i] - Pax_1[i]) + CN2E * (fF.Ey[location[i][0]][location[i][1]] + Ey_1[i]) * (Pay[i] - Pay_1[i]);
		N1_1[i] = N1[i];
		N1[i] = CN1N1 * N1[i] + CN1N2 * (N2[i] + N2_1[i]) + CN1E * (fF.Ex[location[i][0]][location[i][1]] + Ex_1[i]) * (Pbx[i] - Pbx_1[i]) + CN1E * (fF.Ey[location[i][0]][location[i][1]] + Ey_1[i]) * (Pby[i] - Pby_1[i]);
		N0[i] = N0[i] + CN0N2 * (N2[i] + N2_1[i]) + CN0N1 * (N1[i] + N1_1[i]) - CN0Ea * (fF.Ex[location[i][0]][location[i][1]] + Ex_1[i]) * (Pax[i] - Pax_1[i]) - CN0Ea * (fF.Ey[location[i][0]][location[i][1]] + Ey_1[i]) * (Pay[i] - Pay_1[i]) - CN0Ea * (fF.Ex[location[i][0]][location[i][1]] + Ex_1[i]) * (Pbx[i] - Pbx_1[i]) - CN0Eb * (fF.Ey[location[i][0]][location[i][1]] + Ey_1[i]) * (Pby[i] - Pby_1[i]);
	}
}