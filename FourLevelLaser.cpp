//    This file is part of the arhFDTD package
//    Copyright (C) <2022>  <AliReza Hashemi>

//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Affero General Public License as
//    published by the Free Software Foundation, either version 3 of the
//    License, or (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Affero General Public License for more details.

//    You should have received a copy of the GNU Affero General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#include "stdafx.h"
#include "FourLevelLaser.h"

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

FourLevelLaser::FourLevelLaser(Structure &sS, char *fileName)
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
	iss >> lambda_b;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> lambda_a;
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
	iss >> t30;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> t32;
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
	zeta_a = 6.0*pi*eps0*C*C*C / (omega_a*omega_a*t21);
	zeta_b = 6.0*pi*eps0*C*C*C / (omega_b*omega_b*t30);
	CPa1 = (2.0 - omega_a*omega_a*sS.dt*sS.dt) / (1.0 + gamma_a*sS.dt / 2.0);
	CPa2 = (gamma_a*sS.dt / 2.0 - 1.0) / (1.0 + gamma_a*sS.dt / 2.0);
	CPa3 = sS.dt*sS.dt * zeta_a / (1.0 + gamma_a*sS.dt / 2.0);
	CPb1 = (2.0 - omega_b*omega_b*sS.dt*sS.dt) / (1.0 + gamma_b*sS.dt / 2.0);
	CPb2 = (gamma_b*sS.dt / 2.0 - 1.0) / (1.0 + gamma_b*sS.dt / 2.0);
	CPb3 = sS.dt*sS.dt * zeta_b / (1.0 + gamma_b*sS.dt / 2.0);

	//Locating
	idif = i2 - i1;
	jdif = j2 - j1;
	if (distri == random_ditri){
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
	N3 = new double[NumSources];

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
		N0[i] = 1.0;// N_0;
		N1[i] = 1.0;// N_1;
		N2[i] = 0.0;
		N3[i] = 0.0;
		Ex_1[i] = 0.0;
		Ey_1[i] = 0.0;
		//Ez_1[i] = 0.0;
	}
	CEP = 2.0 * N_0 * sS.dx*sS.dy*sS.dx / eps0 / sS.dt;
	t30t32dt = t30*t32*sS.dt; t30t32_2 = 2.0*t30*t32; t30dt = t30*sS.dt; t32dt = t32*sS.dt; dt_2 = 2.0 / sS.dt; homega_bdt = hPlanck*omega_b*sS.dt / 2.0 / pi;
	t10t21dt = t10*t21*sS.dt; t10t21_2 = 2.0*t10*t21; t10dt = t10*sS.dt; t21dt = t21*sS.dt; t21_2 = t21 / 2.0; homega_adt = hPlanck*omega_a*sS.dt / 2.0 / pi;
	t32t21dt = t32*t21*sS.dt; t32t21_2 = 2.0*t32*t21; t32_2 = t32 / 2.0;
	t10t30dt = t10*t30*sS.dt; t10t30_2 = 2.0*t10*t30;

	dtt32 = sS.dt / t32; dtt30 = sS.dt / t30; homega_b_2 = 2 * hPlanck*omega_b;
	dtt21 = sS.dt / t21; dtt10 = sS.dt / t10; homega_a_2 = 2 * hPlanck*omega_a;


	cout << "..." << endl;
}


FourLevelLaser::~FourLevelLaser()
{
	cout << "Distructing gain material";
	delete[] Pax;
	Pax = 0;
	delete[] Pay;
	Pay = 0;
	delete[] Paz;
	Paz = 0;
	delete[] Pbx;
	Pbx = 0;
	delete[] Pby;
	Pby = 0;
	delete[] Pbz;
	Pbz = 0;
	delete[] Pax_1;
	Pax_1 = 0;
	delete[] Pay_1;
	Pay_1 = 0;
	delete[] Paz_1;
	Paz_1 = 0;
	delete[] Pbx_1;
	Pbx_1 = 0;
	delete[] Pby_1;
	Pby_1 = 0;
	delete[] Pbz_1;
	Pbz_1 = 0;
	delete[] Pax_2;
	Pax_2 = 0;
	delete[] Pay_2;
	Pay_2 = 0;
	delete[] Paz_2;
	Paz_2 = 0;
	delete[] Pbx_2;
	Pbx_2 = 0;
	delete[] Pby_2;
	Pby_2 = 0;
	delete[] Pbz_2;
	Pbz_2 = 0;
	delete[] N0;
	N0 = 0;
	delete[] N1;
	N1 = 0;
	delete[] N2;
	N2 = 0;
	delete[] N3;
	N3 = 0;
	for (int i = 0; i < NumSources; i++){
		delete[] location[i];
	}
	delete[] location;
	location = 0;
	cout << "..." << endl;
}
void FourLevelLaser::UpdatePolarizationDensities(Field &fF)
{
	for (i = 0; i < NumSources; ++i)
	{
		Pax_2[i] = Pax_1[i];
		Pax_1[i] = Pax[i];
		Pax[i] = CPa1*Pax[i] + CPa2*Pax_2[i] + fF.Ex[location[i][0]][location[i][1]] * CPa3 * (N2[i] - N1[i]);
		Pay_2[i] = Pay_1[i];
		Pay_1[i] = Pay[i];
		Pay[i] = CPa1*Pay[i] + CPa2*Pay_2[i] + fF.Ey[location[i][0]][location[i][1]] * CPa3 * (N2[i] - N1[i]);
//		Paz_2[i] = Paz_1[i];
//		Paz_1[i] = Paz[i];
//		Paz[i] = CPa1*Paz[i] + CPa2*Paz_2[i] + CPa3*(N2[i] - N1[i])*fF.Ez[location[i][0]][location[i][1]];

		Pbx_2[i] = Pbx_1[i];
		Pbx_1[i] = Pbx[i];
		Pbx[i] = CPb1*Pbx[i] + CPb2*Pbx_2[i] + fF.Ex[location[i][0]][location[i][1]] * CPb3 * (N3[i] - N0[i]);
		Pby_2[i] = Pby_1[i];
		Pby_1[i] = Pby[i];
		Pby[i] = CPb1*Pby[i] + CPb2*Pby_2[i] + fF.Ey[location[i][0]][location[i][1]] * CPb3*(N3[i] - N0[i]);
//		Pbz_2[i] = Pbz_1[i];
//		Pbz_1[i] = Pbz[i];
//		Pbz[i] = CPb1*Pbz[i] + CPb2*Pbz_2[i] + CPb3*(N3[i] - N0[i])*fF.Ez[location[i][0]][location[i][1]];

		Ex_1[i] = fF.Ex[location[i][0]][location[i][1]];
		Ey_1[i] = fF.Ey[location[i][0]][location[i][1]];
		//		Ez_1[i] = fF.Ez[location[i][0]][location[i][1]];
	}
}

void FourLevelLaser::UpdateElectricField(Field &fF)
{
	for (i = 0; i < NumSources; ++i){		
		fF.Ex[location[i][0]][location[i][1]] = fF.Ex[location[i][0]][location[i][1]] - CEP*(Pax[i] - Pax_1[i] + Pbx[i] - Pbx_1[i]);
		fF.Ey[location[i][0]][location[i][1]] = fF.Ey[location[i][0]][location[i][1]] - CEP*(Pay[i] - Pay_1[i] + Pby[i] - Pby_1[i]);
//		fF.Ez[location[i][0]][location[i][1]] = fF.Ez[location[i][0]][location[i][1]] - CEP*(Paz[i] - Paz_1[i] + Pbz[i] - Pbz_1[i]);
	}
	
}

void FourLevelLaser::UpdateStatesPopulations(Field &fF)
{
	for (i = 0; i < NumSources; ++i){
		N3[i] = t30t32dt / (t30t32_2 + t30dt*(1 - N2[i]) + t32dt*(1 - N0[i])) * (N3[i] * ((N0[i] - 1) / t30 + (N2[i] - 1) / t32 + dt_2) + (Ex_1[i] * (Pbx[i] - Pbx_2[i]) + Ey_1[i] * (Pby[i] - Pby_2[i])) / homega_bdt);
		N1[i] = t10t21dt / (t10t21_2 + t10dt*N2[i] + t21dt*(1 - N0[i])) * (N1[i] * ((N0[i] - 1) / t10 - N2[i] / t21 + dt_2) + N2[i] / t21_2 - (Ex_1[i] * (Pax[i] - Pax_2[i]) + Ey_1[i] * (Pay[i] - Pay_2[i])) / homega_adt);

		N2[i] = t32t21dt / (t32t21_2 + t21dt*N3[i] + t32dt*(1 - N1[i])) * (N2[i] * ((N1[i] - 1) / t21 - N3[i] / t32 + dt_2) + N3[i] / t32_2 + ((fF.Ex[location[i][0]][location[i][1]] + Ex_1[i]) * (Pax[i] - Pax_2[i]) + (fF.Ey[location[i][0]][location[i][1]] + Ey_1[i]) * (Pay[i] - Pay_1[i])) / homega_adt);
		//N0[i] = t10t30dt / (t10t30_2 + t10dt*N3[i] + t30dt*N1[i]) * (N0[i] * (-N1[i] / t10 - N3[i] / t30 + dt_2) + N1[i] / t10 + N3[i] / t30 - ((fF.Ex[location[i][0]][location[i][1]] + Ex_1[i]) * (Pbx[i] - Pbx_2[i]) + (fF.Ey[location[i][0]][location[i][1]] + Ey_1[i]) * (Pby[i] - Pby_1[i])) / homega_bdt);
		

		//N3[i] = N3[i] - N3[i] * (1 - N2[i])*dtt32 - N3[i] * (1 - N0[i])*dtt30 + fF.Ex[location[i][0]][location[i][1]] * (Pbx[i] - Pbx_2[i]) / homega_b_2 + fF.Ey[location[i][0]][location[i][1]] * (Pby[i] - Pby_2[i]) / homega_b_2;
		//N1[i] = N1[i] + N2[i] * (1 - N1[i])*dtt21 - N1[i] * (1 - N0[i])*dtt10 - fF.Ex[location[i][0]][location[i][1]] * (Pax[i] - Pax_2[i]) / homega_a_2 + fF.Ey[location[i][0]][location[i][1]] * (Pay[i] - Pay_2[i]) / homega_a_2;
		//N2[i] = N2[i] + N3[i] * (1 - N2[i])*dtt32 - N2[i] * (1 - N1[i])*dtt21 + (fF.Ex[location[i][0]][location[i][1]] + Ex_1[i]) * (Pax[i] - Pax_2[i]) / homega_a_2 / 2 + (fF.Ey[location[i][0]][location[i][1]] + Ey_1[i]) * (Pay[i] - Pay_2[i]) / homega_a_2 / 2;
		////N0[i] = N0[i] + N3[i] * (1 - N0[i])*dtt30 + N1[i] * (1 - N0[i])*dtt10 - (fF.Ex[location[i][0]][location[i][1]] + Ex_1[i]) * (Pbx[i] - Pbx_2[i]) / homega_b_2 / 2 + (fF.Ey[location[i][0]][location[i][1]] + Ey_1[i]) * (Pby[i] - Pby_2[i]) / homega_b_2 / 2;
		N0[i] = 2.0 - N1[i] - N2[i] - N3[i];
	}
}
