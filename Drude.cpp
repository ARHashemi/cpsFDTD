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
#include "Drude.h"

#include "Constants.h"
#include "Structure.h"
#include "Fields.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>
using namespace std;

Drude::Drude(Structure &sS, char *fileName)
{
	cout << "Constructing Drude Material";
	ifstream InputFile;
	InputFile.open(fileName);
	string sline, id;
	getline(InputFile, sline);
	istringstream iss(sline);
	iss >> epsInf;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> omegaP;
	omegaP = omegaP * 2.0 * pi;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> sigma;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> gamma;
	gamma = gamma * 2.0 * pi;
	
	beta = (omegaP*omegaP*eps0*sS.dt / 2) / (1 + gamma*sS.dt / 2);
	kappa = (1 - gamma*sS.dt / 2) / (1 + gamma*sS.dt / 2);
	C1Drude = (2 * eps0*epsInf - sS.dt*beta - sigma*sS.dt) / (2 * eps0*epsInf + sS.dt*beta + sigma*sS.dt);
	C2Drude = (2 * sS.dt) / (2 * eps0*epsInf + sS.dt*beta + sigma*sS.dt);
	C3Drude = ((1 + kappa)*sS.dt) / (2 * eps0*epsInf + sS.dt*beta + sigma*sS.dt);

	Jx = new double*[sS.Nx + 1];
	Jy = new double*[sS.Nx + 1];
	int i, j;
	for (i = 0; i <= sS.Nx; i++){
		Jx[i] = new double[sS.Ny + 1];
		Jy[i] = new double[sS.Ny + 1];
		for (j = 0; j <= sS.Ny; j++){
			Jx[i][j] = 0.0;
			Jy[i][j] = 0.0;			
		}
	}
	cout << " ..." << endl;
}


Drude::~Drude()
{
}

void Drude::EFieldUpdate(Field &fF, Structure &sS)
{	
	unsigned short id;
	int i, j;
	double temp;
#pragma omp parallel for default(shared) private(i,j,id)
	for (i = 1; i < sS.Nx; ++i){
		for (j = 1; j < sS.Ny; ++j){
			id = sS.ID[i][j];
			if (id == 2){
				temp = fF.Ex[i][j];
				fF.Ex[i][j] = C1Drude * fF.Ex[i][j] + C2Drude * (fF.Hz[i][j] - fF.Hz[i][j - 1]) * sS.cExy[j] - C3Drude * Jx[i][j];
				Jx[i][j] = kappa * Jx[i][j] + beta * (fF.Ex[i][j] + temp);
				temp = fF.Ey[i][j];
				fF.Ey[i][j] = C1Drude * fF.Ey[i][j] - C2Drude * (fF.Hz[i][j] - fF.Hz[i - 1][j]) * sS.cEyx[i] - C3Drude * Jy[i][j];
				Jy[i][j] = kappa * Jy[i][j] + beta * (fF.Ey[i][j] + temp);
			}
			else if (id == 4){
				fF.Ex[i][j] = 0.0;
				fF.Ey[i][j] = 0.0;
			}
			else{
				fF.Ex[i][j] = sS.aE[id] * fF.Ex[i][j] + sS.bE[id] * (fF.Hz[i][j] - fF.Hz[i][j - 1]) * sS.cExy[j];
				fF.Ey[i][j] = sS.aE[id] * fF.Ey[i][j] - sS.bE[id] * (fF.Hz[i][j] - fF.Hz[i - 1][j]) * sS.cEyx[i];
			}
		}
	}
	i = 0;
	for (j = 1; j<sS.Ny; ++j){
		fF.Ex[i][j] = sS.aE[0] * fF.Ex[i][j] + sS.bE[0] * (fF.Hz[i][j] - fF.Hz[i][j - 1]) * sS.cExy[j];
	}
	j = 0;
	for (i = 1; i<sS.Nx; ++i){
		fF.Ey[i][j] = sS.aE[0] * fF.Ey[i][j] - sS.bE[0] * (fF.Hz[i][j] - fF.Hz[i - 1][j]) * sS.cEyx[i];
	}
}
