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
#include "Drude_Lorentz.h"

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

Drude_Lorentz::Drude_Lorentz(Structure &sS, char *fileName)
{
	int p, i, j, ii;
	cout << "Constructing Drude Material"<< endl;
	ifstream InputFile;
	InputFile.open(fileName);
	string sline;
	getline(InputFile, sline);
	istringstream iss(sline);
	iss >> num_of_D_poles;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> num_of_L_poles;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> sigma;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> epsInf;
	omega_D = new double[num_of_D_poles];
	gamma_D = new double[num_of_D_poles];
	delta_eps_L = new double[num_of_L_poles];
	omega_L = new double[num_of_L_poles];
	delta_L = new double[num_of_L_poles];
	beta = new double[num_of_D_poles];
	kappa = new double[num_of_D_poles];
	alpha = new double[num_of_L_poles];
	xi = new double[num_of_L_poles];
	gamma_L = new double[num_of_L_poles];
	temp1=0;
	temp2=0;
	cout << "omega_D\t Gamma_D" << endl;
	for (p=0;p<num_of_D_poles;p++){
		getline(InputFile, sline);
		iss.clear();iss.str(sline);
		iss >> omega_D[p];
		omega_D[p] = omega_D[p] * 2.0 * pi;
		cout << omega_D[p]<<"\t";
		iss >> gamma_D[p];
		gamma_D[p] = gamma_D[p] * 2.0 * pi;
		cout << gamma_D[p]<<endl;
		kappa[p] = (1.0 - gamma_D[p]*sS.dt / 2.0) / (1.0 + gamma_D[p]*sS.dt / 2.0);
		beta[p] = (omega_D[p]*omega_D[p]*eps0*sS.dt/2.0) / (1.0 + gamma_D[p]*sS.dt / 2.0);
		temp1=temp1+beta[p];
	}
	cout << "omega_L\t Gamma_L\t Delta_epsilon" << endl;
	for (p=0;p<num_of_L_poles;p++){
		getline(InputFile, sline);
		iss.clear();iss.str(sline);
		iss >> omega_L[p];
		omega_L[p] = omega_L[p] * 2.0 * pi;
		cout << omega_L[p]<<"\t";
		iss >> delta_L[p];
		delta_L[p] = delta_L[p] * 2.0 * pi;
		cout << delta_L[p]<<"\t";
		iss >> delta_eps_L[p];		
		cout << delta_eps_L[p]<<endl;
		alpha[p] = (2.0 - omega_L[p]*omega_L[p]*sS.dt*sS.dt) / (1.0 + delta_L[p]*sS.dt);
		xi[p] = (delta_L[p]*sS.dt - 1.0) / (1.0 + delta_L[p]*sS.dt);
		gamma_L[p] = (eps0*delta_eps_L[p]*omega_L[p]*omega_L[p]*sS.dt*sS.dt) / (1.0 + delta_L[p]*sS.dt);
		temp2=temp2+gamma_L[p];
	}
	denominator=2.0*eps0*epsInf + sigma*sS.dt + sS.dt*temp1 + temp2/2.0;
	C1 = (2.0*eps0*epsInf + sigma*sS.dt + sS.dt*temp1)/denominator;
	C2 = temp2/2.0/denominator;
	C3 = 2.0*sS.dt/denominator;
	C4 = new double[num_of_D_poles];
	for (p=0;p<num_of_D_poles;p++){
		C4[p] = (1.0+kappa[p])*sS.dt/denominator;
	}
	C5 = new double[num_of_L_poles];
	C6 = new double[num_of_L_poles];
	for (p=0;p<num_of_L_poles;p++){
		C5[p] = (1.0+alpha[p])/denominator;
		C6[p] = xi[p]*sS.dt/denominator;
	}
	ii = 0;
	for (i = 0; i < sS.Nx; ++i){
			for (j = 0; j < sS.Ny; ++j){
				if (sS.ID[i][j] == 2){
					ii = ii + 1;
				}
			}
	}
	num_of_DLM_gridcells = ii + 1;
	DLM = new DL_material[num_of_DLM_gridcells];
	ii = 0;
	for (i = 0; i < sS.Nx; ++i){
		for (j = 0; j < sS.Ny; ++j){
			if (sS.ID[i][j] == 2){
				DLM[ii].ix = i;
				DLM[ii].jy = j;
				DLM[ii].JDx = new double[num_of_D_poles];
				DLM[ii].JDy = new double[num_of_D_poles];
				for (p=0;p<num_of_D_poles;p++){
					DLM[ii].JDx[p] = 0.0;
					DLM[ii].JDy[p] = 0.0;
				}
				DLM[ii].JLx = new double[num_of_L_poles];
				DLM[ii].JLy = new double[num_of_L_poles];
				DLM[ii].JLx_1 = new double[num_of_L_poles];
				DLM[ii].JLy_1 = new double[num_of_L_poles];
				for (p=0;p<num_of_L_poles;p++){
					DLM[ii].JLx[p] = 0.0;
					DLM[ii].JLy[p] = 0.0;
					DLM[ii].JLx_1[p] = 0.0;
					DLM[ii].JLy_1[p] = 0.0;
				}
				DLM[ii].Ex_1 = 0.0;
				DLM[ii].Ey_1 = 0.0;
				ii = ii + 1;
			}
		}
	}
	
}
Drude_Lorentz::~Drude_Lorentz()
{
}

void Drude_Lorentz::EFieldUpdate(Field &fF, Structure &sS)	
{	unsigned short id;
	int p, i, j, ii;
	double temp, temp0;	
	ii = 0;
//#pragma omp parallel for default(shared) private(i,j,id,p,temp)
	for (i = 1; i < sS.Nx; ++i){
		for (j = 1; j < sS.Ny; ++j){
			id = sS.ID[i][j];
			if (id == 2){
				temp = fF.Ex[i][j];
				DLM[ii].temp1x = 0.0;
				for (p=0;p<num_of_D_poles;p++){
					DLM[ii].temp1x = DLM[ii].temp1x + C4[p]*DLM[ii].JDx[p];
				}
				DLM[ii].temp2x = 0.0;
				for (p=0;p<num_of_L_poles;p++){
					DLM[ii].temp2x = DLM[ii].temp2x + C5[p]*DLM[ii].JLx[p] + C6[p]*DLM[ii].JLx_1[p];
				}
				fF.Ex[i][j] = C1 * fF.Ex[i][j] + C2 * DLM[ii].Ex_1 + C3 * (fF.Hz[i][j] - fF.Hz[i][j - 1]) * sS.cExy[j] - DLM[ii].temp1x - DLM[ii].temp2x;				
				
				for (p=0;p<num_of_D_poles;p++){
					DLM[ii].JDx[p] = kappa[p] * DLM[ii].JDx[p] + beta[p] * (fF.Ex[i][j] + temp);					
				}
				for (p=0;p<num_of_L_poles;p++){
					//temp0 = DLM[ii].JLx_1[p];
					//DLM[ii].JLx_1[p] = DLM[ii].JLx[p];
					temp0 = DLM[ii].JLx[p];
					DLM[ii].JLx[p] = alpha[p] * DLM[ii].JLx[p] + xi[p] * DLM[ii].JLx_1[p] + gamma_L[p] * (fF.Ex[i][j] - DLM[ii].Ex_1) / 2.0 / sS.dt;
					DLM[ii].JLx_1[p] = temp0;
				}				
				DLM[ii].Ex_1 = temp;
				
				temp = fF.Ey[i][j];
				DLM[ii].temp1y = 0.0;
				for (p=0;p<num_of_D_poles;p++){
					DLM[ii].temp1y = DLM[ii].temp1y + C4[p]*DLM[ii].JDy[p];
				}
				DLM[ii].temp2y = 0.0;
				for (p=0;p<num_of_L_poles;p++){
					DLM[ii].temp2y = DLM[ii].temp2y + C5[p]*DLM[ii].JLy[p] + C6[p]*DLM[ii].JLy_1[p];
				}
				fF.Ey[i][j] = C1 * fF.Ey[i][j] + C2 * DLM[ii].Ey_1 - C3 * (fF.Hz[i][j] - fF.Hz[i - 1][j]) * sS.cEyx[i] - DLM[ii].temp1y - DLM[ii].temp2y;				
				for (p=0;p<num_of_D_poles;p++){
					DLM[ii].JDy[p] = kappa[p] * DLM[ii].JDy[p] + beta[p] * (fF.Ey[i][j] + temp);
				}				
				for (p=0;p<num_of_L_poles;p++){
					//temp0 = DLM[ii].JLy_1[p];
					//DLM[ii].JLy_1[p] = DLM[ii].JLy[p];
					temp0 = DLM[ii].JLy[p];
					DLM[ii].JLy[p] = alpha[p] * DLM[ii].JLy[p] + xi[p] * DLM[ii].JLy_1[p] + gamma_L[p] * (fF.Ey[i][j] - DLM[ii].Ey_1) / 2.0 / sS.dt;
					DLM[ii].JLy_1[p] = temp0;
				}				
				DLM[ii].Ey_1 = temp;				
				ii = ii + 1;
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
