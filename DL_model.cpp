#include "stdafx.h"
#include "DL_model.h"

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

DL_model::DL_model(Structure &sS, char *fileName)
{
	int i, j, k, ii;
	cout << BOLDRED << ">" << RESET << "Constructing " << BOLDBLUE << "Drude–Lorentz " << RESET << "material...";
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
	iss >> epsInf;
	//cout<<epsInf<<"\t";
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> omega_D;
	//cout<<omega_D<<"\t";
	//omega_D = 2.0*pi*omega_D;
	iss >> gamma_D;
	//cout<<gamma_D<<"\t";
	//gamma_D = 2.0*pi*gamma_D;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> omega_L;
	//cout<<omega_L<<"\t";
	//omega_L = 2.0*pi*omega_L;
	iss >> gamma_L;
	//cout<<gamma_L<<endl;
	//gamma_L = 2.0*pi*gamma_L;
	iss >> delta_eps_L;

	alpha = (1.0-sS.dt*gamma_D/2.0)/(1.0+sS.dt*gamma_D/2.0);
	beta = sS.dt*omega_D*omega_D*eps0/2.0/(1.0+sS.dt*gamma_D/2.0);
	zeta = sS.dt*sS.dt*delta_eps_L*omega_L*omega_L/2.0/(1.0+sS.dt*gamma_L+sS.dt*sS.dt*omega_L*omega_L/2.0);
	rho = 1.0/(1.0+sS.dt*gamma_L+sS.dt*sS.dt*omega_L*omega_L/2.0);
	tau = (2.0+sS.dt*gamma_L-sS.dt*sS.dt*omega_L*omega_L/2.0)/(1.0+sS.dt*gamma_L+sS.dt*sS.dt*omega_L*omega_L/2.0);
	OMEGA = eps0*zeta/(eps0*epsInf) + 1.0 + sS.dt*beta/(2.0*eps0*epsInf);

	C1 = 1.0/OMEGA;
	C2 = sS.dt/(OMEGA*eps0*epsInf);
	C3 = -sS.dt*alpha/(2.0*OMEGA*eps0*epsInf);
	C4 = -sS.dt*beta/(2.0*OMEGA*eps0*epsInf);
	C5 = -sS.dt/(2.0*OMEGA*eps0*epsInf);
	C6 = -eps0*zeta/(OMEGA*eps0*epsInf);
	C7 = -eps0*tau/(OMEGA*eps0*epsInf);
	C8 = eps0*rho/(OMEGA*eps0*epsInf);
	C9 = eps0/(OMEGA*eps0*epsInf);

	ii = 0;
	for (i = 0; i < sS.Nx; ++i){
		for (j = 0; j < sS.Ny; ++j){
			for (k = 0; k < sS.Nz; ++k){
				if (sS.ID[i][j][k] == 2){
					ii = ii + 1;
				}
			}
		}
	}
	num_of_DLM_gridcells = ii;
	DLM = new DL_material[num_of_DLM_gridcells];
	ii = 0;
	for (i = 0; i < sS.Nx; ++i){
		for (j = 0; j < sS.Ny; ++j){
			for (k = 0; k < sS.Nz; ++k){
				if (sS.ID[i][j][k] == 2){
					DLM[ii].ix = i;
					DLM[ii].jy = j;
					DLM[ii].kz = k;
					DLM[ii].Jx = 0.0;
					DLM[ii].Jy = 0.0;
					DLM[ii].Jz = 0.0;
					DLM[ii].Px = 0.0;
					DLM[ii].Py = 0.0;
					DLM[ii].Pz = 0.0;
					DLM[ii].Px_1 = 0.0;
					DLM[ii].Py_1 = 0.0;
					DLM[ii].Pz_1 = 0.0;
					DLM[ii].Ex = 0.0;
					DLM[ii].Ey = 0.0;
					DLM[ii].Ez = 0.0;
					DLM[ii].Ex_1 = 0.0;
					DLM[ii].Ey_1 = 0.0;
					DLM[ii].Ez_1 = 0.0;
					ii = ii + 1;
				}
			}
		}
	}
	cout << " Done." << endl;
}

DL_model::~DL_model()
{
}

void DL_model::EFieldUpdate(Field &fF, Structure &sS)
{	//short id;
	int i;//, j, k;
	double tem1, tem2, tem3;


//#pragma omp parallel for default(shared) private(i,j,k,id)
//	for (i = 1; i < sS.Nx; ++i){
//		for (j = 1; j < sS.Ny; ++j){
//			for (k = 0; k < sS.Nz; ++k){
//				id = sS.ID[i][j][k];
//				if (id == 4){
//					fF.Ex[i][j][k] = 0.0;
//					fF.Ey[i][j][k] = 0.0;
//					fF.Ez[i][j][k] = 0.0;
//				}
//				else if(id != 2){
//					fF.Ex[i][j][k] = sS.aE[id] * fF.Ex[i][j][k] + sS.bE[id] * (fF.Hz[i][j][k] - fF.Hz[i][j - 1][k]) * sS.cExz[j] - sS.bE[id] * (fF.Hy[i][j][k] - fF.Hy[i][j][k - 1]) * sS.cExy[k];
//					fF.Ey[i][j][k] = sS.aE[id] * fF.Ey[i][j][k] + sS.bE[id] * (fF.Hx[i][j][k] - fF.Hx[i][j][k - 1]) * sS.cEyx[k] - sS.bE[id] * (fF.Hz[i][j][k] - fF.Hz[i - 1][j][k]) * sS.cEyz[i];
//					fF.Ez[i][j][k] = sS.aE[id] * fF.Ez[i][j][k] + sS.bE[id] * (fF.Hy[i][j][k] - fF.Hy[i - 1][j][k]) * sS.cEzy[i] - sS.bE[id] * (fF.Hx[i][j][k] - fF.Hx[i][j - 1][k]) * sS.cEzx[j];
//				}
//			}
//		}
//	}
#pragma omp parallel for default(shared) private(i,tem1,tem2,tem3)
	for(i=0; i<num_of_DLM_gridcells; ++i)
	{
		DLM[i].Ex = fF.Ex[DLM[i].ix][DLM[i].jy][DLM[i].kz];
		DLM[i].Ey = fF.Ey[DLM[i].ix][DLM[i].jy][DLM[i].kz];
		DLM[i].Ez = fF.Ez[DLM[i].ix][DLM[i].jy][DLM[i].kz];
		fF.Ex[DLM[i].ix][DLM[i].jy][DLM[i].kz] = C1 * fF.Ex[DLM[i].ix][DLM[i].jy][DLM[i].kz] + C2 * (fF.Hz[DLM[i].ix][DLM[i].jy][DLM[i].kz] - fF.Hz[DLM[i].ix][DLM[i].jy - 1][DLM[i].kz]) * sS.cExz[DLM[i].jy] - C2 * (fF.Hy[DLM[i].ix][DLM[i].jy][DLM[i].kz] - fF.Hy[DLM[i].ix][DLM[i].jy][DLM[i].kz - 1]) * sS.cExy[DLM[i].kz] + C3 * DLM[i].Jx + C4 * fF.Ex[DLM[i].ix][DLM[i].jy][DLM[i].kz] + C5 * DLM[i].Jx + C6 * fF.Ex[DLM[i].ix][DLM[i].jy][DLM[i].kz] + C7 * DLM[i].Px + C8 * DLM[i].Px_1 + C9 * DLM[i].Px;
		fF.Ey[DLM[i].ix][DLM[i].jy][DLM[i].kz] = C1 * fF.Ey[DLM[i].ix][DLM[i].jy][DLM[i].kz] + C2 * (fF.Hx[DLM[i].ix][DLM[i].jy][DLM[i].kz] - fF.Hx[DLM[i].ix][DLM[i].jy][DLM[i].kz - 1]) * sS.cEyx[DLM[i].kz] - C2 * (fF.Hz[DLM[i].ix][DLM[i].jy][DLM[i].kz] - fF.Hz[DLM[i].ix - 1][DLM[i].jy][DLM[i].kz]) * sS.cEyz[DLM[i].ix] + C3 * DLM[i].Jy + C4 * fF.Ey[DLM[i].ix][DLM[i].jy][DLM[i].kz] + C5 * DLM[i].Jy + C6 * fF.Ey[DLM[i].ix][DLM[i].jy][DLM[i].kz] + C7 * DLM[i].Py + C8 * DLM[i].Py_1 + C9 * DLM[i].Py;
		fF.Ez[DLM[i].ix][DLM[i].jy][DLM[i].kz] = C1 * fF.Ez[DLM[i].ix][DLM[i].jy][DLM[i].kz] + C2 * (fF.Hy[DLM[i].ix][DLM[i].jy][DLM[i].kz] - fF.Hy[DLM[i].ix - 1][DLM[i].jy][DLM[i].kz]) * sS.cEzy[DLM[i].ix] - C2 * (fF.Hx[DLM[i].ix][DLM[i].jy][DLM[i].kz] - fF.Hx[DLM[i].ix][DLM[i].jy - 1][DLM[i].kz]) * sS.cEzx[DLM[i].jy] + C3 * DLM[i].Jz + C4 * fF.Ez[DLM[i].ix][DLM[i].jy][DLM[i].kz] + C5 * DLM[i].Jz + C6 * fF.Ez[DLM[i].ix][DLM[i].jy][DLM[i].kz] + C7 * DLM[i].Pz + C8 * DLM[i].Pz_1 + C9 * DLM[i].Pz;
		DLM[i].Jx = alpha * DLM[i].Jx + beta * (fF.Ex[DLM[i].ix][DLM[i].jy][DLM[i].kz] + DLM[i].Ex_1);
		DLM[i].Jy = alpha * DLM[i].Jy + beta * (fF.Ey[DLM[i].ix][DLM[i].jy][DLM[i].kz] + DLM[i].Ey_1);
		DLM[i].Jz = alpha * DLM[i].Jz + beta * (fF.Ez[DLM[i].ix][DLM[i].jy][DLM[i].kz] + DLM[i].Ez_1);
		tem1 = DLM[i].Px;
		tem2 = DLM[i].Py;
		tem3 = DLM[i].Pz;
		DLM[i].Px = zeta * (fF.Ex[DLM[i].ix][DLM[i].jy][DLM[i].kz] + DLM[i].Ex) + tau * DLM[i].Px - rho * DLM[i].Px_1;
		DLM[i].Py = zeta * (fF.Ey[DLM[i].ix][DLM[i].jy][DLM[i].kz] + DLM[i].Ey) + tau * DLM[i].Py - rho * DLM[i].Py_1;
		DLM[i].Pz = zeta * (fF.Ez[DLM[i].ix][DLM[i].jy][DLM[i].kz] + DLM[i].Ez) + tau * DLM[i].Pz - rho * DLM[i].Pz_1;
		DLM[i].Px_1 = tem1;
		DLM[i].Py_1 = tem2;
		DLM[i].Pz_1 = tem3;
	}

//	id = 0;
//	i = 0;
//	//#pragma omp parallel for default(shared) private(j,k)
//	for (k = 1; k<(sS.Nz-1); ++k){
//		for (j = 1; j<(sS.Ny-1); ++j){
//			fF.Ex[i][j][k] = sS.aE[id] * fF.Ex[i][j][k] + sS.bE[id] * (fF.Hz[i][j][k] - fF.Hz[i][j - 1][k]) * sS.cExz[j] - sS.bE[id] * (fF.Hy[i][j][k] - fF.Hy[i][j][k - 1]) * sS.cExy[k];
//		}
//	}
//	j = 0;
//	//#pragma omp parallel for default(shared) private(i,k)
//	for (k = 1; k<(sS.Nz-1); ++k){
//		for (i = 1; i<(sS.Nx-1); ++i){
//			fF.Ey[i][j][k] = sS.aE[id] * fF.Ey[i][j][k] + sS.bE[id] * (fF.Hx[i][j][k] - fF.Hx[i][j][k - 1]) * sS.cEyx[k] - sS.bE[id] * (fF.Hz[i][j][k] - fF.Hz[i - 1][j][k]) * sS.cEyz[i];
//		}
//	}
//	k = 0;
//	//#pragma omp parallel for default(shared) private(i,j)
//	for (j = 1; j<(sS.Ny-1); ++j){
//		for (i = 1; i<(sS.Nx-1); ++i){
//			fF.Ez[i][j][k] = sS.aE[id] * fF.Ez[i][j][k] + sS.bE[id] * (fF.Hy[i][j][k] - fF.Hy[i - 1][j][k]) * sS.cEzy[i] - sS.bE[id] * (fF.Hx[i][j][k] - fF.Hx[i][j - 1][k]) * sS.cEzx[j];
//		}
//	}
}
