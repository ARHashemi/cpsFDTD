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
#include "Gaussian.h"
#include "Constants.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;


Gaussian::Gaussian(Structure &sS, char *fileName)
{
	cout << "Constructing Gaussian Source...";
	dt = sS.dt;
	dx = sS.dx;
	dy = sS.dy;
	dz = sS.dz;
	ifstream InputFile;
	InputFile.open(fileName);
	string sline, id;
	getline(InputFile, sline);
	istringstream iss(sline);
	iss >> lambda;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> w0;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> l0;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> i1TFSF; iss >> i2TFSF;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> j1TFSF; iss >> j2TFSF;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> k1TFSF; iss >> k2TFSF;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> polarization; iss >> direction;


	cout << " Done." << endl;
}

void Gaussian::InitializeGaussianSource()
{
	cout << "Initializing Gaussian Source...";
	//pulse data (cw data)
	nu = C / lambda;
	T = 1.0 / nu;
	nT = T / dt;
	pinudt = 2.0 * pi * nu * dt;
	nw = 2 * nT;
	nd = 4 * nw;
	nTrans = nd + nw;

	amp = 10.0; ampH = -amp*sqrt(eps0 / mu0);

	//Gaussian data
	K = 2.0 * pi / lambda;
	z0 = pi * w0*w0 / lambda;
	Cl1 = sin(K*l0) * z0 * l0 / (z0*z0 + l0*l0);
	Cl2 = K * l0 / (2 * (z0* z0 + l0*l0));
	Cl3 = -z0*z0 / w0*w0 * (z0*z0 + l0*l0);


	if (direction == 'z'){
		Er = new double*[i2TFSF - i1TFSF + 1];
		r01 = double(i1TFSF + i2TFSF) * dx / 2.0;
		r02 = double(j1TFSF + j2TFSF) * dy / 2.0;
		for (int i = i1TFSF; i <= i2TFSF; ++i){
			Er[i - i1TFSF] = new double[j2TFSF - j1TFSF + 1];
			for (int j = j1TFSF; j <= j2TFSF; ++j){
				Er[i - i1TFSF][j - j1TFSF] = GaussianWaveSpatialDistribution(sqrt((i*dx - r01)*(i*dx - r01) + (j*dy - r02)*(j*dy - r02)), Cl1, Cl2, Cl3);
			}
		}
	} else if (direction == 'y'){
		Er = new double*[i2TFSF - i1TFSF + 1];
		r01 = double(i1TFSF + i2TFSF) * dx / 2.0;
		r02 = double(k1TFSF + k2TFSF) * dz / 2.0;
		for (int i = i1TFSF; i <= i2TFSF; ++i){
			Er[i - i1TFSF] = new double[k2TFSF - k1TFSF + 1];
			for (int k = k1TFSF; k <= j2TFSF; ++k){
				Er[i - i1TFSF][k - k1TFSF] = GaussianWaveSpatialDistribution(sqrt((i*dx - r01)*(i*dx - r01) + (k*dz - r02)*(k*dz - r02)), Cl1, Cl2, Cl3);
			}
		}
	} else if (direction == 'x'){
		Er = new double*[j2TFSF - j1TFSF + 1];
		r01 = double(j1TFSF + j2TFSF) * dy / 2.0;
		r02 = double(k1TFSF + k2TFSF) * dz / 2.0;
		for (int j = j1TFSF; j <= j2TFSF; ++j){
			Er[j - j1TFSF] = new double[k2TFSF - k1TFSF + 1];
			for (int k = k1TFSF; k <= j2TFSF; ++k){
				Er[j - j1TFSF][k - k1TFSF] = GaussianWaveSpatialDistribution(sqrt((j*dy - r01)*(j*dy - r01) + (k*dz - r02)*(k*dz - r02)), Cl1, Cl2, Cl3);
			}
		}
	}
	cout << " Done." << endl;
}

Gaussian::~Gaussian()
{
	cout << "Deconstructing Gaussian Source...";
	delete[] Er;
	Er = 0;
	cout << " Done." << endl;
}


double ContinuousWave(double Amplitude, double retardedtime, int nd, int nw, double pinudt)
{
	//retard = ntimestep - location / C / dt;
	if (retardedtime < nd){
		return (Amplitude * exp(-(retardedtime - nd)*(retardedtime - nd) / (nw*nw)) * sin(pinudt*retardedtime));
	}
	else{
		return (Amplitude * sin(pinudt*retardedtime));
	}
}

double Pulse(double Amplitude, double retardedtime, int nd, int nw, double pinudt)
{
	//retard = ntimestep - location / C / dt;
	return (Amplitude * exp(-double(retardedtime - nd)*(retardedtime - nd) / (nw*nw)) * sin(pinudt*(retardedtime - nd)));
}

double GaussianWaveSpatialDistribution(double r, double C1, double C2, double C3)
{
	return (C1 * sin(C2 * r*r) * exp(C3 * r*r));
}

void Gaussian::GaussianCWUpdate_z_x(int timestep, Field &fF, Structure &sS)
{
	retard = timestep - 0.5 - (j1TFSF - 0.5)*sS.dy / C / dt;
	Hinc = ContinuousWave(ampH, retard, nd, nw, pinudt);
	retard = timestep - (j1TFSF)*sS.dy / C / dt;
	Einc = ContinuousWave(amp, retard, nd, nw, pinudt);
	for (int i = i1TFSF; i <= i2TFSF; ++i){
		for (int j = j1TFSF; j <= j2TFSF; ++j){
			fF.Ex[i][j][k1TFSF] = fF.Ex[i][j][k1TFSF] - sS.bE[0] * Hinc * Er[i - i1TFSF][j - j1TFSF] / sS.dy;
			//Ex[i][j2TFSF] = Ex[i][j2TFSF] + bE[0] * Hzinc[j2TFSF] / dy; * Er[i - i1TFSF]
			fF.Hz[i][j - 1][k1TFSF] = fF.Hz[i][j - 1][k1TFSF] - sS.bH[0] * Einc * Er[i - i1TFSF][j - j1TFSF] / sS.dy;
			//Hz[i][j2TFSF] = Hz[i][j2TFSF] + bH[0] * Exinc[j2TFSF] / dy; * Er[i - i1TFSF]
		}
	}
}

void Gaussian::GaussianPulseUpdate_z_x(int timestep, Field &fF, Structure &sS)
{
	nstep = timestep % (7 * nw);//
	retard = nstep - 0.5 - (j1TFSF - 0.5)*sS.dy / C / dt;
	Hinc = Pulse(ampH, retard, nd, nw, pinudt);
	retard = nstep - (j1TFSF)*sS.dy / C / dt;
	Einc = Pulse(amp, retard, nd, nw, pinudt);
	for (int i = i1TFSF; i <= i2TFSF; ++i){
		for (int j = j1TFSF; j <= j2TFSF; ++j){
			fF.Ex[i][j][k1TFSF] = fF.Ex[i][j][k1TFSF] - sS.bE[0] * Hinc * Er[i - i1TFSF][j - j1TFSF] / sS.dy;
			//Ex[i][j2TFSF] = Ex[i][j2TFSF] + bE[0] * Hzinc[j2TFSF] / dy; * Er[i - i1TFSF]
			fF.Hz[i][j - 1][k1TFSF] = fF.Hz[i][j - 1][k1TFSF] - sS.bH[0] * Einc * Er[i - i1TFSF][j - j1TFSF] / sS.dy;
			//Hz[i][j2TFSF] = Hz[i][j2TFSF] + bH[0] * Exinc[j2TFSF] / dy; * Er[i - i1TFSF]
		}
	}
}

void Gaussian::GaussianCWUpdate_z_y(int timestep, Field &fF, Structure &sS)
{
	retard = timestep - 0.5 - (j1TFSF - 0.5)*sS.dy / C / dt;
	Hinc = ContinuousWave(ampH, retard, nd, nw, pinudt);
	retard = timestep - (j1TFSF)*sS.dy / C / dt;
	Einc = ContinuousWave(amp, retard, nd, nw, pinudt);
	for (int i = i1TFSF; i <= i2TFSF; ++i){
		for (int j = j1TFSF; j <= j2TFSF; ++j){
			fF.Ex[i][j][k1TFSF] = fF.Ex[i][j][k1TFSF] - sS.bE[0] * Hinc * Er[i - i1TFSF][j - j1TFSF] / sS.dy;
			//Ex[i][j2TFSF] = Ex[i][j2TFSF] + bE[0] * Hzinc[j2TFSF] / dy; * Er[i - i1TFSF]
			fF.Hz[i][j - 1][k1TFSF] = fF.Hz[i][j - 1][k1TFSF] - sS.bH[0] * Einc * Er[i - i1TFSF][j - j1TFSF] / sS.dy;
			//Hz[i][j2TFSF] = Hz[i][j2TFSF] + bH[0] * Exinc[j2TFSF] / dy; * Er[i - i1TFSF]
		}
	}
}
