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
#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "Structure.h"
#include "Fields.h"



double ContinuousWave(double Amplitude, double retardedtime, int nd, int nw, double pinudt);
double Pulse(double Amplitude, double retardedtime, int nd, int nw, double pinudt);
double GaussianWaveSpatialDistribution(double r, double C1, double C2, double C3);

class Gaussian
{
public:
	Gaussian(Structure &sS, char *fileName);
	~Gaussian();


	void GaussianPulseUpdate_z_x(int timestep, Field &fF, Structure &sS);
	void GaussianCWUpdate_z_x(int timestep, Field &fF, Structure &sS);
	void GaussianPulseUpdate_z_y(int timestep, Field &fF, Structure &sS);
	void GaussianCWUpdate_z_y(int timestep, Field &fF, Structure &sS);
	void GaussianPulseUpdate_y_x(int timestep, Field &fF, Structure &sS);
	void GaussianCWUpdate_y_x(int timestep, Field &fF, Structure &sS);
	void GaussianPulseUpdate_y_z(int timestep, Field &fF, Structure &sS);
	void GaussianCWUpdate_y_z(int timestep, Field &fF, Structure &sS);
	void GaussianPulseUpdate_x_y(int timestep, Field &fF, Structure &sS);
	void GaussianCWUpdate_x_y(int timestep, Field &fF, Structure &sS);
	void GaussianPulseUpdate_x_z(int timestep, Field &fF, Structure &sS);
	void GaussianCWUpdate_x_z(int timestep, Field &fF, Structure &sS);
	void InitializeGaussianSource();

	double lambda, nu, T, K;
	double amp, ampH;
	double l0, w0, z0;
	int i1TFSF, i2TFSF, j1TFSF, j2TFSF, k1TFSF, k2TFSF;
	char direction, polarization;

private:
	double dt, dx, dy, dz, r, r01, r02;
	double retard, pinudt, Cl1, Cl2, Cl3, P;
	int nd, nw, nT, nTrans, nstep;
	double **Er, Einc, Hinc;
};

#endif
