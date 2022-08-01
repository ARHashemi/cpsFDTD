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
#ifndef TFSF_H
#define TFSF_H
#include "Structure.h"
#include "Fields.h"
#include <string.h>

struct TFSFpoint{
	int i, j, k;
	//double x0E, y0E, z0E;
	//double xpE, ypE;
    double r0E1, r0E2;
    //double x0H, y0H, z0H;
	//double xpH, ypH;
    double r0H1, r0H2;
    //double kr;
    //double sinkr, coskr;
    double fpE1, fpH1, fpE2, fpH2;
};

double GaussianWaveSpatialDistribution(double r2, double C1, double C2, double C3);



class TFSF
{
public:
	TFSF(Structure &sS, const char *fileName);
	~TFSF();
	void InitializeSource();
	void UpdateFields(int timestep, Field &fF, Structure &sS);

	std::string temporal_mode, spatial_mode;
	int Nt;
	double lambda, dlambda, nu, dnu, T, K, eta;
	double amp, ampEx, ampEy, ampEz, ampHx, ampHy, ampHz;
	double l0, w0, z0;
	int i1TFSF, i2TFSF, j1TFSF, j2TFSF, k1TFSF, k2TFSF;
	double polarization_angle, theta, phi;
	int ixc, jyc, kzc;
	int numPoints;
	TFSFpoint *TFSFp;
	double pinudt, *Ft;

private:
	double Fout(int n, double r);
	double dt, dx, dy, dz, r, r01, r02;
	double retard, Cl1, Cl2, Cl3, P, sig;
	int nd, nw, nT, nTrans, Np, nFloor;
	double nstep, nprime;
	int ii1, ii2, ii3, ii4, ii5, ii6;
	double **Er, Einc, Hinc;
};

#endif // TFSF_H
