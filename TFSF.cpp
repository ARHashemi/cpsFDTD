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
#include "TFSF.h"
#include "stdafx.h"
#include "Constants.h"
#include "Structure.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

TFSF::TFSF(Structure &sS, const char *fileName)
{
	cout << "Constructing TFSF...";
	dt = sS.dt;
	dx = sS.dx;
	dy = sS.dy;
	dz = sS.dz;
	ifstream InputFile;
	InputFile.open(fileName);
	string sline, id;
	getline(InputFile, sline);
	istringstream iss(sline);
	iss >> Nt;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> lambda;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> dlambda;
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
	iss >> amp;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> temporal_mode;
	iss >> spatial_mode;
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
	iss >> polarization_angle; iss >> theta; iss >> phi;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> ixc; iss >> jyc; iss >> kzc;

	K = 2.0*pi/lambda;
	nu = C / lambda;
	T = 1.0 / nu;
	nT = T / dt;
	pinudt = 2.0 * pi * nu * dt;
	cout << " Done." << endl;
}

TFSF::~TFSF()
{
	//dtor
}

void TFSF::InitializeSource()
{
	cout << "Initilizing Source...";
	int i, j, k;



	dnu = C / (lambda - dlambda/2.0) - C / (lambda + dlambda/2.0);
	sig = dnu/(2.0*sqrt(2.0*log(2.0)));
	nd = 4.0/sqrt(2.0*sig*sig*pi*pi*dt*dt);
	Np = 2.0*nd;

	Ft = new double[Nt];
	if (temporal_mode == "singlepulse"){
		for (i=0; i<Nt; ++i){
			Ft[i] = exp(-2.0*sig*sig*pi*pi*dt*dt*(i-nd)*(i-nd)) * sin(pinudt * (i));//-nd
		}
	}else if (temporal_mode == "pulsetrain"){
		for (i=0; i<Nt; ++i){
			Ft[i] = exp(-2.0*sig*sig*pi*pi*dt*dt*((i%Np)-nd)*((i%Np)-nd)) * sin(pinudt * ((i%Np)));//-nd
		}
	}else if (temporal_mode == "cw"){
		for (i=0; i<nd; ++i){
			Ft[i] = exp(-((double)i-8.0*nT)*((double)i-8.0*nT)/(4.0*nT*nT)) * sin(pinudt * (i));//-nd  exp(-2.0*sig*sig*pi*pi*dt*dt*(n-nd)*(n-nd))
		}
		for (i=nd; i<Nt; ++i){
			Ft[i] = sin(pinudt * (i));//-nd
		}
	}

	//nw = 2 * nT;
	//nd = 4 * nw;
	//nTrans = nd + nw;

	z0 = pi * w0*w0 / lambda;
	Cl1 = sin(K*l0) * z0 * l0 / (z0*z0 + l0*l0);
	Cl2 = K * l0 / (2 * (z0* z0 + l0*l0));
	Cl3 = -z0*z0 / w0*w0 * (z0*z0 + l0*l0);

	numPoints = 2*(i2TFSF-i1TFSF+1)*(j2TFSF-j1TFSF+1) + 2*(k2TFSF-k1TFSF+1)*(j2TFSF-j1TFSF+1) + 2*(i2TFSF-i1TFSF+1)*(k2TFSF-k1TFSF+1);
	TFSFp = new TFSFpoint[numPoints+1];
	double x0, y0, z0, xp, yp;
	int ii = 0;
	k = k1TFSF;
	for (i=i1TFSF; i<i2TFSF; i++){
		for (j=j1TFSF; j<j2TFSF; j++){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			TFSFp[ii].r0E1 = sin(theta)*cos(phi)*(i+0.5-ixc)*dx + sin(theta)*sin(phi)*(j-jyc)*dy + cos(theta)*(k-kzc)*dz;
            x0 = sin(theta)*cos(phi)*TFSFp[ii].r0E1;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0E1;
            z0 = cos(theta)*TFSFp[ii].r0E1;
            xp = -((i+0.5-ixc)*dx - x0)*sin(phi) + ((j-jyc)*dy - y0)*cos(phi);
            yp = -((i+0.5-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpE1 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0E2 = sin(theta)*cos(phi)*(i-ixc)*dx + sin(theta)*sin(phi)*(j+0.5-jyc)*dy + cos(theta)*(k-kzc)*dz;
			x0 = sin(theta)*cos(phi)*TFSFp[ii].r0E2;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0E2;
            z0 = cos(theta)*TFSFp[ii].r0E2;
            xp = -((i-ixc)*dx - x0)*sin(phi) + ((j+0.5-jyc)*dy - y0)*cos(phi);
            yp = -((i-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j+0.5-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpE2 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0H1 = sin(theta)*cos(phi)*(i+0.5-ixc)*dx + sin(theta)*sin(phi)*(j-jyc)*dy + cos(theta)*(k-0.5-kzc)*dz;
            x0 = sin(theta)*cos(phi)*TFSFp[ii].r0H1;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0H1;
            z0 = cos(theta)*TFSFp[ii].r0H1;
            xp = -((i+0.5-ixc)*dx - x0)*sin(phi) + ((j-jyc)*dy - y0)*cos(phi);
            yp = -((i+0.5-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-0.5-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpH1 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0H2 = sin(theta)*cos(phi)*(i-ixc)*dx + sin(theta)*sin(phi)*(j+0.5-jyc)*dy + cos(theta)*(k-0.5-kzc)*dz;
			x0 = sin(theta)*cos(phi)*TFSFp[ii].r0H2;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0H2;
            z0 = cos(theta)*TFSFp[ii].r0H2;
            xp = -((i-ixc)*dx - x0)*sin(phi) + ((j+0.5-jyc)*dy - y0)*cos(phi);
            yp = -((i-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j+0.5-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-0.5-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpH2 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);
            //TFSFp[ii].kr = K*TFSFp[ii].r0;
            //TFSFp[ii].coskr = cos(K*TFSFp[ii].r0);
            //TFSFp[ii].sinkr = sin(K*TFSFp[ii].r0);
            //if (spatial_mode=="Gaussian"){
			//	TFSFp[ii].fpE = GaussianWaveSpatialDistribution((TFSFp[ii].xpE*TFSFp[ii].xpE+TFSFp[ii].ypE*TFSFp[ii].ypE),Cl1,Cl2,Cl3);
            //}
            ii = ii + 1;
		}
	}
	ii1 = ii;
	k = k2TFSF;
	for (i=i1TFSF; i<i2TFSF; i++){
		for (j=j1TFSF; j<j2TFSF; j++){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			TFSFp[ii].r0E1 = sin(theta)*cos(phi)*(i+0.5-ixc)*dx + sin(theta)*sin(phi)*(j-jyc)*dy + cos(theta)*(k-kzc)*dz;
            x0 = sin(theta)*cos(phi)*TFSFp[ii].r0E1;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0E1;
            z0 = cos(theta)*TFSFp[ii].r0E1;
            xp = -((i+0.5-ixc)*dx - x0)*sin(phi) + ((j-jyc)*dy - y0)*cos(phi);
            yp = -((i+0.5-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpE1 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0E2 = sin(theta)*cos(phi)*(i-ixc)*dx + sin(theta)*sin(phi)*(j+0.5-jyc)*dy + cos(theta)*(k-kzc)*dz;
			x0 = sin(theta)*cos(phi)*TFSFp[ii].r0E2;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0E2;
            z0 = cos(theta)*TFSFp[ii].r0E2;
            xp = -((i-ixc)*dx - x0)*sin(phi) + ((j+0.5-jyc)*dy - y0)*cos(phi);
            yp = -((i-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j+0.5-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpE2 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0H1 = sin(theta)*cos(phi)*(i+0.5-ixc)*dx + sin(theta)*sin(phi)*(j-jyc)*dy + cos(theta)*(k+0.5-kzc)*dz;
            x0 = sin(theta)*cos(phi)*TFSFp[ii].r0H1;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0H1;
            z0 = cos(theta)*TFSFp[ii].r0H1;
            xp = -((i+0.5-ixc)*dx - x0)*sin(phi) + ((j-jyc)*dy - y0)*cos(phi);
            yp = -((i+0.5-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k+0.5-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpH1 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0H2 = sin(theta)*cos(phi)*(i-ixc)*dx + sin(theta)*sin(phi)*(j+0.5-jyc)*dy + cos(theta)*(k+0.5-kzc)*dz;
			x0 = sin(theta)*cos(phi)*TFSFp[ii].r0H2;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0H2;
            z0 = cos(theta)*TFSFp[ii].r0H2;
            xp = -((i-ixc)*dx - x0)*sin(phi) + ((j+0.5-jyc)*dy - y0)*cos(phi);
            yp = -((i-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j+0.5-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k+0.5-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpH2 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);
            ii = ii + 1;
		}
	}
	ii2 = ii;
	j = j1TFSF;
	for (i=i1TFSF; i<i2TFSF; i++){
		for (k=k1TFSF; k<k2TFSF; k++){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			TFSFp[ii].r0E1 = sin(theta)*cos(phi)*(i+0.5-ixc)*dx + sin(theta)*sin(phi)*(j-jyc)*dy + cos(theta)*(k-kzc)*dz;
            x0 = sin(theta)*cos(phi)*TFSFp[ii].r0E1;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0E1;
            z0 = cos(theta)*TFSFp[ii].r0E1;
            xp = -((i+0.5-ixc)*dx - x0)*sin(phi) + ((j-jyc)*dy - y0)*cos(phi);
            yp = -((i+0.5-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpE1 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0E2 = sin(theta)*cos(phi)*(i-ixc)*dx + sin(theta)*sin(phi)*(j-jyc)*dy + cos(theta)*(k+0.5-kzc)*dz;
			x0 = sin(theta)*cos(phi)*TFSFp[ii].r0E2;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0E2;
            z0 = cos(theta)*TFSFp[ii].r0E2;
            xp = -((i-ixc)*dx - x0)*sin(phi) + ((j-jyc)*dy - y0)*cos(phi);
            yp = -((i-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k+0.5-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpE2 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0H1 = sin(theta)*cos(phi)*(i+0.5-ixc)*dx + sin(theta)*sin(phi)*(j-0.5-jyc)*dy + cos(theta)*(k-kzc)*dz;
            x0 = sin(theta)*cos(phi)*TFSFp[ii].r0H1;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0H1;
            z0 = cos(theta)*TFSFp[ii].r0H1;
            xp = -((i+0.5-ixc)*dx - x0)*sin(phi) + ((j-0.5-jyc)*dy - y0)*cos(phi);
            yp = -((i+0.5-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-0.5-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpH1 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0H2 = sin(theta)*cos(phi)*(i-ixc)*dx + sin(theta)*sin(phi)*(j-0.5-jyc)*dy + cos(theta)*(k+0.5-kzc)*dz;
			x0 = sin(theta)*cos(phi)*TFSFp[ii].r0H2;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0H2;
            z0 = cos(theta)*TFSFp[ii].r0H2;
            xp = -((i-ixc)*dx - x0)*sin(phi) + ((j-0.5-jyc)*dy - y0)*cos(phi);
            yp = -((i-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-0.5-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k+0.5-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpH2 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);
            ii = ii + 1;
		}
	}
	ii3 = ii;
	j = j2TFSF;
	for (i=i1TFSF; i<i2TFSF; i++){
		for (k=k1TFSF; k<k2TFSF; k++){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			TFSFp[ii].r0E1 = sin(theta)*cos(phi)*(i+0.5-ixc)*dx + sin(theta)*sin(phi)*(j-jyc)*dy + cos(theta)*(k-kzc)*dz;
            x0 = sin(theta)*cos(phi)*TFSFp[ii].r0E1;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0E1;
            z0 = cos(theta)*TFSFp[ii].r0E1;
            xp = -((i+0.5-ixc)*dx - x0)*sin(phi) + ((j-jyc)*dy - y0)*cos(phi);
            yp = -((i+0.5-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpE1 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0E2 = sin(theta)*cos(phi)*(i-ixc)*dx + sin(theta)*sin(phi)*(j-jyc)*dy + cos(theta)*(k+0.5-kzc)*dz;
			x0 = sin(theta)*cos(phi)*TFSFp[ii].r0E2;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0E2;
            z0 = cos(theta)*TFSFp[ii].r0E2;
            xp = -((i-ixc)*dx - x0)*sin(phi) + ((j-jyc)*dy - y0)*cos(phi);
            yp = -((i-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k+0.5-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpE2 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0H1 = sin(theta)*cos(phi)*(i+0.5-ixc)*dx + sin(theta)*sin(phi)*(j+0.5-jyc)*dy + cos(theta)*(k-kzc)*dz;
            x0 = sin(theta)*cos(phi)*TFSFp[ii].r0H1;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0H1;
            z0 = cos(theta)*TFSFp[ii].r0H1;
            xp = -((i+0.5-ixc)*dx - x0)*sin(phi) + ((j+0.5-jyc)*dy - y0)*cos(phi);
            yp = -((i+0.5-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j+0.5-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpH1 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0H2 = sin(theta)*cos(phi)*(i-ixc)*dx + sin(theta)*sin(phi)*(j+0.5-jyc)*dy + cos(theta)*(k+0.5-kzc)*dz;
			x0 = sin(theta)*cos(phi)*TFSFp[ii].r0H2;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0H2;
            z0 = cos(theta)*TFSFp[ii].r0H2;
            xp = -((i-ixc)*dx - x0)*sin(phi) + ((j+0.5-jyc)*dy - y0)*cos(phi);
            yp = -((i-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j+0.5-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k+0.5-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpH2 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);
            ii = ii + 1;
		}
	}
	ii4 = ii;
	i = i1TFSF;
	for (j=j1TFSF; j<j2TFSF; j++){
		for (k=k1TFSF; k<k2TFSF; k++){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			TFSFp[ii].r0E1 = sin(theta)*cos(phi)*(i-ixc)*dx + sin(theta)*sin(phi)*(j+0.5-jyc)*dy + cos(theta)*(k-kzc)*dz;
            x0 = sin(theta)*cos(phi)*TFSFp[ii].r0E1;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0E1;
            z0 = cos(theta)*TFSFp[ii].r0E1;
            xp = -((i-ixc)*dx - x0)*sin(phi) + ((j+0.5-jyc)*dy - y0)*cos(phi);
            yp = -((i-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j+0.5-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpE1 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0E2 = sin(theta)*cos(phi)*(i-ixc)*dx + sin(theta)*sin(phi)*(j-jyc)*dy + cos(theta)*(k+0.5-kzc)*dz;
			x0 = sin(theta)*cos(phi)*TFSFp[ii].r0E2;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0E2;
            z0 = cos(theta)*TFSFp[ii].r0E2;
            xp = -((i-ixc)*dx - x0)*sin(phi) + ((j-jyc)*dy - y0)*cos(phi);
            yp = -((i-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k+0.5-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpE2 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0H1 = sin(theta)*cos(phi)*(i-0.5-ixc)*dx + sin(theta)*sin(phi)*(j+0.5-jyc)*dy + cos(theta)*(k-kzc)*dz;
            x0 = sin(theta)*cos(phi)*TFSFp[ii].r0H1;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0H1;
            z0 = cos(theta)*TFSFp[ii].r0H1;
            xp = -((i-0.5-ixc)*dx - x0)*sin(phi) + ((j+0.5-jyc)*dy - y0)*cos(phi);
            yp = -((i-0.5-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j+0.5-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpH1 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0H2 = sin(theta)*cos(phi)*(i-0.5-ixc)*dx + sin(theta)*sin(phi)*(j-jyc)*dy + cos(theta)*(k+0.5-kzc)*dz;
			x0 = sin(theta)*cos(phi)*TFSFp[ii].r0H2;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0H2;
            z0 = cos(theta)*TFSFp[ii].r0H2;
            xp = -((i-0.5-ixc)*dx - x0)*sin(phi) + ((j-jyc)*dy - y0)*cos(phi);
            yp = -((i-0.5-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k+0.5-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpH2 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);
            ii = ii + 1;
		}
	}
	ii5 = ii;
	i = i2TFSF;
	for (j=j1TFSF; j<j2TFSF; j++){
		for (k=k1TFSF; k<k2TFSF; k++){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			TFSFp[ii].r0E1 = sin(theta)*cos(phi)*(i-ixc)*dx + sin(theta)*sin(phi)*(j+0.5-jyc)*dy + cos(theta)*(k-kzc)*dz;
            x0 = sin(theta)*cos(phi)*TFSFp[ii].r0E1;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0E1;
            z0 = cos(theta)*TFSFp[ii].r0E1;
            xp = -((i-ixc)*dx - x0)*sin(phi) + ((j+0.5-jyc)*dy - y0)*cos(phi);
            yp = -((i-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j+0.5-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpE1 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0E2 = sin(theta)*cos(phi)*(i-ixc)*dx + sin(theta)*sin(phi)*(j-jyc)*dy + cos(theta)*(k+0.5-kzc)*dz;
			x0 = sin(theta)*cos(phi)*TFSFp[ii].r0E2;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0E2;
            z0 = cos(theta)*TFSFp[ii].r0E2;
            xp = -((i-ixc)*dx - x0)*sin(phi) + ((j-jyc)*dy - y0)*cos(phi);
            yp = -((i-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k+0.5-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpE2 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0H1 = sin(theta)*cos(phi)*(i+0.5-ixc)*dx + sin(theta)*sin(phi)*(j+0.5-jyc)*dy + cos(theta)*(k-kzc)*dz;
            x0 = sin(theta)*cos(phi)*TFSFp[ii].r0H1;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0H1;
            z0 = cos(theta)*TFSFp[ii].r0H1;
            xp = -((i+0.5-ixc)*dx - x0)*sin(phi) + ((j+0.5-jyc)*dy - y0)*cos(phi);
            yp = -((i+0.5-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j+0.5-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpH1 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);

            TFSFp[ii].r0H2 = sin(theta)*cos(phi)*(i+0.5-ixc)*dx + sin(theta)*sin(phi)*(j-jyc)*dy + cos(theta)*(k+0.5-kzc)*dz;
			x0 = sin(theta)*cos(phi)*TFSFp[ii].r0H2;
            y0 = sin(theta)*sin(phi)*TFSFp[ii].r0H2;
            z0 = cos(theta)*TFSFp[ii].r0H2;
            xp = -((i+0.5-ixc)*dx - x0)*sin(phi) + ((j-jyc)*dy - y0)*cos(phi);
            yp = -((i+0.5-ixc)*dx - x0)*cos(phi)*cos(theta) - ((j-jyc)*dy - y0)*sin(phi)*cos(theta) + ((k+0.5-kzc)*dz - z0)*sin(theta);
            TFSFp[ii].fpH2 = GaussianWaveSpatialDistribution((xp*xp+yp*yp),Cl1,Cl2,Cl3);
            ii = ii + 1;
		}
	}
	ii6 = ii;
	ampEx = -amp*(sin(phi)*cos(polarization_angle)+cos(phi)*cos(theta)*sin(polarization_angle));
	ampEy = amp*(cos(phi)*cos(polarization_angle)-sin(phi)*cos(theta)*sin(polarization_angle));
	ampEz = amp*sin(theta)*sin(polarization_angle);
	eta = sqrt(mu0/eps0);
	ampHx = -eta*amp*(-sin(phi)*sin(polarization_angle)+cos(phi)*cos(theta)*cos(polarization_angle));
	ampHy = eta*amp*(-cos(phi)*sin(polarization_angle)-sin(phi)*cos(theta)*cos(polarization_angle));
	ampHz = eta*amp*sin(theta)*cos(polarization_angle);


	cout << " Done." << endl;
}

void TFSF::UpdateFields(int timestep, Field &fF, Structure &sS)
{
	int ii;
	for (ii=0; ii<ii1; ++ii){
        fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + ampHy * sS.bE[0] * Fout(timestep, TFSFp[ii].r0H1) / dz;
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - ampHx * sS.bE[0] * Fout(timestep, TFSFp[ii].r0H2) / dz;

		fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] + ampEx * sS.bH[0] * Fout(timestep, TFSFp[ii].r0E1) / dz;
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] - ampEy * sS.bH[0] * Fout(timestep, TFSFp[ii].r0E2) / dz;
	}
	for (ii=ii1; ii<ii2; ++ii){
        fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - ampHy * sS.bE[0] * Fout(timestep, TFSFp[ii].r0H1) / dz;
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + ampHx * sS.bE[0] * Fout(timestep, TFSFp[ii].r0H2) / dz;

		fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - ampEx * sS.bH[0] * Fout(timestep, TFSFp[ii].r0E1) / dz;
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + ampEy * sS.bH[0] * Fout(timestep, TFSFp[ii].r0E2) / dz;
	}
	for (ii=ii2; ii<ii3; ++ii){
        fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - ampHz * sS.bE[0] * Fout(timestep, TFSFp[ii].r0H1) / dy;
		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + ampHx * sS.bE[0] * Fout(timestep, TFSFp[ii].r0H2) / dy;

		fF.Hz[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] - ampEx * sS.bH[0] * Fout(timestep, TFSFp[ii].r0E1) / dy;
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] + ampEz * sS.bH[0] * Fout(timestep, TFSFp[ii].r0E2) / dy;
	}
	for (ii=ii3; ii<ii4; ++ii){
        fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + ampHz * sS.bE[0] * Fout(timestep, TFSFp[ii].r0H1) / dy;
		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - ampHx * sS.bE[0] * Fout(timestep, TFSFp[ii].r0H2) / dy;

		fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + ampEx * sS.bH[0] * Fout(timestep, TFSFp[ii].r0E1) / dy;
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - ampEz * sS.bH[0] * Fout(timestep, TFSFp[ii].r0E2) / dy;
	}
	for (ii=ii4; ii<ii5; ++ii){
        fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + ampHz * sS.bE[0] * Fout(timestep, TFSFp[ii].r0H1) / dx;
		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - ampHy * sS.bE[0] * Fout(timestep, TFSFp[ii].r0H2) / dx;

		fF.Hz[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] + ampEy * sS.bH[0] * Fout(timestep, TFSFp[ii].r0E1) / dx;
		fF.Hy[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] - ampEz * sS.bH[0] * Fout(timestep, TFSFp[ii].r0E2) / dx;
	}
	for (ii=ii5; ii<ii6; ++ii){
        fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - ampHz * sS.bE[0] * Fout(timestep, TFSFp[ii].r0H1) / dx;
		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + ampHy * sS.bE[0] * Fout(timestep, TFSFp[ii].r0H2) / dx;

		fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - ampEy * sS.bH[0] * Fout(timestep, TFSFp[ii].r0E1) / dx;
		fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + ampEz * sS.bH[0] * Fout(timestep, TFSFp[ii].r0E2) / dx;
	}

}

double TFSF::Fout(int n, double r)
{
	nstep = (double)n - r/C/dt;
	nFloor = floor(nstep);
	nprime = nstep - nFloor;
//	if (nstep < 0){
//		return 0;
//	}else{
//		return ((1-nprime)*Ft[nFloor] + nprime*Ft[nFloor+1]);
//	}
	return sin(pinudt*nstep);
}

double GaussianWaveSpatialDistribution(double r2, double C1, double C2, double C3)
{
	return (C1 * sin(C2 * r2) * exp(C3 * r2));
}


//double ContinuousWave(double Amplitude, double retardedtime, int nd, int nw, double pinudt)
//{
//	//retard = ntimestep - location / C / dt;
//	if (retardedtime < nd){
//		return (Amplitude * exp(-(retardedtime - nd)*(retardedtime - nd) / (nw*nw)) * sin(pinudt*retardedtime));
//	}
//	else{
//		return (Amplitude * sin(pinudt*retardedtime));
//	}
//}
//
//double Pulse(double Amplitude, double retardedtime, int nd, int nw, double pinudt)
//{
//	//retard = ntimestep - location / C / dt;
//	return (Amplitude * exp(-double(retardedtime - nd)*(retardedtime - nd) / (nw*nw)) * sin(pinudt*(retardedtime - nd)));
//}
