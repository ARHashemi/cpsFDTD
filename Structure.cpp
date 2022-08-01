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
#include "Structure.h"
#include "Constants.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
using namespace std;

Structure::Structure(const char *fileName)
{
    cout << BOLDRED << ">" << RESET << "Constructing " << BOLDBLUE << "Structure" << RESET << "...\n" << BOLDYELLOW << " \u251C" << RESET << "Reading " << BOLDBLUE << "material " << RESET << "data file \"" << fileName << "\"... ";
	ifstream InputFile;
	istringstream iss;
	string sline;
	InputFile.open(fileName);
	getline(InputFile, sline);
	//iss.clear();
	iss.str(sline);
	iss >> NoMedia;
	//cout << NoMedia << endl;
	epsR = new double[NoMedia];
	muR = new double[NoMedia];
	sig = new double[NoMedia];
	sigM = new double[NoMedia];
	for (int i = 0; i < NoMedia; i++){
		getline(InputFile, sline);
		iss.clear();
		iss.str(sline);
		iss >> epsR[i];
		iss >> muR[i];
		iss >> sig[i];
		iss >> sigM[i];
	}
	cout << "Done." << endl;
}

Structure::~Structure()
{
	cout << "Deconstructing Structure...";
	for (int i = 0; i < Nx; i++){
		for (int j = 0; j < Ny; j++){
			delete[] ID[i][j];
		}
		delete[] ID[i];
	}
	delete[] ID;
	ID = 0;
	cout << " Done." << endl;
}

void Structure::generate_mesh(const char *fileName)
{
	cout << BOLDYELLOW << " \u251C" << RESET << "Reading " << BOLDBLUE << "mesh " << RESET << "data file \"" << fileName << "\"... ";
	ifstream InputFile;
	istringstream iss;
	string sline;
	InputFile.open(fileName);
	getline(InputFile, sline);
	iss.str(sline);
	iss >> Lx;
	iss >> Ly;
	iss >> Lz;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> dx;
	iss >> dy;
	iss >> dz;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> freespace;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> NoShapes;
	shapes = new shape[NoShapes];
	for(int i=0; i<NoShapes; i++)
	{
		getline(InputFile, sline);
		iss.clear();
		iss.str(sline);
		iss >> shapes[i].index;
		iss >> shapes[i].material_index;
		iss >> shapes[i].type;
		iss >> shapes[i].theta;
		iss >> shapes[i].phi;
		iss >> shapes[i].psi;
		iss >> shapes[i].x0;
		iss >> shapes[i].y0;
		iss >> shapes[i].z0;
		iss >> shapes[i].a;
		iss >> shapes[i].b;
		iss >> shapes[i].c;

		shapes[i].ix0 = shapes[i].x0/dx;
		shapes[i].jy0 = shapes[i].y0/dy;
		shapes[i].kz0 = shapes[i].z0/dz;
		shapes[i].ia = shapes[i].a/dx;
		shapes[i].jb = shapes[i].b/dy;
		shapes[i].kc = shapes[i].c/dz;
	}
	cout << "Done.\n" << BOLDYELLOW << " \u251C" << RESET << "Initializing " << BOLDBLUE << "mesh" << RESET << "...";
	Nx = Lx/dx;
	Ny = Ly/dy;
	Nz = Lz/dz;
	dt = 0.99 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));
	ID = new short**[Nx+1];
	for (int i = 0; i <= Nx; i++){
		ID[i] = new short*[Ny+1];
		for (int j = 0; j <= Ny; j++){
			ID[i][j] = new short[Nz+1];
			for (int k = 0; k <= Nz; k++){
				ID[i][j][k] = 0;
			}
		}
	}
	if (freespace){
		cout << "Done. " << RED << "Free Space!" << RESET << endl;
	}else{
		cout << "Done.\n" << BOLDYELLOW << " \u251C" << RESET << "Creating " << BOLDBLUE << "objects" << RESET << "...";
		for(int i=0; i<NoShapes; i++)
		{
			if(shapes[i].type=="cuboid")
			{
				cuboid(ID,Nx,Ny,Nz,dx,dy,dz,shapes[i]);
			}else if(shapes[i].type=="ellipsoid")
			{
				ellipsoid(ID,shapes[i],Nx,Ny,Nz,dx,dy,dz);
			}else if(shapes[i].type=="cylinder")
			{
				cylinder(ID,shapes[i],Nx,Ny,Nz,dx,dy,dz);
			}else if(shapes[i].type=="multi")
			{
				multishape(ID,Nx,Ny,Nz,dx,dy,dz,shapes[i],shapes[shapes[i].material_index]);
			}
		}
		cout << "Done.\n";
	}
}

void cuboid(short ***&ID, int Nx, int Ny, int Nz, double dx, double dy, double dz, shape cuboid_shape)
{
	if ((cuboid_shape.theta == 0) && (cuboid_shape.phi == 0) && (cuboid_shape.psi == 0)){
		int i,j,k;
		for (i=0; i<cuboid_shape.ia; ++i){
			for (j=0; j<cuboid_shape.jb; ++j){
				for (k=0; k<cuboid_shape.kc; ++k){
					ID[i+cuboid_shape.ix0][j+cuboid_shape.jy0][k+cuboid_shape.kz0]=cuboid_shape.material_index;
				}
			}
		}
	}else{
		double x1, y1, z1;
		double x2, y2, z2;
		double x3, y3, z3;
		int i3, j3, k3;
		for(double x=0; x<cuboid_shape.a; x=x+(dx/2.0))//
		{
			for(double y=0; y<cuboid_shape.b; y=y+(dy/2.0))///3.0
			{
				for(double z=0; z<cuboid_shape.c; z=z+(dz/2.0))///3.0
				{
					//x-rot
					x1  = x;
					y1 = cos(cuboid_shape.theta)*(y) - sin(cuboid_shape.theta)*(z);
					z1 = sin(cuboid_shape.theta)*(y) + cos(cuboid_shape.theta)*(z);
					//y-rot
					x2 = cos(cuboid_shape.phi)*(x1) + sin(cuboid_shape.phi)*(z1);
					y2 = y1;
					z2 = -sin(cuboid_shape.phi)*(x1) + cos(cuboid_shape.phi)*(z1);
					//z-rot
					x3 = cos(cuboid_shape.psi)*(x2) - sin(cuboid_shape.psi)*(y2);
					y3 = sin(cuboid_shape.psi)*(x2) + cos(cuboid_shape.psi)*(y2);
					z3 = z2;

					i3 = (x3+cuboid_shape.x0)/dx;
					j3 = (y3+cuboid_shape.y0)/dy;
					k3 = (z3+cuboid_shape.z0)/dz;
					if ((i3>=0) && (i3<=Nx) && (j3>=0) && (j3<=Ny) && (k3>=0) && (k3<=Nz))
						ID[i3][j3][k3]=cuboid_shape.material_index;
				}
			}
		}
	}
}

void ellipsoid(short ***&ID, shape ellipsoid_shape, int Nx, int Ny, int Nz, double dx, double dy, double dz)
{
	double temp;
	for(int i=0; i<Nx; i++)
	{
		for(int j=0; j<Ny; j++)
		{
			for(int k=0; k<Nz; k++)
			{
				temp = ((double)i*dx-ellipsoid_shape.x0)*((double)i*dx-ellipsoid_shape.x0)/(ellipsoid_shape.a*ellipsoid_shape.a)+((double)j*dy-ellipsoid_shape.y0)*((double)j*dy-ellipsoid_shape.y0)/(ellipsoid_shape.b*ellipsoid_shape.b)+((double)k*dz-ellipsoid_shape.z0)*((double)k*dz-ellipsoid_shape.z0)/(ellipsoid_shape.c*ellipsoid_shape.c);
				if(temp <= 1.0){
					ID[i][j][k] = ellipsoid_shape.material_index;
				}
			}
		}
	}
}

void cylinder(short ***&ID, shape cylinder_shape, int Nx, int Ny, int Nz, double dx, double dy, double dz)
{
	double temp;
	for(int k=0; k<cylinder_shape.kc; k++){
		for(int i=0; i<Nx; i++){
			for(int j=0; j<Ny; j++){
				temp = ((double)i*dx-cylinder_shape.x0)*((double)i*dx-cylinder_shape.x0)/(cylinder_shape.a*cylinder_shape.a)+((double)j*dy-cylinder_shape.y0)*((double)j*dy-cylinder_shape.y0)/(cylinder_shape.b*cylinder_shape.b);
				if(temp <= 1.0){
					ID[i][j][k+cylinder_shape.kz0] = cylinder_shape.material_index;
				}
			}
		}
	}
}

void multishape(short ***&ID, int Nx, int Ny, int Nz, double dx, double dy, double dz, shape multi_shape, shape single_shape)
{
	shape my_shape;
	my_shape.a = single_shape.a;
	my_shape.b = single_shape.b;
	my_shape.c = single_shape.c;
//	my_shape.ia = single_shape.ia;
	my_shape.index = single_shape.index;
//	my_shape.ix0 = single_shape.ix0;
//	my_shape.jb = single_shape.jb;
//	my_shape.jy0 = single_shape.jy0;
//	my_shape.kc = single_shape.kc;
//	my_shape.kz0 = single_shape.kz0;
	my_shape.material_index = single_shape.material_index;
	my_shape.phi = single_shape.phi;
	my_shape.psi = single_shape.psi;
	my_shape.theta = single_shape.theta;
	my_shape.type = single_shape.type;
	my_shape.ix0 = my_shape.x0/dx;
	my_shape.jy0 = my_shape.y0/dy;
	my_shape.kz0 = my_shape.z0/dz;
	my_shape.ia = my_shape.a/dx;
	my_shape.jb = my_shape.b/dy;
	my_shape.kc = my_shape.c/dz;

	if (single_shape.type == "cuboid"){
		for (int i=0; i<multi_shape.theta; ++i){
			my_shape.x0 = multi_shape.x0 + i * multi_shape.a;
			my_shape.ix0 = my_shape.x0/dx;
			for (int j=0; j<multi_shape.phi; ++j){
				my_shape.y0 = multi_shape.y0 + j * multi_shape.b;
				my_shape.jy0 = my_shape.y0/dy;
				for (int k=0; k<multi_shape.psi; ++k){
					my_shape.z0 = multi_shape.z0 + k * multi_shape.c;
					my_shape.kz0 = my_shape.z0/dz;
					cuboid(ID,Nx,Ny,Nz,dx,dy,dz,my_shape);
				}
			}
		}
	}else if (single_shape.type == "ellipsoid"){
		for (int i=0; i<multi_shape.theta; ++i){
			my_shape.x0 = multi_shape.x0 + i * multi_shape.a;
			my_shape.ix0 = my_shape.x0/dx;
			for (int j=0; j<multi_shape.phi; ++j){
				my_shape.y0 = multi_shape.y0 + j * multi_shape.b;
				my_shape.jy0 = my_shape.y0/dy;
				for (int k=0; k<multi_shape.psi; ++k){
					my_shape.z0 = multi_shape.z0 + k * multi_shape.c;
					my_shape.kz0 = my_shape.z0/dz;
					ellipsoid(ID,my_shape,Nx,Ny,Nz,dx,dy,dz);
				}
			}
		}
	}else if (single_shape.type == "cylinder"){
		for (int i=0; i<multi_shape.theta; ++i){
			my_shape.x0 = multi_shape.x0 + i * multi_shape.a;
			my_shape.ix0 = my_shape.x0/dx;
			for (int j=0; j<multi_shape.phi; ++j){
				my_shape.y0 = multi_shape.y0 + j * multi_shape.b;
				my_shape.jy0 = my_shape.y0/dy;
				for (int k=0; k<multi_shape.psi; ++k){
					my_shape.z0 = multi_shape.z0 + k * multi_shape.c;
					my_shape.kz0 = my_shape.z0/dz;
					cylinder(ID,my_shape,Nx,Ny,Nz,dx,dy,dz);
				}
			}
		}
	}
}


void Structure::InitializeBasicFDTDCoefficients()
{
	cout << BOLDYELLOW << " \u2514" << RESET << "Initializing " << BOLDBLUE << "basic coefficients" << RESET << "...";
	aE = new double[NoMedia];
	bE = new double[NoMedia];
	aH = new double[NoMedia];
	bH = new double[NoMedia];

	cHzy = new double[Nx];
	cHyz = new double[Nx];
	cHzx = new double[Ny];
	cHxz = new double[Ny];
	cHxy = new double[Nz];
	cHyx = new double[Nz];

	cExy = new double[Nz];
	cEyx = new double[Nz];
	cExz = new double[Ny];
	cEzx = new double[Ny];
	cEyz = new double[Nx];
	cEzy = new double[Nx];

	double temp = 1.0 / dx;
	for (int i = 0; i < Nx; i++){
		cEyz[i] = temp;
		cEzy[i] = temp;
		cHyz[i] = temp;
		cHzy[i] = temp;
	}
	temp = 1.0 / dy;
	for (int i = 0; i < Ny; i++){
		cExz[i] = temp;
		cEzx[i] = temp;
		cHxz[i] = temp;
		cHzx[i] = temp;
	}
	temp = 1.0 / dz;
	for (int i = 0; i < Nz; i++){
		cEyx[i] = temp;
		cExy[i] = temp;
		cHyx[i] = temp;
		cHxy[i] = temp;
	}
	for (int i = 0; i<NoMedia; ++i){
		aE[i] = (2.0 * eps0*epsR[i] - sig[i] * dt) / (2.0 * eps0*epsR[i] + sig[i] * dt);
		bE[i] = (2.0 * dt) / (2.0 * eps0*epsR[i] + sig[i] * dt);
		aH[i] = (2.0 * mu0*muR[i] - sigM[i] * dt) / (2.0 * mu0*muR[i] + sigM[i] * dt);
		bH[i] = (2.0 * dt) / (2.0 * mu0*muR[i] + sigM[i] * dt);
	}
	cout << " Done." << endl;
}
