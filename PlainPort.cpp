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
#include "PlainPort.h"
#include "stdafx.h"
#include "Constants.h"
#include "Structure.h"
#include "Fields.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

PlainPort::PlainPort(Structure &sS, const char *fileName)
{
	cout << "Constructing Port...";
	dx = sS.dx;
	dy = sS.dy;
	dz = sS.dz;
	Nx = sS.Nx;
	Ny = sS.Ny;
	Nz = sS.Nz;
	ifstream InputFile;
	InputFile.open(fileName);
	string sline, id;
	getline(InputFile, sline);
	istringstream iss(sline);
	iss >> normal_direction;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> position;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> w1;
	iss >> w2;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> h1;
	iss >> h2;

	ii = 0;
	if (normal_direction == 'x'){
		i1 = position/dx; i2 = i1;
		j1 = w1/dy; j2 = w2/dy;
		k1 = h1/dz; k2 = h2/dz;
		il = k2 - k1 + 1;
		Num_of_Port_Points = (j2-j1+1)*((k2-k1+1) + 1;
		port_point = new Port_Points[Num_of_Port_Points];
		port_out = new double*[j2-j1+1];
		for (int j=j1; j<=j2; j++){
			port_out[j-j1] = new double[k2-k1+1];
			for (int k=k1; k<=k2; k++){
				port_out[j-j1][k-k1] = 0.0;
				port_point[ii].i = i1;
				port_point[ii].j = j;
				port_point[ii].k = k;
				port_point[ii].id = sS.ID[i1][j][k];
				ii = ii + 1;
			}
		}
	} else if (normal_direction == 'y'){
		j1 = position/dy; j2 = j1;
		i1 = w1/dx; i2 = w2/dx;
		k1 = h1/dz; k2 = h2/dz;
		il = k2 - k1 + 1;
		Num_of_Port_Points = (i2-i1+1)*((k2-k1+1) + 1;
		port_point = new Port_Points[Num_of_Port_Points];
		port_out = new double*[i2-i1+1];
		for (int i=i1; i<=i2; i++){
			port_out[i-i1] = new double[k2-k1+1];
			for (int k=k1; k<=k2; k++){
				port_out[i-i1][k-k1] = 0.0;
				port_point[ii].i = i;
				port_point[ii].j = j1;
				port_point[ii].k = k;
				port_point[ii].id = sS.ID[i][j1][k];
				ii = ii + 1;
			}
		}
	} else if (normal_direction == 'z'){
		k1 = position/dz; k2 = k1;
		i1 = w1/dx; i2 = w2/dx;
		j1 = h1/dy; j2 = h2/dy;
		il = j2 - j1 + 1;
		Num_of_Port_Points = (i2-i1+1)*((j2-j1+1) + 1;
		port_point = new Port_Points[Num_of_Port_Points];
		port_out = new double*[i2-i1+1];
		for (int i=i1; i<=i2; i++){
			port_out[i-i1] = new double[j2-j1+1];
			for (int j=j1; j<=j2; j++){
				port_out[i-i1][j-j1] = 0.0;
				port_point[ii].i = i;
				port_point[ii].j = j;
				port_point[ii].k = k1;
				port_point[ii].id = sS.ID[i][j1][k];
				ii = ii + 1;
			}
		}
	}

	cout << " Done." << endl;
}

void Export_ID_onPort(const char *fileName)
{
	ofstream Out(fileName);
	ii = 0;
	for (int i=0; ii<Num_of_Port_Points; ++i){
		for (j=0; j<il; ++j){
			Out << port_point[ii].id << "\t";
			ii = ii+1;
		}
		Out << endl;
    }
    Out.close();
}

void PlainPort::Update_Fields_on_Port(Field &fF)
{
	for (int ii = 0; ii < Num_of_Port_Points; ++ii){
		port_point[ii].Ex = fF.Ex[port_point[ii].i][port_point[ii].j][port_point[ii].k];
		port_point[ii].Ey = fF.Ey[port_point[ii].i][port_point[ii].j][port_point[ii].k];
		port_point[ii].Ez = fF.Ez[port_point[ii].i][port_point[ii].j][port_point[ii].k];

		port_point[ii].Hx = fF.Hx[port_point[ii].i][port_point[ii].j][port_point[ii].k];
		port_point[ii].Hy = fF.Hy[port_point[ii].i][port_point[ii].j][port_point[ii].k];
		port_point[ii].Hz = fF.Hz[port_point[ii].i][port_point[ii].j][port_point[ii].k];
	}
}

void PlainPort::Export_Port(const char *fileName)
{
    ofstream Out(fileName);
    for (int ii=0; ii<Num_of_Port_Points; ++ii){
		Out << port_point[ii].i << "\t" << port_point[ii].j << "\t" << port_point[ii].k << endl;
    }
    Out.close();
}

void PlainPort::Export_E_onPort(const char *fileName)
{
	ofstream Out(fileName);
    for (int ii=0; ii<Num_of_Port_Points; ++ii){
		Out << port_point[ii].Ex << "\t" << port_point[ii].Ey << "\t" << port_point[ii].Ez << endl;
    }
    Out.close();
}
void PlainPort::Export_E_onPort(ofstream OutFile)
{
	for (int ii=0; ii<Num_of_Port_Points; ++ii){
		OutFile << port_point[ii].Ex << "\t" << port_point[ii].Ey << "\t" << port_point[ii].Ez << endl;
    }
    //OutFile << endl;
}

void PlainPort::Export_Ex_onPort(const char *fileName)
{
	ofstream Out(fileName);
	ii = 0;
	for (int i=0; ii<Num_of_Port_Points; ++i){
		for (j=0; j<il; ++j){
			Out << port_point[ii].Ex << "\t";
			ii = ii+1;
		}
		Out << endl;
    }
    Out.close();
}
void PlainPort::Export_Ex_onPort(ofstream OutFile)
{
	ii = 0;
	for (int i=0; ii<Num_of_Port_Points; ++i){
		for (j=0; j<il; ++j){
			OutFile << port_point[ii].Ex << "\t";
			ii = ii+1;
		}
		OutFile << endl;
    }
    //OutFile << endl;
}

void PlainPort::Export_Ey_onPort(const char *fileName)
{
	ofstream Out(fileName);
	ii = 0;
	for (int i=0; ii<Num_of_Port_Points; ++i){
		for (j=0; j<il; ++j){
			Out << port_point[ii].Ey << "\t";
			ii = ii+1;
		}
		Out << endl;
    }
    Out.close();
}
void PlainPort::Export_Ey_onPort(ofstream OutFile)
{
	ii = 0;
    for (int i=0; ii<Num_of_Port_Points; ++i){
		for (j=0; j<il; ++j){
			OutFile << port_point[ii].Ey << "\t";
			ii = ii+1;
		}
		OutFile << endl;
    }
    //OutFile << endl;
}

void PlainPort::Export_Ez_onPort(const char *fileName)
{
	ofstream Out(fileName);
    ii = 0;
    for (int i=0; ii<Num_of_Port_Points; ++i){
		for (j=0; j<il; ++j){
			Out << port_point[ii].Ez << "\t";
			ii = ii+1;
		}
		Out << endl;
    }
    Out.close();
}
void PlainPort::Export_Ez_onPort(ofstream OutFile)
{
	ii = 0;
    for (int i=0; ii<Num_of_Port_Points; ++i){
		for (j=0; j<il; ++j){
			OutFile << port_point[ii].Ez << "\t";
			ii = ii+1;
		}
		OutFile << endl;
    }
    OutFile << endl;
}

void PlainPort::Export_EIntensity_onPort(const char *fileName)
{
	ofstream Out(fileName);
    ii = 0;
    for (int i=0; ii<Num_of_Port_Points; ++i){
		for (j=0; j<il; ++j){
			Out << port_point[ii].Ex*port_point[ii].Ex + port_point[ii].Ey*port_point[ii].Ey + port_point[ii].Ez*port_point[ii].Ez << "\t";
			ii = ii+1;
		}
		Out << endl;
    }
}
void PlainPort::Export_EIntensity_onPort(ofstream OutFile)
{
	ii = 0;
    for (int i=0; ii<Num_of_Port_Points; ++i){
		for (j=0; j<il; ++j){
			Out << port_point[ii].Ex*port_point[ii].Ex + port_point[ii].Ey*port_point[ii].Ey + port_point[ii].Ez*port_point[ii].Ez << "\t";
			ii = ii+1;
		}
		Out << endl;
    }
}

PlainPort::~PlainPort()
{

}
