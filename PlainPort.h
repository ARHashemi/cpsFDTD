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
#ifndef PLAINPORT_H
#define PLAINPORT_H

#include "Fields.h"
#include "Structure.h"
#include <string.h>
#include <iostream>
#include <fstream>
using namespace std;


struct Port_Points
{
	int i, j, k;
	double Ex, Ey, Ez;
	double Hx, Hy, Hz;
	unsigned short id;
};

class PlainPort
{
public:
	PlainPort(Structure &sS, const char *fileName);
	~PlainPort();
	void Update_Fields_on_Port(Field &fF);
	void Export_Port(const char *fileName);
	void Export_E_onPort(const char *fileName);
	void Export_E_onPort(ofstream OutFile);
	void Export_Ex_onPort(const char *fileName);
	void Export_Ex_onPort(ofstream OutFile);
	void Export_Ey_onPort(const char *fileName);
	void Export_Ey_onPort(ofstream OutFile);
	void Export_Ez_onPort(const char *fileName);
	void Export_Ez_onPort(ofstream OutFile);
	void Export_EIntensity_onPort(const char *fileName);
	void Export_EIntensity_onPort(ofstream OutFile);
	void Export_ID_onPort(const char *fileName);
	//std::string port_type;
	char normal_direction;
	double position, w1, w2, h1, h2;
	double **port_out;
    Port_Points *port_point;
    int Num_of_Port_Points;


protected:

private:
	int Nx, Ny, Nz;
	int i1, i2, j1, j2, k1, k2, ii, il;
	double dx, dy, dz;
};

#endif // PLAINPORT_H
