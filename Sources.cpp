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
#include "Constants.h"
#include "Sources.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;



Source::Source(double deltat, char *fileName)
{
	dt = deltat;
	cout << "Constructing Source";
	ifstream InputFile;
	InputFile.open(fileName);
	string sline, id;
	getline(InputFile, sline);
	istringstream iss(sline);
	iss >> lambda;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> FWHM;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> delay;
	nu = C / lambda;
	T = 1.0 / nu;
	nT = T / dt;
	pinudt = 2.0 * pi * nu * dt;
	nw = 2 * nT;
	nd = 4 * nw;
	nTrans = nd + nw;
	cout << " ..." << endl;
}

Source::~Source()
{
}

double Source::ContinuousWave(double Amplitude, double ntimestep, double location)
{
	retard = ntimestep - location / C / dt;
	if (retard<nd){
		return (Amplitude * exp(-(retard - nd)*(retard - nd) / (nw*nw)) * sin(pinudt*retard));
	}
	else{
		return (Amplitude * sin(pinudt*retard));
	}
}

double Source::Pulse(double Amplitude, double ntimestep, double location)
{
	retard = ntimestep - location / C / dt;
	return (Amplitude * exp(-double(retard - nd)*(retard - nd) / (nw*nw)) * sin(pinudt*(retard - nd)));
}
