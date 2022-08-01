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
#ifndef SOURCES_H
#define SOURCES_H

#include "Constants.h"




class Source
{
public:
	Source(double deltat, char *fileName);
	~Source();

	double Pulse(double Amplitude, double ntimestep, double location);
	double lambda, nu, T, delay, FWHM;

	double ContinuousWave(double Amplitude, double ntimestep, double location);

private:
	double dt, retard, pinudt;
	int nd, nw, nT, nTrans;
	double *OF;

};




#endif
