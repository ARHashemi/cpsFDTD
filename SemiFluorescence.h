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

#ifndef SEMIFLUORESCENCE_H
#define SEMIFLUORESCENCE_H

#include "Structure.h"
#include "Fields.h"



class SemiFluorescence
{
public:
	SemiFluorescence(Structure &sS, char *fileName);
	~SemiFluorescence();
	void UpdateEmitters(Field &fF);
	double lambda, nu, T, EmissionThershold, EmissionAmp;
	int nT, nw, nd;
	int NumSources, i, i1, i2, j1, j2, idif, jdif;
	int **location;
	double *Amps;
	
private:
	double dt, pinudt;
};

#endif
