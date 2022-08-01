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
#ifndef CPML_ABC_H
#define CPML_ABC_H

#include "Structure.h"
#include "Fields.h"

struct PML_point
{
	int i, j, k, grade;
	double QH1, QH2, QE1, QE2;
};

class CPML_ABC
{
public:
	CPML_ABC(Structure &sS, char *fileName);
	~CPML_ABC();
	void UpdateCPMLEFields(Field &fF);
	void UpdateCPMLHFields(Field &fF);

	PML_point *PMLpoint;

	double *bxE, *cxE, *byE, *cyE, *bzE, *czE, *bxH, *cxH, *byH, *cyH, *bzH, *czH;
	double cpmlE, cpmlH;
	double m, malp, kapmax, alpmax, sigmax;
	int Npmlx, Npmly, Npmlz;
	int Num_of_PML_points;
	int iE1, iE2, iE3;//, iH1, iH2, iH3;

private:
	int xN, yN, zN;
	double *sigmaEx, *kappaEx, *alphaEx;
	double *sigmaEy, *kappaEy, *alphaEy;
	double *sigmaEz, *kappaEz, *alphaEz;
	double *sigmaHx, *kappaHx, *alphaHx;
	double *sigmaHy, *kappaHy, *alphaHy;
	double *sigmaHz, *kappaHz, *alphaHz;
};

#endif // CPML_ABC_H
