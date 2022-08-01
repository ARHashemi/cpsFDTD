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
#ifndef DRUDE_LORENTZ_H
#define DRUDE_LORENTZ_H

#include "Structure.h"
#include "Fields.h"

struct DL_material
{
	double *JDx, *JDy, *JLx, *JLy;
	double *JLx_1, *JLy_1, Ex_1, Ey_1;
	int ix, jy;
	double temp1x, temp2x, temp1y, temp2y;
};

class Drude_Lorentz
{
public:
	Drude_Lorentz(Structure &sS, char *fileName);
	~Drude_Lorentz();
	double sigma, epsInf, *omega_D, *gamma_D, *delta_eps_L, *omega_L, *delta_L;
	double **Jx, **Jy;
	int num_of_D_poles, num_of_L_poles, num_of_DLM_gridcells;
	DL_material *DLM;
	void EFieldUpdate(Field &fF, Structure &sS);

private:
	double *beta, *kappa, *alpha, *xi, *gamma_L;
	double C1, C2, C3, *C4, *C5, *C6;
	double denominator, temp1, temp2;

};



#endif
