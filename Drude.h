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
#ifndef DRUDE_H
#define DRUDE_H

#include "Structure.h"
#include "Fields.h"

class Drude
{
public:
	Drude(Structure &sS, char *fileName);
	~Drude();
	double epsInf, omegaP, sigma, gamma;
	double **Jx, **Jy;
	void EFieldUpdate(Field &fF, Structure &sS);

private:
	double beta, kappa, C1Drude, C2Drude, C3Drude;
	
};

#endif

