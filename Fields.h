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
#ifndef FIELDS_H
#define FIELDS_H

#include "Structure.h"
void DynamicMemoryAllocate(double ***D, int sizex, int sizey, int sizez, double initialValue);
void DynamicMemoryDeAllocate(double ***D, int sizex, int sizey);

class Field
{
public:
	Field(Structure &sS, char *smode);
	~Field();
	double ***Ex;
	double ***Ey;
	double ***Ez;
	double ***Hx;
	double ***Hy;
	double ***Hz;
	char *mode;
	char const *TE_mode = "TE";
	char const *TM_mode = "TM";

private:
	int xN;
	int yN;
	int zN;

	//void setSpaceDimensions(int nNx, int nNy);

};



#endif
