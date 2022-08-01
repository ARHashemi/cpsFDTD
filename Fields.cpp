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
#include "Fields.h"
#include "Structure.h"
#include <iostream>
using namespace std;

void DynamicMemoryAllocate(double ***D, int sizex, int sizey, int sizez, double initialValue){
	D = new double**[sizex];
	for (int i = 0; i < sizex; i++){
		D[i] = new double*[sizey];
		for (int j = 0; j < sizey; j++){
			D[i][j] = new double[sizez];
			for (int k = 0; k < sizez; k++){
				D[i][j][k] = initialValue;
			}
		}
	}
}

void DynamicMemoryDeAllocate(double ***D, int sizex, int sizey){
	for (int i = 0; i < sizex; i++){
		for (int j = 0; j < sizey; j++){
			delete[] D[i][j];
		}
		delete[] D[i];
	}
	delete[] D;
	//D = 0;
}

//void Field::setSpaceDimensions(int nNx, int nyN)
//{
//	Nx = nNx;
//	Ny = nNy;
//}

Field::Field(Structure &sS, char *smode)
{
	cout << BOLDRED << ">" << RESET << "Constructing " << BOLDBLUE << "Fields" << RESET << "... ";
	xN = sS.Nx; yN = sS.Ny; zN = sS.Nz;
	mode = smode;
	if (mode == TE_mode)
	{
		Ex = new double**[xN+1];
		Ey = new double**[xN+1];
		Hz = new double**[xN+1];
		for (int i = 0; i <= xN; i++){
			Ex[i] = new double*[yN+1];
			Ey[i] = new double*[yN+1];
			Hz[i] = new double*[yN+1];
			for (int j = 0; j <= yN; j++){
				Ex[i][j] = new double[zN+1];
				Ey[i][j] = new double[zN+1];
				Hz[i][j] = new double[zN+1];
				for (int k = 0; k <= zN; k++){
					Ex[i][j][k] = 0.0;
					Ey[i][j][k] = 0.0;
					Hz[i][j][k] = 0.0;
				}
			}
		}
	}
	else if (mode == TM_mode)
	{
		Hx = new double**[xN+1];
		Hy = new double**[xN+1];
		Ez = new double**[xN+1];
		for (int i = 0; i <= xN; i++){
			Hx[i] = new double*[yN+1];
			Hy[i] = new double*[yN+1];
			Ez[i] = new double*[yN+1];
			for (int j = 0; j <= yN; j++){
				Hx[i][j] = new double[zN+1];
				Hy[i][j] = new double[zN+1];
				Ez[i][j] = new double[zN+1];
				for (int k = 0; k <= zN; k++){
					Hx[i][j][k] = 0.0;
					Hy[i][j][k] = 0.0;
					Ez[i][j][k] = 0.0;
				}
			}
		}
	}
	else
	{
		Ex = new double**[xN+1];
		Ey = new double**[xN+1];
		Hz = new double**[xN+1];
		Hx = new double**[xN+1];
		Hy = new double**[xN+1];
		Ez = new double**[xN+1];
		for (int i = 0; i <= xN; i++){
			Ex[i] = new double*[yN+1];
			Ey[i] = new double*[yN+1];
			Hz[i] = new double*[yN+1];
			Hx[i] = new double*[yN+1];
			Hy[i] = new double*[yN+1];
			Ez[i] = new double*[yN+1];
			for (int j = 0; j <= yN; j++){
				Ex[i][j] = new double[zN+1];
				Ey[i][j] = new double[zN+1];
				Hz[i][j] = new double[zN+1];
				Hx[i][j] = new double[zN+1];
				Hy[i][j] = new double[zN+1];
				Ez[i][j] = new double[zN+1];
				for (int k = 0; k <= zN; k++){
					Ex[i][j][k] = 0.0;
					Ey[i][j][k] = 0.0;
					Hz[i][j][k] = 0.0;
					Hx[i][j][k] = 0.0;
					Hy[i][j][k] = 0.0;
					Ez[i][j][k] = 0.0;
				}
			}
		}
	}

/*	DynamicMemoryAllocate(Ex, xN, yN, 0.0);
	DynamicMemoryAllocate(Ey, xN, yN, 0.0);
	DynamicMemoryAllocate(Ez, xN, yN, 0.0)*/;
	cout << " Done." << endl;
}

Field::~Field()
{
	cout << "Destructing Fields ...";
	if (mode == TE_mode)
	{
		for (int i = 0; i <= xN; i++){
			for (int j = 0; j <= yN; j++){
				delete[] Ex[i][j];
				delete[] Ey[i][j];
				delete[] Hz[i][j];
			}
			delete[] Ex[i];
			delete[] Ey[i];
			delete[] Hz[i];
		}
		delete[] Ex;
		delete[] Ey;
		delete[] Hz;
		Ex = 0;
		Ey = 0;
		Hz = 0;
	}
	else if (mode == TM_mode)
	{
		for (int i = 0; i <= xN; i++){
			for (int j = 0; j <= yN; j++){
				delete[] Hx[i][j];
				delete[] Hy[i][j];
				delete[] Ez[i][j];
			}
			delete[] Hx[i];
			delete[] Hy[i];
			delete[] Ez[i];
		}
		delete[] Hx;
		delete[] Hy;
		delete[] Ez;
		Hx = 0;
		Hy = 0;
		Ez = 0;
	}
	else
	{
		for (int i = 0; i <= xN; i++){
			for (int j = 0; j <= yN; j++){
				delete[] Ex[i][j];
				delete[] Ey[i][j];
				delete[] Hz[i][j];
				delete[] Hx[i][j];
				delete[] Hy[i][j];
				delete[] Ez[i][j];
			}
			delete[] Ex[i];
			delete[] Ey[i];
			delete[] Hz[i];
			delete[] Hx[i];
			delete[] Hy[i];
			delete[] Ez[i];
		}
		delete[] Ex;
		delete[] Ey;
		delete[] Hz;
		delete[] Hx;
		delete[] Hy;
		delete[] Ez;
		Ex = 0;
		Ey = 0;
		Hz = 0;
		Hx = 0;
		Hy = 0;
		Ez = 0;
	}

	/*DynamicMemoryDeAllocate(Ex, xN);
	DynamicMemoryDeAllocate(Ey, xN);
	DynamicMemoryDeAllocate(Ez, xN);*/
	cout << " Done." << endl;
}
