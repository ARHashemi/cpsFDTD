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
#include "BasicFDTD.h"
#include "Fields.h"
#include "Structure.h"
#include <omp.h>



void BasicEFieldUpdate(Field &fF, Structure &sS)
{
	short id;
	int i, j, k;
#pragma omp parallel for default(shared) private(i,j,k,id)
	for (i = 1; i < (sS.Nx); ++i){//-1
		for (j = 1; j < (sS.Ny); ++j){//-1
			for (k = 1; k < (sS.Nz); ++k){//-1
				id = sS.ID[i][j][k];
				if (id == 4){
					fF.Ex[i][j][k] = 0.0;
					fF.Ey[i][j][k] = 0.0;
					fF.Ez[i][j][k] = 0.0;
				}else if (id != 2){
					fF.Ex[i][j][k] = sS.aE[id] * fF.Ex[i][j][k] + sS.bE[id] * (fF.Hz[i][j][k] - fF.Hz[i][j - 1][k]) * sS.cExz[j] - sS.bE[id] * (fF.Hy[i][j][k] - fF.Hy[i][j][k - 1]) * sS.cExy[k];
					fF.Ey[i][j][k] = sS.aE[id] * fF.Ey[i][j][k] + sS.bE[id] * (fF.Hx[i][j][k] - fF.Hx[i][j][k - 1]) * sS.cEyx[k] - sS.bE[id] * (fF.Hz[i][j][k] - fF.Hz[i - 1][j][k]) * sS.cEyz[i];
					fF.Ez[i][j][k] = sS.aE[id] * fF.Ez[i][j][k] + sS.bE[id] * (fF.Hy[i][j][k] - fF.Hy[i - 1][j][k]) * sS.cEzy[i] - sS.bE[id] * (fF.Hx[i][j][k] - fF.Hx[i][j - 1][k]) * sS.cEzx[j];
				}
			}
		}
	}
//	id = 0;
//	i = 0;
//	for (k = 1; k<(sS.Nz); ++k){//-1
//		for (j = 1; j<(sS.Ny); ++j){//-1
//			fF.Ex[i][j][k] = sS.aE[id] * fF.Ex[i][j][k] + sS.bE[id] * (fF.Hz[i][j][k] - fF.Hz[i][j - 1][k]) * sS.cExz[j] - sS.bE[id] * (fF.Hy[i][j][k] - fF.Hy[i][j][k - 1]) * sS.cExy[k];
//		}
//	}
//	j = 0;
//	for (k = 1; k<(sS.Nz); ++k){//-1
//		for (i = 1; i<(sS.Nx); ++i){//-1
//			fF.Ey[i][j][k] = sS.aE[id] * fF.Ey[i][j][k] + sS.bE[id] * (fF.Hx[i][j][k] - fF.Hx[i][j][k - 1]) * sS.cEyx[k] - sS.bE[id] * (fF.Hz[i][j][k] - fF.Hz[i - 1][j][k]) * sS.cEyz[i];
//		}
//	}
//	k = 0;
//	for (j = 1; j<(sS.Ny); ++j){//-1
//		for (i = 1; i<(sS.Nx); ++i){//-1
//			fF.Ez[i][j][k] = sS.aE[id] * fF.Ez[i][j][k] + sS.bE[id] * (fF.Hy[i][j][k] - fF.Hy[i - 1][j][k]) * sS.cEzy[i] - sS.bE[id] * (fF.Hx[i][j][k] - fF.Hx[i][j - 1][k]) * sS.cEzx[j];
//		}
//	}
}

void BasicHFieldUpdate(Field &fF, Structure &sS)
{
	short id;
	int i, j, k;
#pragma omp parallel for default(shared) private(i,j,k,id)
	for (i = 0; i < (sS.Nx); ++i){//-1
		for (j = 0; j < (sS.Ny); ++j){//-1
			for (k = 0; k < (sS.Nz); ++k){//-1
				id = sS.ID[i][j][k];
				fF.Hx[i][j][k] = sS.aH[id] * fF.Hx[i][j][k] - sS.bH[id] * (fF.Ez[i][j + 1][k] - fF.Ez[i][j][k]) * sS.cHxz[j] + sS.bH[id] * (fF.Ey[i][j][k + 1] - fF.Ey[i][j][k]) * sS.cHxy[k];
				fF.Hy[i][j][k] = sS.aH[id] * fF.Hy[i][j][k] - sS.bH[id] * (fF.Ex[i][j][k + 1] - fF.Ex[i][j][k]) * sS.cHyx[k] + sS.bH[id] * (fF.Ez[i + 1][j][k] - fF.Ez[i][j][k]) * sS.cHyz[i];
				fF.Hz[i][j][k] = sS.aH[id] * fF.Hz[i][j][k] + sS.bH[id] * (fF.Ex[i][j + 1][k] - fF.Ex[i][j][k]) * sS.cHzx[j] - sS.bH[id] * (fF.Ey[i + 1][j][k] - fF.Ey[i][j][k]) * sS.cHzy[i];
			}
		}
	}
//	id = 0;
//	i = sS.Nx;// - 1
//	for (j = 0; j < (sS.Ny); ++j){//-1
//		for (k = 0; k < (sS.Nz); ++k){//-1
//			fF.Hx[i][j][k] = sS.aH[id] * fF.Hx[i][j][k] - sS.bH[id] * (fF.Ez[i][j + 1][k] - fF.Ez[i][j][k]) * sS.cHxz[j] + sS.bH[id] * (fF.Ey[i][j][k + 1] - fF.Ey[i][j][k]) * sS.cHxy[k];
//		}
//	}
//	j = sS.Ny;// - 1
//	for (i = 0; i < (sS.Nx); ++i){//-1
//		for (k = 0; k < (sS.Nz); ++k){//-1
//			fF.Hy[i][j][k] = sS.aH[id] * fF.Hy[i][j][k] - sS.bH[id] * (fF.Ex[i][j][k + 1] - fF.Ex[i][j][k]) * sS.cHyx[k] + sS.bH[id] * (fF.Ez[i + 1][j][k] - fF.Ez[i][j][k]) * sS.cHyz[i];
//		}
//	}
//	k = sS.Nz;// - 1
//	for (i = 0; i < (sS.Nx); ++i){//-1
//		for (j = 0; j < (sS.Ny); ++j){//-1
//			fF.Hz[i][j][k] = sS.aH[id] * fF.Hz[i][j][k] + sS.bH[id] * (fF.Ex[i][j + 1][k] - fF.Ex[i][j][k]) * sS.cHzx[j] - sS.bH[id] * (fF.Ey[i + 1][j][k] - fF.Ey[i][j][k]) * sS.cHzy[i];
//		}
//	}
}

//void BasicEFieldUpdateTE(Field &fF, Structure &sS)
//{
//	unsigned short id;
//	int i, j, k;
//#pragma omp parallel for default(shared) private(i,j,id)
//	for (i = 1; i < sS.Nx; ++i){
//		for (j = 1; j < sS.Ny; ++j){
//			for (k = 1; k < sS.Nz; ++k){
//			id = sS.ID[i][j][k];
//			fF.Ex[i][j][k] = sS.aE[id] * fF.Ex[i][j][k] + sS.bE[id] * (fF.Hz[i][j][k] - fF.Hz[i][j - 1][k]) * sS.cExy[j];
//			fF.Ey[i][j][k] = sS.aE[id] * fF.Ey[i][j][k] - sS.bE[id] * (fF.Hz[i][j][k] - fF.Hz[i - 1][j][k]) * sS.cEyx[i];
//		}
//	}
//	i = 0;
//	for (j = 1; j<sS.Ny; ++j){
//		fF.Ex[i][j] = sS.aE[0] * fF.Ex[i][j] + sS.bE[0] * (fF.Hz[i][j] - fF.Hz[i][j - 1]) * sS.cExy[j];
//	}
//	j = 0;
//	for (i = 1; i<sS.Nx; ++i){
//		fF.Ey[i][j] = sS.aE[0] * fF.Ey[i][j] - sS.bE[0] * (fF.Hz[i][j] - fF.Hz[i - 1][j]) * sS.cEyx[i];
//	}
//}
//
//
//void BasicHFieldUpdateTE(Field &fF, Structure &sS)
//{
//	unsigned short id;
//	int i, j;
//#pragma omp parallel for default(shared) private(i,j,id)
//	for (i = 0; i < sS.Nx; ++i){
//		for (j = 0; j < sS.Ny; ++j){
//			id = sS.ID[i][j];
//			fF.Hz[i][j] = sS.aH[id] * fF.Hz[i][j] - sS.bH[id] * (fF.Ey[i + 1][j] - fF.Ey[i][j]) * sS.cHzx[i] + sS.bH[id] * (fF.Ex[i][j + 1] - fF.Ex[i][j]) * sS.cHzy[j];
//		}
//	}
//}
